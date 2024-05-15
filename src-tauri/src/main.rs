#![cfg_attr(
    all(not(debug_assertions), target_os = "windows"),
    windows_subsystem = "windows"
)]

use itertools::Itertools;
use ordered_float::OrderedFloat;
use render::{display_formula, display_mass};
use rustyms::{
    align::{align, matrix::BLOSUM62, Alignment},
    error::*,
    identification::*,
    model::*,
    modification::{
        GnoComposition, LinkerSpecificity, ModificationSearchResult, Ontology, SimpleModification,
    },
    placement_rule::PlacementRule,
    spectrum::*,
    system::{da, dalton, e, mz, usize::Charge, Mass, MassOverCharge},
    ReturnModification, *,
};
use state::State;
use std::sync::Mutex;

use crate::{
    html_builder::{HtmlElement, HtmlTag},
    metadata_render::RenderToHtml, render::display_neutral_loss,
};
use serde::{Deserialize, Serialize};

mod html_builder;
mod metadata_render;
mod render;
mod state;

type ModifiableState<'a> = tauri::State<'a, std::sync::Mutex<State>>;

#[tauri::command]
fn refresh(state: ModifiableState) -> (usize, usize) {
    let state = state.lock().unwrap();
    (state.spectra.len(), state.peptides.len())
}

#[tauri::command]
async fn load_mgf<'a>(path: &'a str, state: ModifiableState<'a>) -> Result<usize, String> {
    match rustyms::rawfile::mgf::open(path) {
        Ok(v) => {
            let count = v.len();
            state.lock().unwrap().spectra = v;
            Ok(count)
        }
        Err(err) => Err(err.to_string()),
    }
}

#[tauri::command]
async fn load_identified_peptides<'a>(
    path: &'a str,
    state: ModifiableState<'a>,
) -> Result<usize, CustomError> {
    let actual_extension = path
        .rsplit('.')
        .next()
        .map(|ex| {
            (ex == "gz")
                .then(|| path.rsplit('.').nth(1))
                .flatten()
                .unwrap_or(ex)
        })
        .map(|ex| ex.to_lowercase());
    match actual_extension.as_deref() {
        Some("csv") => PeaksData::parse_file(path)
            .map(|peptides| {
                state.lock().unwrap().peptides =
                    peptides.filter_map(|p| p.ok()).map(|p| p.into()).collect()
            })
            .or_else(|_| {
                NovorData::parse_file(path).map(|peptides| {
                    state.lock().unwrap().peptides =
                        peptides.filter_map(|p| p.ok()).map(|p| p.into()).collect()
                })
            })
            .map_err(|_| CustomError::error("Unknown file", "Could not be recognised as either a Peaks or Novor file", Context::None)),
        Some("tsv") => SageData::parse_file(path).map(|peptides| {
            state.lock().unwrap().peptides =
                peptides.filter_map(|p| p.ok()).map(|p| p.into()).collect()
        }),
        Some("psmtsv") => OpairData::parse_file(path).map(|peptides| {
            state.lock().unwrap().peptides =
                peptides.filter_map(|p| p.ok()).map(|p| p.into()).collect()
        }),
        Some("fasta") => FastaData::parse_file(path)
            .map(|peptides| {
                state.lock().unwrap().peptides = peptides.into_iter().map(|p| p.into()).collect()
            }),
        _ => Err(CustomError::error("Unknown extension", "Use CSV, TSV, PSMTSV, or Fasta, or any of these as a gzipped file (eg csv.gz).", Context::None)),
    }?;
    Ok(state.lock().unwrap().peptides.len())
}

#[tauri::command]
async fn search_peptide<'a>(
    text: &'a str,
    minimal_match_score: f64,
    minimal_peptide_score: f64,
    state: ModifiableState<'a>,
) -> Result<String, CustomError> {
    let state = state
        .lock()
        .map_err(|_| CustomError::error("Cannot search", "The state is locked, are you trying to do many things at the same time?", Context::None))?;
    let search = LinearPeptide::<Linked>::pro_forma(text, None)?
        .simple()
        .ok_or_else(|| CustomError::error("Invalid search peptide", "A search peptide should be simple", Context::None))?;
    let data = state
        .peptides
        .iter()
        .filter(|p| p.score.map_or(true, |score| score >= minimal_peptide_score))
        .enumerate()
        .map(|(index, peptide)| {
            (
                index,
                align::<4, VerySimple, Simple>(
                    &peptide.peptide,
                    &search,
                    BLOSUM62,
                    Tolerance::new_absolute(da(0.1)),
                    align::AlignType::GLOBAL_B,
                ),
                peptide,
            )
        })
        .sorted_unstable_by(|a, b| b.1.score().normalised.cmp(&a.1.score().normalised))
        .filter(|(_, alignment, _)| alignment.normalised_score() >= minimal_match_score)
        .take(25)
        .map(|(index, alignment, peptide)| {
            let start = alignment.start_a();
            let end = alignment.start_a() + alignment.len_a();
            vec![
                format!("<a onclick=\"document.getElementById('details-identified-peptide-index').value={0};document.getElementById('details-identified-peptide-index').dispatchEvent(new FocusEvent('focus'))\">{0}</a>", index.to_string()),
                format!(
                    "{}<span class='match'>{}</span>{}",
                    peptide.peptide.sub_peptide(..start).to_string(),
                    peptide.peptide.sub_peptide(start..end).to_string(),
                    peptide.peptide.sub_peptide(end..).to_string(),
                ),
                format!("{:.3}", alignment.normalised_score()),
                peptide
                    .score
                    .map(|score| format!("{:.3}", score))
                    .unwrap_or_default(),
            ]
        })
        .collect::<Vec<_>>();
    Ok(html_builder::HtmlElement::table(
        Some(&[
            "Index".to_string(),
            "Sequence".to_string(),
            "Match Score".to_string(),
            "Peptide Score".to_string(),
        ]),
        data,
    )
    .to_string())
}

#[tauri::command]
async fn details_formula(text: &str) -> Result<String, CustomError> {
    let mut formula = if text.is_empty() {
        Err(CustomError::error("Invalid molecular formula", "The test is empty", Context::None))
    } else {
        MolecularFormula::from_pro_forma(text)
    };
    let formula = formula?;
    let isotopes = formula.isotopic_distribution(0.001);
    let (max, max_occurrence) = isotopes
        .iter()
        .enumerate()
        .max_by_key(|f| OrderedFloat(*f.1))
        .unwrap();

    let isotopes_display = if formula.elements().len() == 1
        && formula.elements()[0].1.is_none()
        && formula.elements()[0].2 == 1
    {
        html_builder::HtmlElement::table(
            Some(&["N", "Mass", "Occurrence"]),
            formula.elements()[0]
                .0
                .isotopes()
                .iter()
                .map(|(n, mass, occurrence)| {
                    [
                        n.to_string(),
                        display_mass(*mass).to_string(),
                        if *occurrence == 0.0 {
                            "-".to_string()
                        } else {
                            format!("{:.4}%", occurrence * 100.0)
                        },
                    ]
                }),
        )
        .to_string()
    } else {
        let start = isotopes
            .iter()
            .take_while(|i| **i / *max_occurrence < 0.001)
            .count();
        let end = isotopes
            .iter()
            .rev()
            .take_while(|i| **i / *max_occurrence < 0.001)
            .count();
        let middle = isotopes.len() - start - end;

        format!("<div class='isotopes-distribution'>{}</div>",
            isotopes.iter()
                .copied()
                .enumerate()
                .skip(start)
                .take(middle)
                .map(|(offset, i)| format!("<span class='{} {}' style='--intensity:{}' title='mono + {} Da {:.4}% of total intensity {:.4}% of highest intensity'></span>", 
                    if offset == 0 {"mono"} else {""}, 
                    if offset == max {"most-abundant"} else {""}, 
                    i / *max_occurrence, 
                    offset, 
                    i * 100.0, 
                    i / * max_occurrence * 100.0))
                .join(""))
    };

    Ok(format!(
        "<p>Details on {}</p><p><span style='color:var(--color-red)'>Monoisotopic mass</span> {}, average weight {}, <span style='color:var(--color-green)'>most abundant isotope</span> offset {max} Da</p>{}", 
            display_formula(&formula), 
            display_mass(formula.monoisotopic_mass()), 
            display_mass(formula.average_weight()), 
            isotopes_display,
        ))
}

#[tauri::command]
async fn search_modification(text: &str, tolerance: f64) -> Result<String, CustomError> {
    fn render_places(places: &[PlacementRule]) -> String {
        if places.is_empty() {
            "Anywhere".to_string()
        } else {
        places
            .iter()
            .map(|rule| match rule {
                PlacementRule::AminoAcid(aa, pos) => {
                    format!(
                        "{}@{}",
                        aa.iter().map(|a| a.char()).collect::<String>(),
                        pos
                    )
                }
                PlacementRule::PsiModification(index, pos) => {
                    format!(
                        "{}@{}",
                        Ontology::Psimod.find_id(*index, None).unwrap(),
                        pos
                    )
                }
                PlacementRule::Terminal(pos) => {
                    format!("{}", pos)
                }
            })
            .join(",")
        }
    }

    let modification = if text.is_empty() {
        Err(CustomError::error("Invalid modification", "The modification is empty", Context::None))
    } else {
        Modification::try_from(text, 0..text.len(), &mut Vec::new(), &mut Vec::new(), None)
            .map(|m| match m {
                ReturnModification::Defined(d) => Ok(d),
                _ => Err(
                    CustomError::error("Invalid modification", "Can not define ambiguous modifications for the modifications parameter", Context::None)
                ),
            })
    }??;
    let tolerance = Tolerance::new_absolute(Mass::new::<dalton>(tolerance));
    let result = Modification::search(&modification, tolerance, None);

    match result {
        ModificationSearchResult::Single(modification) => {
            let mut output = HtmlElement::new(HtmlTag::div);

            output = output.content(HtmlElement::new(HtmlTag::p).content(format!(
                "Formula {} monoisotopic mass {} average mass {} most abundant mass {}",
                display_formula(&modification.formula()), 
                display_mass(modification.formula().monoisotopic_mass()),
                display_mass(modification.formula().average_weight()),
                display_mass(modification.formula().most_abundant_mass()),
            )));

            if let Modification::Simple(modification::SimpleModification::Predefined(_, rules, ontology, name, index)) = &modification {
                output = output.content(
                    HtmlElement::new(HtmlTag::p)
                        .content(format!(
                            "Ontology: {ontology}, name: {name}, index: {index}, "
                        ))
                        .maybe_content(
                            modification.ontology_url().map(|url| HtmlElement::new(HtmlTag::a)
                                .content("ontology link")
                                .header("href", url)
                                .header("target", "_blank")),
                        ),
                );
                output = output.content(HtmlElement::new(HtmlTag::p).content("Placement rules"));
                let mut ul = HtmlElement::new(HtmlTag::ul);

                for rule in rules {
                    ul = ul.content(HtmlElement::new(HtmlTag::li).content(
                        format!("Positions: {}{}{}{}{}", render_places(&rule.0),
                            if rule.1.is_empty() { ""} else {", Neutral losses: "},
                            rule.1.iter().map(display_neutral_loss).join(","),
                            if rule.2.is_empty() { ""} else {", Diagnostic ions: "},
                            rule.2.iter().map(|d| &d.0).map(display_formula).join(","))
                        ));                    
                }
                output = output.content(ul);
            } else if let Modification::Simple(SimpleModification::Gno(composition, name)) = &modification {
                output = output.content(
                    HtmlElement::new(HtmlTag::p)
                        .content(format!("Ontology: Gnome, name: {name}, "))
                        .maybe_content(
                            modification.ontology_url().map(|url| HtmlElement::new(HtmlTag::a)
                                .content("ontology link")
                                .header("href", url)
                                .header("target", "_blank")),
                        ),
                );
                match composition {
                    GnoComposition::Mass(_) => {
                        output =
                            output.content(HtmlElement::new(HtmlTag::p).content("Only mass known"));
                    }
                    GnoComposition::Structure(structure) => {
                        output.extend([
                            HtmlElement::new(HtmlTag::p).content(format!(
                                "Composition: {}",
                                structure
                                    .composition()
                                    .iter()
                                    .fold(String::new(), |acc, m| acc + &format!("{}{}", m.0, m.1))
                            )),
                            HtmlElement::new(HtmlTag::p).content(format!("Structure: {structure}")),
                        ]);
                    }
                }
            } else if let Modification::Simple(SimpleModification::Linker { specificities, name, id, length, ontology, diagnostic_ions, .. }) = &modification {
                output = output.content(
                    HtmlElement::new(HtmlTag::p)
                        .content(format!(
                            "Ontology: {ontology}, name: {name}, index: {id}"
                        ))
                        .maybe_content(
                            modification.ontology_url().map(|url| HtmlElement::new(HtmlTag::a)
                                .content("ontology link")
                                .header("href", url)
                                .header("target", "_blank")),
                        ),
                );
                output = output.content(
                    HtmlElement::new(HtmlTag::p)
                        .content(format!(
                            "Length: {}{}{}", length.map_or("-".to_string(), |l|l.to_string()), if diagnostic_ions.is_empty() {""} else {", diagnostic ions: "}, diagnostic_ions.iter().map(|d| &d.0).map(display_formula).join(",")
                        ))
                );
                output = output.content(HtmlElement::new(HtmlTag::p).content("Placement rules"));
                let mut ul = HtmlElement::new(HtmlTag::ul);

                match specificities {
                    LinkerSpecificity::Symmetric(places, stubs) => {
                        ul = ul.content(HtmlElement::new(HtmlTag::li).content(format!("Positions: {}{}{}", render_places(places), if stubs.is_empty() {""} else {", stubs: "}, stubs.iter().map(display_formula).join(","))));
                    }
                    LinkerSpecificity::Asymmetric((left_places, left_stubs),(right_places,right_stubs)) => {
                        ul = ul.content(HtmlElement::new(HtmlTag::li).content(format!("Positions: {}{}{}", render_places(left_places), if left_stubs.is_empty() {""} else {", stubs: "}, left_stubs.iter().map(display_formula).join(","))));
                        ul = ul.content(HtmlElement::new(HtmlTag::li).content(format!("Positions: {}{}{}", render_places(right_places), if right_stubs.is_empty() {""} else {", stubs: "}, right_stubs.iter().map(display_formula).join(","))));
                    }
                }
                output = output.content(ul);
            }

            Ok(output.to_string())
        }
        ModificationSearchResult::Mass(_, _, modifications) => {
            Ok(html_builder::HtmlElement::table(
                Some(&["Name", "Id", "Monoisotopic mass", "Formula"]),
                modifications
                    .iter()
                    .map(|(ontology, id, _, modification)| {
                        [
                            modification.to_string(),
                            format!("<a onclick='document.getElementById(\"search-modification\").value=\"{0}:{1}\";document.getElementById(\"search-modification-button\").click()'>{0}:{1}</a>", ontology.name(), id),
                            display_mass(modification.formula().monoisotopic_mass()).to_string(),
                            display_formula(&modification.formula()),
                        ]
                    }),
            )
            .to_string())
        }
        ModificationSearchResult::Formula(_, modifications) => {
            Ok(html_builder::HtmlElement::table(
                Some(&["Name", "Id"]),
                modifications
                    .iter()
                    .map(|(ontology, id, _, modification)| {
                        [
                            modification.to_string(),
                            format!("<a onclick='document.getElementById(\"search-modification\").value=\"{0}:{1}\";document.getElementById(\"search-modification-button\").click()'>{0}:{1}</a>", ontology.name(), id),
                        ]
                    }),
            )
            .to_string())
        }
        ModificationSearchResult::Glycan(_, modifications) => Ok(html_builder::HtmlElement::table(
            Some(&["Name", "Structure"]),
            modifications
                .iter()
                .map(|(_, _, _, modification)| {
                    [
                        format!("<a onclick='document.getElementById(\"search-modification\").value=\"{0}\";document.getElementById(\"search-modification-button\").click()'>{0}</a>", modification),
                        if let SimpleModification::Gno(GnoComposition::Structure(structure), _) =
                            modification
                        {
                            structure.to_string()
                        } else {
                            "-".to_string()
                        },
                    ]
                }),
        )
        .to_string()),
    }
}

#[derive(Serialize, Deserialize)]
struct Settings {
    peptide: String,
    charge: Option<usize>,
    mode: Option<String>,
    scan_index: Option<usize>,
}

impl Settings {
    fn from_peptide(peptide: &IdentifiedPeptide, scan: Option<usize>) -> Self {
        Self {
            peptide: peptide.peptide.to_string(),
            charge: peptide.metadata.charge().map(|v| v.value),
            mode: peptide
                .metadata
                .mode()
                .map(|mode| {
                    if mode.to_lowercase() == "hcd" || mode.to_lowercase() == "cid" {
                        "CidHcd"
                    } else {
                        mode
                    }
                })
                .map(|mode| mode.to_string()),
            scan_index: scan,
        }
    }
}

#[tauri::command]
fn load_identified_peptide(index: usize, state: ModifiableState) -> Option<Settings> {
    if let Ok(state) = state.lock() {
        state.peptides.get(index).map(|peptide| {
            Settings::from_peptide(
                peptide,
                peptide.metadata.scan_number().and_then(|scan| {
                    state
                        .spectra
                        .iter()
                        .enumerate()
                        .find(|(_, spectrum)| {
                            spectrum
                                .raw_scan_number
                                .map_or(false, |spectrum_scan| scan == spectrum_scan)
                        })
                        .map(|(i, _)| i)
                }),
            )
        })
    } else {
        None
    }
}

#[tauri::command]
fn find_scan_number(scan_number: usize, state: ModifiableState) -> Result<usize, &'static str> {
    state
        .lock()
        .unwrap()
        .spectra
        .iter()
        .position(|scan| scan.raw_scan_number == Some(scan_number))
        .ok_or("Could not find scan number")
}

#[tauri::command]
fn spectrum_details(index: usize, state: ModifiableState) -> String {
    state.lock().unwrap().spectra.get(index).map_or(
        "Spectrum index not valid".to_string(),
        |spectrum| {
            format!(
                "{}\n{:.3}@{:.3}{:+.0}{:.3}{}{}",
                spectrum.title,
                spectrum.mass.value,
                spectrum.rt.value,
                spectrum.charge.value,
                spectrum
                    .intensity
                    .map_or(String::new(), |i| format!(" I:{i}")),
                spectrum
                    .raw_scan_number
                    .map_or(String::new(), |i| format!(" raw scan number:{i}")),
                spectrum
                    .sequence
                    .as_ref()
                    .map(|seq| format!("\n{seq}"))
                    .unwrap_or_default()
            )
        },
    )
}

#[tauri::command]
fn identified_peptide_details(index: usize, state: ModifiableState) -> String {
    state.lock().unwrap().peptides.get(index).map_or(
        "Identified peptide index not valid".to_string(),
        |peptide| peptide.to_html().to_string(),
    )
}

#[tauri::command]
fn load_clipboard(data: &str, state: ModifiableState) -> Result<usize, String> {
    let lines = data.lines().collect_vec();
    if data.is_empty() {
        return Err("Empty clipboard".to_string());
    }
    let spectrum = match lines[0].trim() {
        "#	m/z	Res.	S/N	I	I %	FWHM" => load_bruker_clipboard(&lines),
        "m/z	Charge	Intensity	FragmentType	MassShift	Position" => load_stitch_clipboard(&lines),
        "Mass/Charge Intensity" => load_sciex_clipboard(&lines),
        _ => Err("Not a recognised format (Bruker/Stitch/Sciex)".to_string()),
    }?;

    let mut new_spectrum = RawSpectrum::default();
    new_spectrum.extend(spectrum);
    new_spectrum.title = "Clipboard".to_string();
    new_spectrum.charge = Charge::new::<e>(1);

    state.lock().unwrap().spectra = vec![new_spectrum];
    Ok(1)
}

fn load_bruker_clipboard(lines: &[&str]) -> Result<Vec<RawPeak>, String> {
    let mut spectrum = Vec::new();
    for (line_number, line) in lines.iter().enumerate().skip(1) {
        let line_number = line_number + 1; // Humans like 1 based counting...
        let cells = line.split('\t').collect_vec();
        if cells.len() != 8 {
            return Err(format!("Incorrect number of columns at line {line_number}"));
        }
        if let (Ok(mass_over_charge), Ok(intensity)) = (cells[1].parse(), cells[4].parse()) {
            spectrum.push(RawPeak {
                mz: MassOverCharge::new::<mz>(mass_over_charge),
                intensity,
                charge: Charge::new::<e>(0),
            })
        } else {
            return Err(format!(
                "Could not read numbers at line {line_number} '{line}' '{}' '{}'",
                cells[1], cells[4]
            ));
        }
    }
    Ok(spectrum)
}

fn load_stitch_clipboard(lines: &[&str]) -> Result<Vec<RawPeak>, String> {
    let mut spectrum = Vec::new();
    for (line_number, line) in lines.iter().enumerate().skip(1) {
        let line_number = line_number + 1; // Humans like 1 based counting...
        let cells = line.split('\t').collect_vec();
        if cells.len() != 6 {
            return Err(format!("Incorrect number of columns at line {line_number}"));
        }
        if let (Ok(mass_over_charge), Ok(intensity)) = (cells[0].parse(), cells[2].parse()) {
            spectrum.push(RawPeak {
                mz: MassOverCharge::new::<mz>(mass_over_charge),
                intensity,
                charge: Charge::new::<e>(0),
            })
        } else {
            return Err(format!(
                "Could not read numbers at line {line_number} '{line}' '{}' '{}'",
                cells[0], cells[2]
            ));
        }
    }
    Ok(spectrum)
}

fn load_sciex_clipboard(lines: &[&str]) -> Result<Vec<RawPeak>, String> {
    let mut spectrum = Vec::new();
    for (line_number, line) in lines.iter().enumerate().skip(1) {
        let line_number = line_number + 1; // Humans like 1 based counting...
        let cells = line.split(' ').collect_vec();
        if cells.len() != 2 {
            return Err(format!("Incorrect number of columns at line {line_number}"));
        }
        if let (Ok(mass_over_charge), Ok(intensity)) = (cells[0].parse(), cells[1].parse()) {
            spectrum.push(RawPeak {
                mz: MassOverCharge::new::<mz>(mass_over_charge),
                intensity,
                charge: Charge::new::<e>(0),
            })
        } else {
            return Err(format!(
                "Could not read numbers at line {line_number} '{line}' '{}' '{}'",
                cells[0], cells[1]
            ));
        }
    }
    Ok(spectrum)
}

#[derive(Debug, PartialEq, PartialOrd, Default, Serialize, Deserialize)]
pub enum NoiseFilter {
    #[default]
    None,
    Relative(f64),
    Absolute(f64),
    TopX(f64, usize),
}

#[derive(Debug, PartialEq, PartialOrd, Default, Serialize, Deserialize)]
pub struct ModelParameters {
    pub a: (Location, String),
    pub b: (Location, String),
    pub c: (Location, String),
    pub d: (Location, String),
    pub v: (Location, String),
    pub w: (Location, String),
    pub x: (Location, String),
    pub y: (Location, String),
    pub z: (Location, String),
    pub precursor: String,
    pub immonium: bool,
    pub m: bool,
    pub modification_diagnostic: bool,
    pub modification_neutral: bool,
    pub glycan: (bool, String),
}

#[derive(Debug, PartialEq, PartialOrd, Default, Serialize, Deserialize)]
pub struct AnnotationResult {
    pub spectrum: String,
    pub fragment_table: String,
    pub logs: String,
    pub mz_max: f64,
    pub intensity_max: f64,
}

#[allow(clippy::too_many_arguments)]
#[tauri::command]
async fn annotate_spectrum<'a>(
    index: usize,
    tolerance: (f64, &'a str),
    charge: Option<usize>,
    filter: NoiseFilter,
    model: &'a str,
    peptide: &'a str,
    cmodel: ModelParameters,
    state: ModifiableState<'a>,
    mass_mode: &'a str,
) -> Result<AnnotationResult, CustomError> {
    let state = state.lock().unwrap();
    if index >= state.spectra.len() {
        return Err(CustomError::error(
            "Invalid settings",
            "Non existent spectrum index",
            Context::none(),
        ));
    }
    let get_model_param = |neutral_losses: &String| {
        neutral_losses
            .split(',')
            .filter(|n| !n.is_empty())
            .map(|n| n.parse::<NeutralLoss>())
            .collect::<Result<Vec<_>, _>>()
    };
    let mut model = match model {
        "all" => Model::all(),
        "ethcd" => Model::ethcd(),
        "cidhcd" => Model::cid_hcd(),
        "etd" => Model::etd(),
        "none" => Model::none(),
        "custom" => Model::none()
            .a(cmodel.a.0, get_model_param(&cmodel.a.1)?)
            .b(cmodel.b.0, get_model_param(&cmodel.b.1)?)
            .c(cmodel.c.0, get_model_param(&cmodel.c.1)?)
            .d(cmodel.d.0, get_model_param(&cmodel.d.1)?)
            .v(cmodel.v.0, get_model_param(&cmodel.v.1)?)
            .w(cmodel.w.0, get_model_param(&cmodel.w.1)?)
            .x(cmodel.x.0, get_model_param(&cmodel.x.1)?)
            .y(cmodel.y.0, get_model_param(&cmodel.y.1)?)
            .z(cmodel.z.0, get_model_param(&cmodel.z.1)?)
            .precursor(get_model_param(&cmodel.precursor)?)
            .immonium(cmodel.immonium)
            .m(cmodel.m)
            .modification_specific_diagnostic_ions(cmodel.modification_diagnostic)
            .modification_specific_neutral_losses(cmodel.modification_neutral)
            .glycan(
                cmodel
                    .glycan
                    .0
                    .then(|| get_model_param(&cmodel.glycan.1))
                    .invert()?,
            ),
        _ => Model::all(),
    };
    if tolerance.1 == "ppm" {
        model.tolerance = Tolerance::new_ppm(tolerance.0);
    } else if tolerance.1 == "th" {
        model.tolerance =
            Tolerance::new_absolute(MassOverCharge::new::<rustyms::system::mz>(tolerance.0));
    } else {
        return Err(CustomError::error(
            "Invalid tolerance unit",
            "",
            Context::None,
        ));
    }
    let mass_mode = match mass_mode {
        "monoisotopic" => MassMode::Monoisotopic,
        "average_weight" => MassMode::Average,
        "most_abundant" => MassMode::MostAbundant,
        _ => return Err(CustomError::error("Invalid mass mode", "", Context::None)),
    };
    let peptide = rustyms::CompoundPeptidoform::pro_forma(peptide, None)?;
    let multiple_peptidoforms = peptide.peptidoforms().len() == 1;
    let multiple_peptides = peptide
        .peptidoforms()
        .iter()
        .flat_map(|p| p.peptides())
        .count()
        == 1;
    let mut spectrum = state.spectra[index].clone();
    match filter {
        NoiseFilter::None => (),
        NoiseFilter::Relative(i) => spectrum.relative_noise_filter(i),
        NoiseFilter::Absolute(i) => spectrum.absolute_noise_filter(i),
        NoiseFilter::TopX(size, t) => spectrum.top_x_filter(size, t),
    }
    let use_charge = charge.map_or(spectrum.charge, Charge::new::<e>);
    let fragments = peptide.generate_theoretical_fragments(use_charge, &model);
    let annotated = spectrum.annotate(peptide, &fragments, &model, mass_mode);
    let (spectrum, limits) =
        render::annotated_spectrum(&annotated, "spectrum", &fragments, &model, mass_mode);
    Ok(AnnotationResult {
        spectrum,
        fragment_table: render::spectrum_table(&annotated, &fragments, multiple_peptides),
        logs: format!("{annotated:#?}\n{model:#?}"),
        mz_max: limits.mz.value,
        intensity_max: limits.intensity,
    })
}

fn main() {
    tauri::Builder::default()
        .manage(Mutex::new(State {
            spectra: Vec::new(),
            peptides: Vec::new(),
        }))
        .invoke_handler(tauri::generate_handler![
            annotate_spectrum,
            details_formula,
            find_scan_number,
            identified_peptide_details,
            load_clipboard,
            load_identified_peptide,
            load_identified_peptides,
            load_mgf,
            refresh,
            search_modification,
            search_peptide,
            spectrum_details,
        ])
        .run(tauri::generate_context!())
        .expect("error while running tauri application");
}

pub trait InvertResult<T, E> {
    /// # Errors
    /// If any of the errors contained within has an error.
    fn invert(self) -> Result<Option<T>, E>;
}

impl<T, E> InvertResult<T, E> for Option<Result<T, E>> {
    fn invert(self) -> Result<Option<T>, E> {
        self.map_or_else(|| Ok(None), |o| o.map(|v| Some(v)))
    }
}
impl<T, E> InvertResult<T, E> for Option<Result<Option<T>, E>> {
    fn invert(self) -> Result<Option<T>, E> {
        self.map_or_else(|| Ok(None), |o| o)
    }
}
