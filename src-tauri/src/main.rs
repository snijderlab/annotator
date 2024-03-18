#![cfg_attr(
    all(not(debug_assertions), target_os = "windows"),
    windows_subsystem = "windows"
)]

use itertools::Itertools;
use render::display_mass;
use rustyms::{
    align::{align, matrix::BLOSUM62, Alignment},
    error::*,
    identification::*,
    model::*,
    modification::{GnoComposition, ModificationSearchResult, Ontology, ReturnModification},
    placement_rule::PlacementRule,
    spectrum::*,
    system::*,
    *,
};
use state::State;
use std::sync::Mutex;

use crate::{
    html_builder::{HtmlElement, HtmlTag},
    metadata_render::RenderToHtml,
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
) -> Result<usize, String> {
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
            .map_err(|_| "Could not be recognised as either a Peaks or Novor file".to_string()),
        Some("psmtsv") => OpairData::parse_file(path).map(|peptides| {
            state.lock().unwrap().peptides =
                peptides.filter_map(|p| p.ok()).map(|p| p.into()).collect()
        }),
        Some("fasta") => FastaData::parse_file(path)
            .map(|peptides| {
                state.lock().unwrap().peptides = peptides.into_iter().map(|p| p.into()).collect()
            })
            .map_err(|err| err.to_string()),
        _ => Err("Not a recognised extension".to_string()),
    }?;
    Ok(state.lock().unwrap().peptides.len())
}

#[tauri::command]
async fn search_peptide<'a>(
    text: &'a str,
    minimal_match_score: f64,
    minimal_peptide_score: f64,
    state: ModifiableState<'a>,
) -> Result<String, String> {
    let state = state
        .lock()
        .map_err(|_| "Search exception: State locked".to_string())?;
    let search = LinearPeptide::pro_forma(text).map_err(|err| err.to_string())?;
    let data = state
        .peptides
        .iter()
        .filter(|p| p.score.map_or(true, |score| score >= minimal_peptide_score))
        .enumerate()
        .map(|(index, peptide)| {
            (
                index,
                align::<4>(
                    &peptide.peptide,
                    &search,
                    BLOSUM62,
                    Tolerance::new_absolute(da(0.1)),
                    align::AlignType::GLOBAL_B,
                ),
                peptide,
            )
        })
        .sorted_unstable_by(|a, b| b.1.cmp(&a.1))
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
        &data,
    )
    .to_string())
}

#[tauri::command]
async fn search_modification(text: &str, tolerance: f64) -> Result<String, String> {
    let modification = if text.is_empty() {
        Err("Empty".to_string())
    } else {
        Modification::try_from(text, 0..text.len(), &mut Vec::new())
            .map(|m| match m {
                ReturnModification::Defined(d) => Ok(d),
                _ => Err(
                    "Can not define ambiguous modifications for the modifications parameter"
                        .to_string(),
                ),
            })
            .map_err(|err| err.to_string())
    }??;
    let tolerance = Tolerance::new_absolute(Mass::new::<dalton>(tolerance));
    let result = Modification::search(&modification, tolerance);

    match result {
        ModificationSearchResult::Single(modification) => {
            let mut output = HtmlElement::new(HtmlTag::div);

            output = output.content(HtmlElement::new(HtmlTag::p).content(format!(
                "Formula <span class='formula'>{}</span> monoisotopic mass {} average mass {}",
                modification.formula().hill_notation_html(),
                display_mass(modification.formula().monoisotopic_mass()),
                display_mass(modification.formula().average_weight()),
            )));

            if let Modification::Predefined(_, rules, ontology, name, index) = &modification {
                output = output.content(
                    HtmlElement::new(HtmlTag::p)
                        .content(format!(
                            "Ontology: {ontology}, name: {name}, index: {index}, "
                        ))
                        .content(
                            HtmlElement::new(HtmlTag::a)
                                .content("ontology link")
                                .header("href", modification.ontology_url().unwrap_or_default())
                                .header("target", "_blank"),
                        ),
                );
                output = output.content(HtmlElement::new(HtmlTag::p).content("Placement rules"));
                let mut ul = HtmlElement::new(HtmlTag::ul);

                for rule in rules {
                    match rule {
                        PlacementRule::AminoAcid(aa, pos) => {
                            ul = ul.content(HtmlElement::new(HtmlTag::li).content(format!(
                                "{}@{}",
                                aa.iter().map(|a| a.char()).collect::<String>(),
                                pos
                            )));
                        }
                        PlacementRule::PsiModification(index, pos) => {
                            ul = ul.content(HtmlElement::new(HtmlTag::li).content(format!(
                                "{}@{}",
                                Ontology::Psimod.find_id(*index).unwrap(),
                                pos
                            )));
                        }
                        PlacementRule::Terminal(pos) => {
                            ul = ul
                                .content(HtmlElement::new(HtmlTag::li).content(format!("{}", pos)));
                        }
                    }
                }
                output = output.content(ul);
            } else if let Modification::Gno(composition, name) = &modification {
                output = output.content(
                    HtmlElement::new(HtmlTag::p)
                        .content(format!("Ontology: Gnome, name: {name}, "))
                        .content(
                            HtmlElement::new(HtmlTag::a)
                                .content("ontology link")
                                .header("href", modification.ontology_url().unwrap_or_default())
                                .header("target", "_blank"),
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
            }

            Ok(output.to_string())
        }
        ModificationSearchResult::Mass(_, _, modifications) => {
            Ok(html_builder::HtmlElement::table(
                Some(&["Name", "Id", "Monoisotopic mass", "Formula"]),
                &modifications
                    .iter()
                    .map(|(ontology, id, _, modification)| {
                        [
                            modification.to_string(),
                            format!("{}:{}", ontology.name(), id),
                            display_mass(modification.formula().monoisotopic_mass()).to_string(),
                            format!(
                                "<span class='formula'>{}</span>",
                                modification.formula().hill_notation_html()
                            ),
                        ]
                    })
                    .collect_vec(),
            )
            .to_string())
        }
        ModificationSearchResult::Formula(_, modifications) => {
            Ok(html_builder::HtmlElement::table(
                Some(&["Name", "Id"]),
                &modifications
                    .iter()
                    .map(|(ontology, id, _, modification)| {
                        [
                            modification.to_string(),
                            format!("{}:{}", ontology.name(), id),
                        ]
                    })
                    .collect_vec(),
            )
            .to_string())
        }
        ModificationSearchResult::Glycan(_, modifications) => Ok(html_builder::HtmlElement::table(
            Some(&["Name", "Id", "Structure"]),
            &modifications
                .iter()
                .map(|(ontology, id, _, modification)| {
                    [
                        modification.to_string(),
                        format!("{}:{}", ontology.name(), id),
                        if let Modification::Gno(GnoComposition::Structure(structure), _) =
                            modification
                        {
                            structure.to_string()
                        } else {
                            "-".to_string()
                        },
                    ]
                })
                .collect_vec(),
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
            charge: peptide.metadata.charge().map(|v| v.value as usize),
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
    // TODO: show local confidence on the sequence (maybe as done in stitch before?)
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
    new_spectrum.charge = Charge::new::<e>(1.0);

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
                charge: Charge::new::<e>(0.0),
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
                charge: Charge::new::<e>(0.0),
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
                charge: Charge::new::<e>(0.0),
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
    ppm: f64,
    charge: Option<f64>,
    filter: NoiseFilter,
    model: &'a str,
    peptide: &'a str,
    cmodel: ModelParameters,
    state: ModifiableState<'a>,
) -> Result<AnnotationResult, CustomError> {
    let state = state.lock().unwrap();
    if index >= state.spectra.len() {
        return Err(CustomError::error(
            "Invalid settings",
            "Non existent spectrum index",
            Context::none(),
        ));
    }
    let get_model_param = |(location, neutral_losses): &(Location, String)| match neutral_losses
        .split(',')
        .filter(|n| !n.is_empty())
        .map(|n| n.parse())
        .collect()
    {
        Ok(n) => Ok((location.clone(), n)),
        Err(err) => Err(err),
    };
    let mut model = match model {
        "all" => Model::all(),
        "ethcd" => Model::ethcd(),
        "cidhcd" => Model::cid_hcd(),
        "etd" => Model::etd(),
        "none" => Model::none(),
        "custom" => Model::new(
            get_model_param(&cmodel.a)?,
            get_model_param(&cmodel.b)?,
            get_model_param(&cmodel.c)?,
            get_model_param(&cmodel.d)?,
            get_model_param(&cmodel.v)?,
            get_model_param(&cmodel.w)?,
            get_model_param(&cmodel.x)?,
            get_model_param(&cmodel.y)?,
            get_model_param(&cmodel.z)?,
            cmodel
                .precursor
                .to_owned()
                .split(',')
                .filter(|n| !n.is_empty())
                .map(|n| n.parse::<NeutralLoss>())
                .collect::<Result<Vec<NeutralLoss>, _>>()?,
            MassOverCharge::new::<mz>(ppm),
            cmodel
                .glycan
                .0
                .then(|| {
                    cmodel
                        .glycan
                        .1
                        .to_owned()
                        .split(',')
                        .filter(|n| !n.is_empty())
                        .map(|n| n.parse::<NeutralLoss>())
                        .collect::<Result<Vec<NeutralLoss>, _>>()
                })
                .invert()?,
        ),
        _ => Model::all(),
    };
    model.ppm = MassOverCharge::new::<mz>(ppm);
    let peptide = rustyms::ComplexPeptide::pro_forma(peptide)?;
    let multiple_peptides = peptide.peptides().len() != 1;
    let mut spectrum = state.spectra[index].clone();
    match filter {
        NoiseFilter::None => (),
        NoiseFilter::Relative(i) => spectrum.relative_noise_filter(i),
        NoiseFilter::Absolute(i) => spectrum.absolute_noise_filter(i),
        NoiseFilter::TopX(size, t) => spectrum.top_x_filter(size, t),
    }
    let use_charge = charge.map_or(spectrum.charge, Charge::new::<e>);
    let fragments = peptide.generate_theoretical_fragments(use_charge, &model);
    let annotated = spectrum.annotate(peptide, &fragments, &model, MassMode::Monoisotopic);
    let (spectrum, limits) = render::annotated_spectrum(&annotated, "spectrum", &fragments, &model);
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
            refresh,
            load_mgf,
            load_identified_peptides,
            load_clipboard,
            find_scan_number,
            spectrum_details,
            search_peptide,
            search_modification,
            identified_peptide_details,
            load_identified_peptide,
            annotate_spectrum,
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
