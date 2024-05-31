#![cfg_attr(
    all(not(debug_assertions), target_os = "windows"),
    windows_subsystem = "windows"
)]

use itertools::Itertools;
use ordered_float::OrderedFloat;
use render::{display_formula, display_mass};
use rustyms::{
    error::*,
    model::*,
    modification::{ Ontology, SimpleModification,
    },
    spectrum::*,
    system::{e, mz, usize::Charge, MassOverCharge},*,
};
use state::State;
use std::sync::Mutex;

use crate::metadata_render::RenderToHtml;
use serde::{Deserialize, Serialize};

mod html_builder;
mod metadata_render;
mod render;
mod state;
mod search_modification;
mod identified_peptides;

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
async fn details_formula(text: &str) -> Result<String, CustomError> {
    let formula = if text.is_empty() {
        Err(CustomError::error("Invalid molecular formula", "The test is empty", Context::None))
    } else {
        MolecularFormula::from_pro_forma(text)
    }?;
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
fn validate_molecular_formula(text: String) -> Result<String, CustomError> {
    text.parse::<f64>().map(MolecularFormula::with_additional_mass).or_else(|_|MolecularFormula::from_pro_forma(&text)).map(|f| f.hill_notation_html())
}

fn edit_modification(state: ModifiableState, id: usize, name: &str, formula: &str) -> Result<(), CustomError> {
    let formula = formula.parse::<f64>().map(|mass| Ok(MolecularFormula::with_additional_mass(mass))).unwrap_or_else(|_| MolecularFormula::from_pro_forma(formula))?;
    if let Ok(mut state) = state.lock() {
        let modification = (id, name.to_string(), SimpleModification::Database{formula, specificities: Vec::new(), id: modification::ModificationId { ontology: Ontology::Custom, name: name.to_string(), id, description: String::new(), synonyms: Vec::new(), cross_ids: Vec::new() }});
        if let Some(index) = state.database.iter().position(|p| p.0 == id) {
            state.database[index] = modification;
        } else {
            state.database.push(modification);
        }
        Ok(())
    } else {
        Err(CustomError::error("State locked", "Cannot unlock the mutable state, are you doing many things in parallel?", Context::None))
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
    custom_model: ModelParameters,
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
            .a(custom_model.a.0, get_model_param(&custom_model.a.1)?)
            .b(custom_model.b.0, get_model_param(&custom_model.b.1)?)
            .c(custom_model.c.0, get_model_param(&custom_model.c.1)?)
            .d(custom_model.d.0, get_model_param(&custom_model.d.1)?)
            .v(custom_model.v.0, get_model_param(&custom_model.v.1)?)
            .w(custom_model.w.0, get_model_param(&custom_model.w.1)?)
            .x(custom_model.x.0, get_model_param(&custom_model.x.1)?)
            .y(custom_model.y.0, get_model_param(&custom_model.y.1)?)
            .z(custom_model.z.0, get_model_param(&custom_model.z.1)?)
            .precursor(get_model_param(&custom_model.precursor)?)
            .immonium(custom_model.immonium)
            .m(custom_model.m)
            .modification_specific_diagnostic_ions(custom_model.modification_diagnostic)
            .modification_specific_neutral_losses(custom_model.modification_neutral)
            .glycan(
                custom_model
                    .glycan
                    .0
                    .then(|| get_model_param(&custom_model.glycan.1))
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
    let peptide = rustyms::CompoundPeptidoform::pro_forma(peptide, Some(&state.database))?;
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
        fragment_table: render::spectrum_table(&annotated, &fragments, multiple_peptidoforms, multiple_peptides),
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
            database: Vec::new(),
        }))
        .invoke_handler(tauri::generate_handler![
            annotate_spectrum,
            details_formula,
            find_scan_number,
            identified_peptide_details,
            load_clipboard,
            identified_peptides::load_identified_peptide,
            identified_peptides::load_identified_peptides,
            load_mgf,
            refresh,
            search_modification::search_modification,
            identified_peptides::search_peptide,
            spectrum_details,
            render::density_graph,
            validate_molecular_formula,
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
