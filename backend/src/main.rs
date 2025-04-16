#![cfg_attr(
    all(not(debug_assertions), target_os = "windows"),
    windows_subsystem = "windows"
)]

use std::sync::Mutex;

use itertools::Itertools;
use mzdata::prelude::{IonProperties, SpectrumLike};
use ordered_float::OrderedFloat;
use render::{display_formula, display_mass};
use rustyms::{
    error::*,
    modification::SimpleModification,
    spectrum::*,
    system::{e, usize::Charge},
    *,
};
use serde::{Deserialize, Serialize};
use tauri::Manager;

mod custom_modifications;
mod html_builder;
mod identified_peptides;
mod metadata_render;
mod model;
mod render;
mod search_modification;
mod spectra;
mod state;
mod validate;

use crate::{metadata_render::RenderToHtml, state::State};

#[derive(Debug, Copy, Clone, Serialize, Deserialize)]
pub enum Theme {
    Light,
    Dark,
}

impl Theme {
    pub fn fg(self) -> [u8; 3] {
        match self {
            Self::Light => [30, 30, 30],
            Self::Dark => [212, 212, 212],
        }
    }
    pub fn bg(self) -> [u8; 3] {
        match self {
            Self::Light => [255, 255, 255],
            Self::Dark => [30, 30, 30],
        }
    }
}

const CUSTOM_MODIFICATIONS_FILE: &str = "custom_modifications.json";
const CUSTOM_MODELS_FILE: &str = "custom_models.json";
type ModifiableState<'a> = tauri::State<'a, std::sync::Mutex<State>>;

#[tauri::command]
fn refresh(state: ModifiableState) -> (usize, usize) {
    let state = state.lock().unwrap();
    let res = (
        0, // state.spectra.as_ref().map(|s| s.len()).unwrap_or_default(),
        state.identified_peptide_files().len(),
    );
    drop(state);
    res
}

#[tauri::command]
async fn details_formula(text: &str) -> Result<String, CustomError> {
    let formula = if text.is_empty() {
        Err(CustomError::error(
            "Invalid molecular formula",
            "The test is empty",
            Context::None,
        ))
    } else {
        MolecularFormula::from_pro_forma(text, .., false, true, true)
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
                        display_mass(*mass, Some(MassMode::Monoisotopic)).to_string(),
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
        display_formula(&formula, true),
        display_mass(formula.monoisotopic_mass(), Some(MassMode::Monoisotopic)),
        display_mass(formula.average_weight(), Some(MassMode::Average)),
        isotopes_display,
    ))
}

#[tauri::command]
fn identified_peptide_details(
    file: usize,
    index: usize,
    state: ModifiableState,
    theme: Theme,
) -> String {
    state
        .lock()
        .unwrap()
        .identified_peptide_files()
        .iter()
        .find(|f| f.id == file)
        .map_or(
            "Identified peptide file index not valid".to_string(),
            |file| {
                file.peptides.get(index).map_or(
                    "Identified peptide index not valid".to_string(),
                    |peptide| peptide.to_html(theme).to_string(),
                )
            },
        )
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
pub struct AnnotationResult {
    pub spectrum: String,
    pub fragment_table: String,
    pub mz_max: f64,
    pub intensity_max: f64,
}

#[allow(clippy::too_many_arguments)]
#[tauri::command]
async fn annotate_spectrum<'a>(
    tolerance: (f64, &'a str),
    charge: Option<usize>,
    filter: f32,
    model: usize,
    peptide: &'a str,
    state: ModifiableState<'a>,
    mass_mode: &'a str,
    mz_range: (Option<f64>, Option<f64>),
    theme: Theme,
) -> Result<AnnotationResult, CustomError> {
    let mut state = state.lock().unwrap();
    let spectrum = crate::spectra::create_selected_spectrum(&mut state, filter)?;
    let model = crate::model::get_models(&state)
        .1
        .get(model)
        .cloned()
        .ok_or_else(|| CustomError::error("Invalid model", "Model does not exist", Context::None))?
        .2;
    let parameters = model::parameters(tolerance, mz_range)?;
    let mass_mode = match mass_mode {
        "monoisotopic" => MassMode::Monoisotopic,
        "average_weight" => MassMode::Average,
        "most_abundant" => MassMode::MostAbundant,
        _ => return Err(CustomError::error("Invalid mass mode", "", Context::None)),
    };
    let peptide = rustyms::CompoundPeptidoformIon::pro_forma(peptide, Some(&state.database))?;
    let multiple_peptidoforms = peptide.peptidoform_ions().len() == 1;
    let multiple_peptides = peptide
        .peptidoform_ions()
        .iter()
        .flat_map(|p| p.peptidoforms())
        .count()
        == 1;

    let use_charge = Charge::new::<e>(
        charge
            .or_else(|| {
                spectrum.precursor().and_then(|p| {
                    p.charge().map(|c| {
                        usize::try_from(c).expect("Can not handle negative mode mass spectrometry")
                    })
                })
            })
            .unwrap_or(1),
    );
    let fragments = peptide.generate_theoretical_fragments(use_charge, model);
    let annotated = spectrum.annotate(peptide, &fragments, &parameters, mass_mode);
    let (spectrum, limits) = render::annotated_spectrum(
        &annotated,
        &spectrum,
        "spectrum",
        &fragments,
        model,
        &parameters,
        mass_mode,
        theme,
    );
    Ok(AnnotationResult {
        spectrum,
        fragment_table: render::spectrum_table(
            &annotated,
            &fragments,
            multiple_peptidoforms,
            multiple_peptides,
        ),
        mz_max: limits.mz.value,
        intensity_max: limits.intensity,
    })
}

fn load_custom_mods_and_models(app: &mut tauri::App) -> Result<(), Box<dyn std::error::Error>> {
    let path = app.path().app_config_dir();
    if let Ok(path) = path {
        let custom_mods = path.join(CUSTOM_MODIFICATIONS_FILE);
        let custom_models = path.join(CUSTOM_MODELS_FILE);
        let app_state = app.state::<Mutex<State>>();
        let mut state = app_state
            .lock()
            .expect("Poisoned mutex at setup of custom mods");
        if let Ok(data) = std::fs::read(custom_mods) {
            let mods: Vec<(Option<usize>, String, SimpleModification)> =
                serde_json::from_slice(&data)?;
            state.database = mods;
        }
        if let Ok(data) = std::fs::read(custom_models) {
            let models: Vec<(String, FragmentationModel)> = serde_json::from_slice(&data)?;
            state.models = models;
        }
        Ok(())
    } else {
        Err("Could not find configuration file".into())
    }
}

#[tauri::command]
fn get_custom_configuration_path(app: tauri::AppHandle) -> (String, String) {
    app.path().app_config_dir().map_or(
        ("not loaded".to_string(), "not loaded".to_string()),
        |dir| {
            (
                dir.join(CUSTOM_MODIFICATIONS_FILE)
                    .to_string_lossy()
                    .to_string(),
                dir.join(CUSTOM_MODELS_FILE).to_string_lossy().to_string(),
            )
        },
    )
}

fn main() {
    tauri::Builder::default()
        .plugin(tauri_plugin_dialog::init())
        .manage(Mutex::new(State {
            spectra: Vec::new(),
            identified_peptide_files: std::cell::RefCell::new(Vec::new()),
            database: Vec::new(),
            models: Vec::new(),
        }))
        .setup(load_custom_mods_and_models)
        .invoke_handler(tauri::generate_handler![
            annotate_spectrum,
            custom_modifications::delete_custom_modification,
            custom_modifications::duplicate_custom_modification,
            custom_modifications::get_custom_modification,
            custom_modifications::get_custom_modifications,
            custom_modifications::update_modification,
            details_formula,
            get_custom_configuration_path,
            identified_peptide_details,
            identified_peptides::close_identified_peptides_file,
            identified_peptides::get_identified_peptides_files,
            identified_peptides::load_identified_peptide,
            identified_peptides::load_identified_peptides_file,
            identified_peptides::search_peptide,
            model::delete_custom_model,
            model::duplicate_custom_model,
            model::get_custom_model,
            model::get_custom_models,
            model::update_model,
            refresh,
            render::density_graph,
            search_modification::search_modification,
            spectra::close_raw_file,
            spectra::deselect_spectrum,
            spectra::get_open_raw_files,
            spectra::get_selected_spectra,
            spectra::load_clipboard,
            spectra::load_raw,
            spectra::load_usi,
            spectra::save_spectrum,
            spectra::select_spectrum_index,
            spectra::select_spectrum_native_id,
            validate::validate_aa_neutral_loss,
            validate::validate_amino_acid,
            validate::validate_custom_linker_specificity,
            validate::validate_custom_single_specificity,
            validate::validate_fragment_kind,
            validate::validate_glycan_fragments,
            validate::validate_molecular_formula,
            validate::validate_monosaccharide_neutral_loss,
            validate::validate_neutral_loss,
            validate::validate_placement_rule,
            validate::validate_satellite_ion,
            validate::validate_stub,
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
