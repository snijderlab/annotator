#![cfg_attr(
    all(not(debug_assertions), target_os = "windows"),
    windows_subsystem = "windows"
)]

use std::sync::Mutex;

use context_error::{BasicKind, BoxedError, CreateError, FullErrorContent};
use itertools::Itertools;
use mzannotate::{mzspeclib::AnalyteTarget, prelude::*};
use mzcore::{
    prelude::*,
    system::{e, isize::Charge},
};
use mzdata::prelude::{IonProperties, SpectrumLike};
use mzident::MetaData;
use ordered_float::OrderedFloat;
use render::{display_formula, display_mass};
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

#[derive(Clone, Copy, Debug, Deserialize, Serialize)]
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
fn refresh(state: ModifiableState, theme: Theme) -> (usize, usize, Option<AnnotationResult>) {
    let state = state.lock().unwrap();
    (
        0, // state.spectra.as_ref().map(|s| s.len()).unwrap_or_default(),
        state.identified_peptide_files().len(),
        state.annotated_spectrum.as_ref().map(|annotated| {
            render_annotated_spectrum(
                annotated,
                &[],
                FragmentationModel::all(),
                &MatchingParameters::default(),
                MassMode::Monoisotopic,
                theme,
            )
        }),
    )
}

#[tauri::command]
async fn details_formula(text: &str) -> Result<String, String> {
    let formula = if text.is_empty() {
        Err(BoxedError::small(
            BasicKind::Error,
            "Invalid molecular formula",
            "The input is empty",
        )
        .to_html(false))
    } else {
        MolecularFormula::from_pro_forma::<false, false>(text, ..).map_err(|err| err.to_html(false))
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
) -> (bool, String) {
    state
        .lock()
        .unwrap()
        .identified_peptide_files()
        .iter()
        .find(|f| f.id == file)
        .map_or(
            (false, "Identified peptide file index not valid".to_string()),
            |file| {
                file.peptides.get(index).map_or(
                    (false, "Identified peptide index not valid".to_string()),
                    |peptide| {
                        (
                            peptide.has_annotated_spectrum(),
                            peptide.to_html(theme).to_string(),
                        )
                    },
                )
            },
        )
}

#[tauri::command]
fn load_annotated_spectrum(
    file: usize,
    index: usize,
    state: ModifiableState,
    theme: Theme,
) -> Result<AnnotationResult, String> {
    let mut state = state.lock().unwrap();
    let annotated = state
        .identified_peptide_files()
        .iter()
        .find(|f| f.id == file)
        .map_or(
            Err("Identified peptide file index not valid".to_string()),
            |file| {
                file.peptides.get(index).map_or(
                    Err("Identified peptide index not valid".to_string()),
                    |peptide| {
                        peptide.annotated_spectrum().map_or(
                            Err(
                                "Identified peptide does not have an associated annotated spectrum"
                                    .to_string(),
                            ),
                            |a| Ok(a.into_owned()),
                        )
                    },
                )
            },
        )?;
    let rendered = render_annotated_spectrum(
        &annotated,
        &[],
        FragmentationModel::all(),
        &MatchingParameters::default(),
        MassMode::Monoisotopic,
        theme,
    );
    state.annotated_spectrum = Some(annotated);
    Ok(rendered)
}

#[derive(Debug, Default, Deserialize, PartialEq, PartialOrd, Serialize)]
pub enum NoiseFilter {
    #[default]
    None,
    Relative(f64),
    Absolute(f64),
    TopX(f64, usize),
}

#[derive(Debug, Default, Deserialize, PartialEq, PartialOrd, Serialize)]
pub struct AnnotationResult {
    pub spectrum: String,
    pub fragment_table: String,
    pub mz_max: f64,
    pub intensity_max: f32,
}

#[allow(clippy::too_many_arguments)]
#[tauri::command]
async fn annotate_spectrum<'a>(
    tolerance: (f64, &'a str),
    charge: Option<isize>,
    filter: f32,
    model: usize,
    peptide: &'a str,
    state: ModifiableState<'a>,
    mass_mode: &'a str,
    mz_range: (Option<f64>, Option<f64>),
    theme: Theme,
) -> Result<(AnnotationResult, Vec<String>), Vec<String>> {
    let mut state = state.lock().unwrap();
    let spectrum = crate::spectra::create_selected_spectrum(&mut state, filter)
        .map_err(|err| vec![err.to_html(false)])?;
    let model = crate::model::get_models(&state)
        .1
        .get(model)
        .cloned()
        .ok_or_else(|| {
            vec![
                BoxedError::small(BasicKind::Error, "Invalid model", "Model does not exist")
                    .to_html(false),
            ]
        })?
        .2;
    let parameters =
        model::parameters(tolerance, mz_range).map_err(|err| vec![err.to_html(false)])?;
    let mass_mode = match mass_mode {
        "monoisotopic" => MassMode::Monoisotopic,
        "average_weight" => MassMode::Average,
        "most_abundant" => MassMode::MostAbundant,
        _ => {
            return Err(vec![
                BoxedError::small(BasicKind::Error, "Invalid mass mode", "").to_html(false),
            ]);
        }
    };
    let (peptide, warnings) =
        CompoundPeptidoformIon::pro_forma(peptide, Some(&state.custom_modifications)).map_err(
            |errs| {
                errs.into_iter()
                    .map(|err| err.to_html(false))
                    .collect::<Vec<_>>()
            },
        )?;

    let use_charge = Charge::new::<e>(
        charge
            .or_else(|| {
                spectrum
                    .precursor()
                    .and_then(|p| p.charge().map(|c| c as isize))
            })
            .unwrap_or(1),
    );
    let fragments = peptide.generate_theoretical_fragments(use_charge, model);
    let annotated = spectrum.annotate(peptide, &fragments, &parameters, mass_mode);
    let rendered =
        render_annotated_spectrum(&annotated, &fragments, model, &parameters, mass_mode, theme);
    state.annotated_spectrum = Some(annotated);
    Ok((
        rendered,
        warnings
            .into_iter()
            .map(|err| err.to_html(false))
            .collect::<Vec<_>>(),
    ))
}

fn render_annotated_spectrum(
    annotated: &AnnotatedSpectrum,
    fragments: &[Fragment],
    model: &FragmentationModel,
    parameters: &MatchingParameters,
    mass_mode: MassMode,
    theme: Theme,
) -> AnnotationResult {
    let multiple_peptidoforms = annotated
        .analytes
        .iter()
        .filter(|a| matches!(&a.target, AnalyteTarget::PeptidoformIon(_)))
        .count()
        == 1;
    let multiple_peptides = annotated
        .analytes
        .iter()
        .filter_map(|a| match &a.target {
            AnalyteTarget::PeptidoformIon(pep) => Some(pep),
            _ => None,
        })
        .flat_map(|p| p.peptidoforms())
        .count()
        == 1;
    let (spectrum, limits) = render::annotated_spectrum(
        annotated, "spectrum", fragments, model, parameters, mass_mode, theme,
    );
    AnnotationResult {
        spectrum,
        fragment_table: render::spectrum_table(
            annotated,
            fragments,
            multiple_peptidoforms,
            multiple_peptides,
        ),
        mz_max: limits.mz.value,
        intensity_max: limits.intensity,
    }
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
        let now = chrono::Local::now().format("%Y%m%d%H%M%S%.3f");
        if custom_mods.exists() {
            match mzcore::sequence::parse_custom_modifications(&custom_mods) {
                Ok(modifications) => state.custom_modifications = modifications,
                Err(error) => {
                    eprintln!("Error while parsing custom modifications:\n{error}");
                    let mut combined_error = (error.to_html(false), Vec::new());
                    if let Err(err) = std::fs::rename(
                        &custom_mods,
                        custom_mods
                            .with_file_name(format!("{CUSTOM_MODELS_FILE}_backup_{now}.json",)),
                    ) {
                        combined_error
                            .1
                            .push(format!("Could not rename modifications file: {err}"));
                    }
                    if let Err(err) = std::fs::write(
                        custom_mods.with_file_name(format!("error_{now}.txt")),
                        error.to_string().as_bytes(),
                    ) {
                        combined_error.1.push(format!(
                            "Could not save modifications parse error file: {err}"
                        ));
                    }
                    state.custom_modifications_error = Some(combined_error);
                }
            }
        }
        if custom_models.exists() {
            match mzannotate::annotation::model::parse_custom_models(&custom_models) {
                Ok(models) => state.custom_models = models,
                Err(error) => {
                    eprintln!("Error while parsing custom models:\n{error}");
                    let mut combined_error = (error.to_html(false), Vec::new());
                    if let Err(err) = std::fs::rename(
                        &custom_models,
                        custom_models
                            .with_file_name(format!("{CUSTOM_MODELS_FILE}_backup_{now}.json")),
                    ) {
                        combined_error
                            .1
                            .push(format!("Could not rename models file: {err}"));
                    }
                    if let Err(err) = std::fs::write(
                        custom_models.with_file_name(format!("error_{now}.txt")),
                        error.to_string().as_bytes(),
                    ) {
                        combined_error
                            .1
                            .push(format!("Could not save models parse error file: {err}"));
                    }
                    state.custom_models_error = Some(combined_error);
                }
            }
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
        .plugin(tauri_plugin_opener::init())
        .manage(Mutex::new(State {
            spectra: Vec::new(),
            identified_peptide_files: std::cell::RefCell::new(Vec::new()),
            annotated_spectrum: None,
            custom_modifications: Vec::new(),
            custom_modifications_error: None,
            custom_models: Vec::new(),
            custom_models_error: None,
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
            load_annotated_spectrum,
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
            spectra::select_retention_time,
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
