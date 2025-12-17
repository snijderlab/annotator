#![cfg_attr(
    all(not(debug_assertions), target_os = "windows"),
    windows_subsystem = "windows"
)]

use std::path::PathBuf;

use tokio::sync::Mutex;

use clap::Parser;
use context_error::{BasicKind, BoxedError, CreateError, FullErrorContent};
use itertools::Itertools;
use mzannotate::{mzspeclib::AnalyteTarget, prelude::*};
use mzcore::{
    ontology::Ontologies,
    prelude::*,
    system::{e, isize::Charge},
};
use mzcv::{CVIndex, CVSource};
use mzdata::prelude::SpectrumLike;
use mzident::PSMMetaData;
use ordered_float::OrderedFloat;
use render::{display_formula, display_mass};
use serde::{Deserialize, Serialize};
use tauri::Manager;

mod custom_modifications;
mod html_builder;
mod identified_peptides;
mod metadata_render;
mod model;
mod psm_file;
mod raw_file;
mod render;
mod search_modification;
mod spectra;
mod state;
mod validate;

use crate::{
    html_builder::{HtmlContent, HtmlElement, HtmlTag},
    metadata_render::{OptionalString, RenderToHtml},
    state::State,
};

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
type ModifiableState<'a> = tauri::State<'a, tokio::sync::Mutex<State>>;

#[tauri::command]
fn refresh(
    state: ModifiableState,
    theme: Theme,
) -> (usize, Option<AnnotationResult>, Vec<String>, String) {
    let state = state.blocking_lock();
    fn get_details<T: CVSource>(index: &CVIndex<T>) -> [HtmlContent; 6] {
        [
            T::cv_name().into(),
            index.version().version.clone().to_optional_string().into(),
            index.version().last_updated().to_optional_string().into(),
            HtmlTag::code
                .new()
                .content(index.version().hash_hex())
                .into(),
            index.len().to_string().into(),
            HtmlTag::div
                .new()
                .content(
                    HtmlTag::button
                        .new()
                        .class("update-ontology-file")
                        .content("File")
                        .data([("ontology", T::cv_name())])
                        .clone(),
                )
                .maybe_content(if T::cv_name() != "RESID" {
                    Some(
                        HtmlTag::button
                            .new()
                            .class("update-ontology-internet")
                            .content("Internet")
                            .data([("ontology", T::cv_name())])
                            .clone(),
                    )
                } else {
                    None
                })
                .into(),
        ]
    }

    (
        state.psm_files().len(),
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
        state.auto_open_errors.clone(),
        HtmlElement::table(
            Some(&[
                "Ontology",
                "Version",
                "Last updated",
                "Hash",
                "Number of mods",
                "Update",
            ]),
            [
                get_details(state.ontologies.unimod()),
                get_details(state.ontologies.psimod()),
                get_details(state.ontologies.gnome()),
                get_details(state.ontologies.xlmod()),
                get_details(state.ontologies.resid()),
            ],
        )
        .to_string(),
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
        MolecularFormula::pro_forma::<false, false>(text).map_err(|err| err.to_html(false))
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
        .blocking_lock()
        .psm_files()
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
    let mut state = state.blocking_lock();
    let annotated = state.psm_files().iter().find(|f| f.id == file).map_or(
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
    let mut state = state.lock().await;
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
    let (peptide, warnings) = CompoundPeptidoformIon::pro_forma(peptide, &state.ontologies)
        .map_err(|errs| {
            errs.into_iter()
                .map(|err| err.to_html(false))
                .collect::<Vec<_>>()
        })?;

    let use_charge = Charge::new::<e>(
        charge
            .or_else(|| {
                spectrum
                    .precursor()
                    .and_then(|p| p.ions.first().and_then(|i| i.charge.map(|c| c as isize)))
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
    let app_state = app.state::<Mutex<State>>();
    let mut state = app_state.blocking_lock();
    let now = chrono::Local::now().format("%Y%m%d%H%M%S%.3f");
    let (ontologies, warnings) = dbg!(Ontologies::init());
    state.ontologies = ontologies;
    if !warnings.is_empty() {
        eprintln!("Error while parsing custom modifications:\n");
        for err in &warnings {
            eprintln!("{err}");
        }

        state.custom_modifications_error = Some((
            BoxedError::small(
                mzcv::CVError::FileDoesNotExist, // Error type does not matter
                "Errors while initialising ontologies",
                "See underlying errors",
            )
            .add_underlying_errors(warnings)
            .to_html(false),
            Vec::new(),
        ))
    }
    let path = app.path().app_config_dir();
    if let Ok(path) = path {
        let custom_models = path.join(CUSTOM_MODELS_FILE);
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

fn auto_open(app: &mut tauri::App, args: Args) -> Result<(), Box<dyn std::error::Error>> {
    const RAW_EXTENSIONS: &[&str] = &["xy", "mgf", "mzml", "imzml", "mzmlb", "raw"];
    let app_state = app.state::<Mutex<State>>();
    let mut state = app_state.blocking_lock();

    for path in &args.paths {
        let actual_extension = path
            .extension()
            .map(|ex| {
                ex.eq_ignore_ascii_case("gz")
                    .then_some(path)
                    .and_then(|p| p.file_stem())
                    .and_then(|p| std::path::Path::new(p).extension())
                    .unwrap_or(ex)
            })
            .map(|ex| ex.to_string_lossy().to_lowercase());
        if let Some(ext) = actual_extension {
            if RAW_EXTENSIONS.contains(&ext.as_str()) {
                match crate::spectra::annotator_open_raw_file(path, &mut state) {
                    Ok(_) => (),
                    Err(error) => state.auto_open_errors.push(error),
                }
            } else {
                match crate::identified_peptides::annotator_open_psm_file(path, &mut state) {
                    Ok(None) => (),
                    Ok(Some(error)) => state.auto_open_errors.push(error),
                    Err(error) => state.auto_open_errors.push(error),
                }
            }
        } else {
            state.auto_open_errors.push(
                BoxedError::new(
                    BasicKind::Error,
                    "Could not identify file",
                    "An extension has to be provided for auto open file to be recognised",
                    context_error::Context::none().lines(0, path.to_string_lossy()),
                )
                .to_html(false),
            )
        }
    }

    Ok(())
}

fn setup(app: &mut tauri::App, args: Args) -> Result<(), Box<dyn std::error::Error>> {
    load_custom_mods_and_models(app)?;
    auto_open(app, args)
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

#[tauri::command]
async fn update_ontology_internet(
    ontology: &str,
    state: ModifiableState<'_>,
) -> Result<(), String> {
    let mut state = state.lock().await;
    match ontology {
        "Unimod" => state
            .ontologies
            .unimod_mut()
            .update_from_url_async(&[])
            .await
            .map_err(|err| err.to_html(true)),
        "PSI-MOD" => state
            .ontologies
            .psimod_mut()
            .update_from_url_async(&[])
            .await
            .map_err(|err| err.to_html(true)),
        "XLMOD" => state
            .ontologies
            .xlmod_mut()
            .update_from_url_async(&[])
            .await
            .map_err(|err| err.to_html(true)),
        "GNOme" => state
            .ontologies
            .gnome_mut()
            .update_from_url_async(&[])
            .await
            .map_err(|err| err.to_html(true)),
        "RESID" => Err("Cannot update RESID from the internet".to_string()),
        "Custom" => Err("Cannot update Custom from the internet".to_string()),
        _ => Err("Invalid ontology".to_string()),
    }
}

#[tauri::command]
async fn update_ontology_file(
    ontology: &str,
    files: Vec<PathBuf>,
    state: ModifiableState<'_>,
) -> Result<(), String> {
    let mut state = state.lock().await;

    // TODO: validate that the picked files make sense (correct extension and amount and put them in the right order if needed)

    let extensions = files
        .iter()
        .map(|path| {
            path.extension().map(|ex| {
                ex.eq_ignore_ascii_case("gz")
                    .then_some(path)
                    .and_then(|p| p.file_stem())
                    .and_then(|p| std::path::Path::new(p).extension())
                    .unwrap_or(ex)
            })
        })
        .collect::<Vec<_>>();

    match ontology {
        "Unimod" => {
            if extensions.len() == 1
                && extensions[0].is_some_and(|ex| ex.eq_ignore_ascii_case("xml"))
            {
                state
                    .ontologies
                    .unimod_mut()
                    .update_from_path(files.iter().map(|p| Some(p.as_ref())), false)
                    .map_err(|err| err.to_html(true))
            } else {
                Err("Provide exactly one .xml file for Unimod".to_string())
            }
        }
        "PSI-MOD" => {
            if extensions.len() == 1
                && extensions[0].is_some_and(|ex| ex.eq_ignore_ascii_case("obo"))
            {
                state
                    .ontologies
                    .psimod_mut()
                    .update_from_path(files.iter().map(|p| Some(p.as_ref())), false)
                    .map_err(|err| err.to_html(true))
            } else {
                Err("Provide exactly one .obo file for PSI-MOD".to_string())
            }
        }
        "XLMOD" => {
            if extensions.len() == 1
                && extensions[0].is_some_and(|ex| ex.eq_ignore_ascii_case("obo"))
            {
                state
                    .ontologies
                    .xlmod_mut()
                    .update_from_path(files.iter().map(|p| Some(p.as_ref())), false)
                    .map_err(|err| err.to_html(true))
            } else {
                Err("Provide exactly one .obo file for XLMOD".to_string())
            }
        }
        "GNOme" => {
            if extensions.len() == 2
                && extensions[0].is_some_and(|ex| ex.eq_ignore_ascii_case("obo"))
                && extensions[1].is_some_and(|ex| ex.eq_ignore_ascii_case("csv"))
            {
                state
                    .ontologies
                    .gnome_mut()
                    .update_from_path(files.iter().map(|p| Some(p.as_ref())), false)
                    .map_err(|err| err.to_html(true))
            } else if extensions.len() == 2
                && extensions[0].is_some_and(|ex| ex.eq_ignore_ascii_case("csv"))
                && extensions[1].is_some_and(|ex| ex.eq_ignore_ascii_case("obo"))
            {
                state
                    .ontologies
                    .gnome_mut()
                    .update_from_path(
                        [
                            Some(std::path::Path::new(&files[1])),
                            Some(std::path::Path::new(&files[0])),
                        ],
                        false,
                    )
                    .map_err(|err| err.to_html(true))
            } else {
                Err("Provide one .obo and one .csv file for GNOme".to_string())
            }
        }
        "RESID" => {
            if extensions.len() == 1
                && extensions[0].is_some_and(|ex| ex.eq_ignore_ascii_case("xml"))
            {
                state
                    .ontologies
                    .resid_mut()
                    .update_from_path(files.iter().map(|p| Some(p.as_ref())), false)
                    .map_err(|err| err.to_html(true))
            } else {
                Err("Provide exactly one .xml file for RESID".to_string())
            }
        }
        "Custom" => {
            if extensions.len() == 1
                && extensions[0].is_some_and(|ex| ex.eq_ignore_ascii_case("xml"))
            {
                state
                    .ontologies
                    .custom_mut()
                    .update_from_path(files.iter().map(|p| Some(p.as_ref())), false)
                    .map_err(|err| err.to_html(true))
            } else {
                Err("Provide exactly one .json file for Custom".to_string())
            }
        }
        _ => Err("Invalid ontology".to_string()),
    }
}

fn main() {
    let args = Args::parse();
    tauri::Builder::default()
        .plugin(tauri_plugin_dialog::init())
        .plugin(tauri_plugin_opener::init())
        .manage(Mutex::new(State {
            spectra: Vec::new(),
            identified_peptide_files: std::cell::RefCell::new(Vec::new()),
            annotated_spectrum: None,
            ontologies: Ontologies::empty(),
            custom_modifications_error: None,
            custom_models: Vec::new(),
            custom_models_error: None,
            auto_open_errors: Vec::new(),
        }))
        .setup(|app| setup(app, args))
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
            load_annotated_spectrum,
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
            spectra::select_retention_time,
            spectra::select_spectrum_index,
            spectra::select_spectrum_native_id,
            update_ontology_internet,
            update_ontology_file,
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

#[derive(Parser)]
struct Args {
    /// The paths to open in the annotator
    paths: Vec<std::path::PathBuf>,
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
