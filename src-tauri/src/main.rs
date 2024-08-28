#![cfg_attr(
    all(not(debug_assertions), target_os = "windows"),
    windows_subsystem = "windows"
)]

use std::sync::Mutex;

use itertools::Itertools;
use mzdata::{
    io::SpectrumSource,
    prelude::{IonProperties, SpectrumLike},
    spectrum::{MultiLayerSpectrum, SpectrumDescription},
};
use mzpeaks::{CentroidPeak, DeconvolutedPeak};
use mzsignal::PeakPicker;
use ordered_float::OrderedFloat;
use render::{display_formula, display_mass};
use rustyms::{
    error::*,
    model::*,
    modification::SimpleModification,
    spectrum::*,
    system::{e, mz, usize::Charge, MassOverCharge},
    *,
};
use serde::{Deserialize, Serialize};
use tauri::Manager;

mod custom_modifications;
mod html_builder;
mod identified_peptides;
mod metadata_render;
mod render;
mod search_modification;
mod spectra;
mod state;

use crate::{metadata_render::RenderToHtml, state::State};

const CUSTOM_MODIFICATIONS_FILE: &str = "custom_modifications.json";
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
        MolecularFormula::from_pro_forma(text, .., false, true)
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
fn identified_peptide_details(file: usize, index: usize, state: ModifiableState) -> String {
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
                    |peptide| peptide.to_html().to_string(),
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

#[derive(Debug, PartialEq, PartialOrd, Serialize, Deserialize)]
pub struct ModelParameters {
    pub a: (Location, Vec<String>, ChargeRange),
    pub b: (Location, Vec<String>, ChargeRange),
    pub c: (Location, Vec<String>, ChargeRange),
    pub d: (Location, Vec<String>, ChargeRange),
    pub v: (Location, Vec<String>, ChargeRange),
    pub w: (Location, Vec<String>, ChargeRange),
    pub x: (Location, Vec<String>, ChargeRange),
    pub y: (Location, Vec<String>, ChargeRange),
    pub z: (Location, Vec<String>, ChargeRange),
    pub precursor: (Vec<String>, ChargeRange),
    pub immonium: (bool, ChargeRange),
    pub m: bool,
    pub modification_diagnostic: (bool, ChargeRange),
    pub modification_neutral: bool,
    pub cleave_cross_links: bool,
    pub glycan: (bool, (usize, usize), Vec<String>, ChargeRange, ChargeRange),
}

impl ModelParameters {
    fn create_model(
        self,
        name: &str,
        tolerance: (f64, &str),
        mz_range: (Option<f64>, Option<f64>),
    ) -> Result<Model, CustomError> {
        let get_model_param = |neutral_losses: &[String]| {
            neutral_losses
                .iter()
                .filter(|n| !n.is_empty())
                .map(|n| n.parse::<NeutralLoss>())
                .collect::<Result<Vec<_>, _>>()
        };
        let mut model = match name {
            "all" => Model::all(),
            "ethcd" => Model::ethcd(),
            "cidhcd" => Model::cid_hcd(),
            "etd" => Model::etd(),
            "none" => Model::none(),
            "custom" => Model::none()
                .a(PrimaryIonSeries {
                    location: self.a.0,
                    neutral_losses: get_model_param(&self.a.1)?,
                    charge_range: self.a.2,
                })
                .b(PrimaryIonSeries {
                    location: self.b.0,
                    neutral_losses: get_model_param(&self.b.1)?,
                    charge_range: self.b.2,
                })
                .c(PrimaryIonSeries {
                    location: self.c.0,
                    neutral_losses: get_model_param(&self.c.1)?,
                    charge_range: self.c.2,
                })
                .d(PrimaryIonSeries {
                    location: self.d.0,
                    neutral_losses: get_model_param(&self.d.1)?,
                    charge_range: self.d.2,
                })
                .v(PrimaryIonSeries {
                    location: self.v.0,
                    neutral_losses: get_model_param(&self.v.1)?,
                    charge_range: self.v.2,
                })
                .w(PrimaryIonSeries {
                    location: self.w.0,
                    neutral_losses: get_model_param(&self.w.1)?,
                    charge_range: self.w.2,
                })
                .x(PrimaryIonSeries {
                    location: self.x.0,
                    neutral_losses: get_model_param(&self.x.1)?,
                    charge_range: self.x.2,
                })
                .y(PrimaryIonSeries {
                    location: self.y.0,
                    neutral_losses: get_model_param(&self.y.1)?,
                    charge_range: self.y.2,
                })
                .z(PrimaryIonSeries {
                    location: self.z.0,
                    neutral_losses: get_model_param(&self.z.1)?,
                    charge_range: self.z.2,
                })
                .precursor(get_model_param(&self.precursor.0)?, self.precursor.1)
                .immonium(self.immonium)
                .m(self.m)
                .modification_specific_diagnostic_ions(self.modification_diagnostic)
                .modification_specific_neutral_losses(self.modification_neutral)
                .allow_cross_link_cleavage(self.cleave_cross_links)
                .glycan(GlycanModel {
                    allow_structural: self.glycan.0,
                    compositional_range: self.glycan.1 .0..=self.glycan.1 .1,
                    neutral_losses: get_model_param(&self.glycan.2)?,
                    oxonium_charge_range: self.glycan.3,
                    other_charge_range: self.glycan.4,
                }),
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
        if let (Some(min), Some(max)) = mz_range {
            if min > max {
                return Err(CustomError::error(
                    "m/z range invalid",
                    "The minimal value is less then the maximal value",
                    Context::None,
                ));
            }
            model.mz_range = MassOverCharge::new::<mz>(min)..=MassOverCharge::new::<mz>(max);
        }
        Ok(model)
    }
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
    model: &'a str,
    peptide: &'a str,
    custom_model: ModelParameters,
    state: ModifiableState<'a>,
    mass_mode: &'a str,
    mz_range: (Option<f64>, Option<f64>),
) -> Result<AnnotationResult, CustomError> {
    let mut state = state.lock().unwrap();
    let model = custom_model.create_model(model, tolerance, mz_range)?;
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

    let mut spectra = Vec::new();
    for file in state.spectra.iter_mut() {
        for index in &file.selected_spectra {
            let spectrum = file.rawfile.get_spectrum_by_index(*index).ok_or_else(|| {
                CustomError::error(
                    "Could not find spectrum",
                    "This spectrum index does not exist",
                    Context::None,
                )
            })?;
            spectra.push(spectrum);
        }
    }
    let mut spectrum = if spectra.is_empty() {
        return Err(CustomError::error(
            "No selected spectra",
            "Select a spectrum from an open raw file, or open a raw file if none are opened yet",
            Context::None,
        ));
    } else if spectra.len() == 1 {
        spectra.pop().unwrap()
    } else {
        let data = mzdata::spectrum::average_spectra(&spectra, 0.001);
        MultiLayerSpectrum::<CentroidPeak, DeconvolutedPeak>::new(
            SpectrumDescription::default(),
            Some(data.into()),
            None,
            None,
        )
    };
    if spectrum.peaks.is_none() {
        spectrum.denoise(filter).map_err(|err| {
            CustomError::error(
                "Spectrum could not be denoised",
                err.to_string(),
                Context::None,
            )
        })?;
        spectrum
            .pick_peaks_with(&PeakPicker::default())
            .map_err(|err| {
                CustomError::error(
                    "Spectrum could not be peak picked",
                    err.to_string(),
                    Context::None,
                )
            })?;
    };

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
    let fragments = peptide.generate_theoretical_fragments(use_charge, &model);
    let annotated = spectrum.annotate(peptide, &fragments, &model, mass_mode);
    let (spectrum, limits) = render::annotated_spectrum(
        &annotated, &spectrum, "spectrum", &fragments, &model, mass_mode,
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

fn load_custom_mods(app: &mut tauri::App) -> Result<(), Box<dyn std::error::Error>> {
    let path = app
        .handle()
        .path_resolver()
        .app_config_dir()
        .map(|dir| dir.join(CUSTOM_MODIFICATIONS_FILE));
    if let Some(path) = path {
        let state = app.state::<Mutex<State>>();
        if let Ok(data) = std::fs::read(path) {
            let mods: Vec<(usize, String, SimpleModification)> = serde_json::from_slice(&data)?;
            state
                .lock()
                .expect("Poisoned mutex at setup of custom mods")
                .database = mods;
        }
        Ok(())
    } else {
        Err("Could not find configuration file".into())
    }
}

fn main() {
    tauri::Builder::default()
        .manage(Mutex::new(State {
            spectra: Vec::new(),
            identified_peptide_files: std::cell::RefCell::new(Vec::new()),
            database: Vec::new(),
        }))
        .setup(load_custom_mods)
        .invoke_handler(tauri::generate_handler![
            annotate_spectrum,
            custom_modifications::delete_custom_modification,
            custom_modifications::get_custom_modification,
            custom_modifications::get_custom_modifications,
            custom_modifications::update_modification,
            custom_modifications::validate_custom_linker_specificity,
            custom_modifications::validate_custom_single_specificity,
            custom_modifications::validate_molecular_formula,
            custom_modifications::validate_neutral_loss,
            custom_modifications::validate_placement_rule,
            custom_modifications::validate_stub,
            details_formula,
            identified_peptide_details,
            identified_peptides::close_identified_peptides_file,
            identified_peptides::get_identified_peptides_files,
            identified_peptides::load_identified_peptide,
            identified_peptides::load_identified_peptides_file,
            identified_peptides::search_peptide,
            refresh,
            render::density_graph,
            search_modification::search_modification,
            spectra::close_raw_file,
            spectra::load_clipboard,
            spectra::load_raw,
            spectra::select_spectrum_index,
            spectra::select_spectrum_native_id,
            spectra::get_open_raw_files,
            spectra::get_selected_spectra,
            spectra::unselect_spectrum,
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
