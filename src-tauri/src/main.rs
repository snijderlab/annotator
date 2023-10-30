#![cfg_attr(
    all(not(debug_assertions), target_os = "windows"),
    windows_subsystem = "windows"
)]

use itertools::Itertools;
use rustyms::{
    align::{align, BLOSUM62},
    error::*,
    identifications::*,
    model::*,
    spectrum::*,
    system::*,
    *,
};
use state::State;
use std::sync::Mutex;

use crate::metadata_render::RenderToHtml;
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
        Err(err) => Err(err),
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
async fn search_peptide<'a>(text: &'a str, state: ModifiableState<'a>) -> Result<String, String> {
    let state = state
        .lock()
        .map_err(|_| "Search exception: State locked".to_string())?;
    let search = ComplexPeptide::pro_forma(text)
        .map_err(|err| err.to_string())?
        .assume_linear();
    let data = state
        .peptides
        .iter()
        .enumerate()
        .map(|(index, peptide)| {
            (
                index,
                align(
                    peptide.peptide.clone(),
                    search.clone(),
                    BLOSUM62,
                    MassTolerance::Absolute(da(0.1)),
                    align::Type::GlobalForB,
                ),
                peptide,
            )
        })
        .sorted_unstable_by(|a, b| {
            b.1.normalised_score
                .partial_cmp(&a.1.normalised_score)
                .unwrap()
        })
        .filter(|(_, alignment, _)| alignment.normalised_score > 0.0)
        .take(25)
        .map(|(index, _, peptide)| {
            vec![
                index.to_string(),
                peptide.peptide.to_string(),
                peptide
                    .score
                    .map(|score| score.to_string())
                    .unwrap_or_default(),
            ]
        })
        .collect::<Vec<_>>();
    Ok(html_builder::HtmlElement::table(
        Some(&[
            "Index".to_string(),
            "Sequence".to_string(),
            "Score".to_string(),
        ]),
        &data,
    )
    .to_string())
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
            charge: peptide.metadata.charge(),
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
fn spectrum_details(index: usize, state: ModifiableState) -> String {
    state.lock().unwrap().spectra.get(index).map_or(
        "Spectrum index not valid".to_string(),
        |spectrum| {
            format!(
                "{}\n{:.3}@{:.3}{:+.0}{:.3}{}",
                spectrum.title,
                spectrum.mass.value,
                spectrum.rt.value,
                spectrum.charge.value,
                spectrum
                    .intensity
                    .map_or(String::new(), |i| format!(" I:{i}")),
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

    state.lock().unwrap().spectra = vec![RawSpectrum {
        title: "Clipboard".to_string(),
        charge: Charge::new::<e>(1.0),
        spectrum,
        ..RawSpectrum::default()
    }];
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

type ModelParameters = Vec<(Location, String)>;

#[allow(clippy::too_many_arguments)]
#[tauri::command]
fn annotate_spectrum(
    index: usize,
    ppm: f64,
    charge: Option<f64>,
    noise_threshold: Option<f64>,
    model: &str,
    peptide: &str,
    state: ModifiableState,
    cmodel: ModelParameters,
) -> Result<(String, String, String), CustomError> {
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
        "etcid" => Model::etcid(),
        "cidhcd" => Model::cid_hcd(),
        "etd" => Model::etd(),
        "custom" => Model::new(
            get_model_param(&cmodel[0])?,
            get_model_param(&cmodel[1])?,
            get_model_param(&cmodel[2])?,
            get_model_param(&cmodel[3])?,
            get_model_param(&cmodel[4])?,
            get_model_param(&cmodel[5])?,
            get_model_param(&cmodel[6])?,
            get_model_param(&cmodel[7])?,
            get_model_param(&cmodel[8])?,
            cmodel[9]
                .1
                .to_owned()
                .split(',')
                .filter(|n| !n.is_empty())
                .map(|n| n.parse::<NeutralLoss>())
                .collect::<Result<Vec<NeutralLoss>, _>>()?,
            MassOverCharge::new::<mz>(ppm),
            None,
        ),
        _ => Model::all(),
    };
    model.ppm = MassOverCharge::new::<mz>(ppm);
    let peptide = rustyms::ComplexPeptide::pro_forma(peptide)?;
    let multiple_peptides = peptide.peptides().len() != 1;
    let mut spectrum = state.spectra[index].clone();
    if let Some(threshold) = noise_threshold {
        spectrum.noise_filter(threshold);
    }
    let use_charge = charge.map_or(spectrum.charge, Charge::new::<e>);
    let fragments = peptide
        .generate_theoretical_fragments(use_charge, &model)
        .ok_or(CustomError::error(
            "Undefined mass",
            "The sequence requested does not have a defined mass (you used B/Z amino acids)",
            Context::none(),
        ))?;
    let annotated = spectrum.annotate(peptide, &fragments, &model, MassMode::Monoisotopic);
    Ok((
        render::annotated_spectrum(&annotated, "spectrum", &fragments),
        render::fragment_table(&fragments, multiple_peptides),
        format!("{annotated:#?}\n{model:#?}"),
    ))
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
            spectrum_details,
            search_peptide,
            identified_peptide_details,
            load_identified_peptide,
            annotate_spectrum,
        ])
        .run(tauri::generate_context!())
        .expect("error while running tauri application");
}
