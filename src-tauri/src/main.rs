#![cfg_attr(
    all(not(debug_assertions), target_os = "windows"),
    windows_subsystem = "windows"
)]

use itertools::Itertools;
use rustyms::error::{Context, CustomError};
use rustyms::MassOverCharge;
use rustyms::{e, Charge, Location, Mass, Model, NeutralLoss, RawPeak, RawSpectrum, Time, Zero};
use state::State;
use std::sync::Mutex;

mod html_builder;
mod render;
mod state;

type ModifiableState<'a> = tauri::State<'a, std::sync::Mutex<State>>;

#[tauri::command]
fn load_mgf(path: &str, state: ModifiableState) -> Result<usize, String> {
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
                    .map(|s| format!("\n{s}"))
                    .unwrap_or_default()
            )
        },
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
        num_scans: 0,
        rt: Time::zero(),
        charge: Charge::new::<e>(1.0),
        mass: Mass::zero(),
        intensity: None,
        spectrum,
        sequence: None,
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
        if let (Ok(mz), Ok(intensity)) = (cells[1].parse(), cells[4].parse()) {
            spectrum.push(RawPeak {
                mz: MassOverCharge::new::<rustyms::mz>(mz),
                intensity,
                charge: Charge::zero(),
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
        if let (Ok(mz), Ok(intensity)) = (cells[0].parse(), cells[2].parse()) {
            spectrum.push(RawPeak {
                mz: MassOverCharge::new::<rustyms::mz>(mz),
                intensity,
                charge: Charge::zero(),
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
        if let (Ok(mz), Ok(intensity)) = (cells[0].parse(), cells[1].parse()) {
            spectrum.push(RawPeak {
                mz: MassOverCharge::new::<rustyms::mz>(mz),
                intensity,
                charge: Charge::zero(),
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
    use rustyms::mz;
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
    let annotated = spectrum.annotate(peptide, &fragments, &model);
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
        }))
        .invoke_handler(tauri::generate_handler![
            load_mgf,
            annotate_spectrum,
            load_clipboard,
            spectrum_details
        ])
        .run(tauri::generate_context!())
        .expect("error while running tauri application");
}
