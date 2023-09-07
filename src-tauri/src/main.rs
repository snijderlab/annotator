#![cfg_attr(
    all(not(debug_assertions), target_os = "windows"),
    windows_subsystem = "windows"
)]

use itertools::Itertools;
use rustyms::MassOverCharge;
use rustyms::{e, Charge, Location, Mass, Model, NeutralLoss, RawPeak, RawSpectrum, Time, Zero};
use state::State;
use std::sync::Mutex;

mod html_builder;
mod render;
mod state;

type ModifiableState<'a> = tauri::State<'a, std::sync::Mutex<State>>;

#[tauri::command]
fn load_mgf(path: &str, state: ModifiableState) -> Result<String, String> {
    match rustyms::rawfile::mgf::open(path) {
        Ok(v) => {
            let count = v.len();
            state.lock().unwrap().spectra = v;
            Ok(format!("Loaded {count} spectra"))
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
                "{}\n{:.3}@{:.3}{:+.0}{:.3}",
                spectrum.title,
                spectrum.mass.value,
                spectrum.rt.value,
                spectrum.charge.value,
                spectrum
                    .intensity
                    .map_or(String::new(), |i| format!(" I:{i}"))
            )
        },
    )
}

#[tauri::command]
fn load_clipboard(data: &str, state: ModifiableState) -> Result<String, String> {
    let lines = data.lines().collect_vec();
    if data.is_empty() {
        return Err("Empty clipboard".to_string());
    }
    let spectrum = match lines[0].trim() {
        "#	m/z	Res.	S/N	I	I %	FWHM" => load_bruker_clipboard(&lines),
        "m/z	Charge	Intensity	FragmentType	MassShift	Position" => load_stitch_clipboard(&lines),
        _ => Err("Not a recognised format (Bruker/Stitch)".to_string()),
    }?;

    state.lock().unwrap().spectra = vec![RawSpectrum {
        title: "Clipboard".to_string(),
        num_scans: 0,
        rt: Time::zero(),
        charge: Charge::new::<e>(1.0),
        mass: Mass::zero(),
        intensity: None,
        spectrum,
    }];
    Ok("Loaded".to_string())
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

type ModelParameters = Vec<(Location, Vec<NeutralLoss>)>;

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
) -> Result<(String, String, String), String> {
    use rustyms::mz;
    let state = state.lock().unwrap();
    if index >= state.spectra.len() {
        return Err("Non existent spectrum index".to_string());
    }
    let mut model = match model {
        "all" => Model::all(),
        "ethcd" => Model::ethcd(),
        "etcid" => Model::etcid(),
        "cidhcd" => Model::cid_hcd(),
        "etd" => Model::etd(),
        "custom" => Model::new(
            cmodel[0].to_owned(),
            cmodel[1].to_owned(),
            cmodel[2].to_owned(),
            cmodel[3].to_owned(),
            cmodel[4].to_owned(),
            cmodel[5].to_owned(),
            cmodel[6].to_owned(),
            cmodel[7].to_owned(),
            cmodel[8].to_owned(),
            cmodel[9].1.to_owned(),
            MassOverCharge::new::<mz>(ppm),
        ),
        _ => Model::all(),
    };
    model.ppm = MassOverCharge::new::<mz>(ppm);
    let peptide = rustyms::Peptide::pro_forma(peptide)?;
    let mut spectrum = state.spectra[index].clone();
    if let Some(threshold) = noise_threshold {
        spectrum.noise_filter(threshold);
    }
    let use_charge = charge.map_or(spectrum.charge, Charge::new::<e>);
    let fragments = peptide
        .generate_theoretical_fragments(use_charge, &model)
        .ok_or("The sequence requested does not have a defined mass (you used B/Z).".to_string())?;
    let annotated = spectrum.annotate(peptide, &fragments, &model);
    Ok((
        render::annotated_spectrum(&annotated, "spectrum", &fragments),
        render::fragment_table(&fragments),
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

#[test]
fn deserialise_test() {
    let _ = dbg!(serde_json::to_string(&vec![
        (Location::SkipN(1), vec![NeutralLoss::Water]),
        (Location::All, vec![NeutralLoss::Water]),
        (Location::SkipNC(1, 2), vec![NeutralLoss::Water]),
        (
            Location::TakeN { take: 1, skip: 2 },
            vec![NeutralLoss::Water]
        )
    ]));
    let text = r###"[[{"SkipN":1},["Water"]],[{"SkipN":1},["Water"]],[{"SkipN":1},["Water"]],[{"SkipN":1},["Water"]],[{"SkipN":1},["Water"]],["All",["Water"]],["All",["Water"]],["All",["Water"]],["All",["Water"]],["All",["Water"]]]"###;
    let params: Result<ModelParameters, _> = dbg!(serde_json::from_str(text));
    assert!(params.is_ok());
}
