#![cfg_attr(
    all(not(debug_assertions), target_os = "windows"),
    windows_subsystem = "windows"
)]

use rustyms::{e, Charge, Location, Model, NeutralLoss};
use rustyms::{mz, MassOverCharge};
use state::State;
use std::sync::Mutex;

mod html_builder;
mod render;
mod state;

type ModifiableState<'a> = tauri::State<'a, std::sync::Mutex<State>>;

// Learn more about Tauri commands at https://tauri.app/v1/guides/features/command
#[tauri::command]
fn load_mgf(path: &str, state: ModifiableState) -> Result<String, String> {
    //match rustyms::thermo::open(path) {
    //    Ok(()) => Ok("loaded thermo".to_string()),
    //    Err(err) => Err(err),
    //}
    match rustyms::rawfile::mgf::open(path) {
        Ok(v) => {
            let count = v.len();
            state.lock().unwrap().spectra = v;
            Ok(format!("Loaded {count} spectra"))
        }
        Err(err) => Err(err),
    }
}

type ModelParameters = Vec<(Location, Vec<NeutralLoss>)>;

#[allow(clippy::too_many_arguments)]
#[tauri::command]
fn annotate_spectrum(
    index: usize,
    ppm: f64,
    _mass: &str,
    charge: Option<f64>,
    noise_threshold: Option<f64>,
    model: &str,
    peptide: &str,
    state: ModifiableState,
    cmodel: ModelParameters,
) -> Result<(String, String, String), String> {
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
        .invoke_handler(tauri::generate_handler![load_mgf, annotate_spectrum])
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
