#![cfg_attr(
    all(not(debug_assertions), target_os = "windows"),
    windows_subsystem = "windows"
)]

use itertools::Itertools;
use mass_alignment::{template::Template, *};
use pdbtbx::*;
use rustyms::{
    e, generate_theoretical_fragments, AverageWeight, Charge, Location, Model, NeutralLoss,
};
use rustyms::{mz, Hecklib, MassOverCharge, MonoIsotopic};
use state::State;
use std::collections::HashMap;
use std::sync::Mutex;

mod html_builder;
mod render;
mod state;

// Learn more about Tauri commands at https://tauri.app/v1/guides/features/command
#[tauri::command]
fn align_sequences(template: &str, reads: &str, alignment_type: &str) -> String {
    let alphabet = Alphabet::default();
    let template = sequence_from_string(template);
    let reads: Vec<Vec<AminoAcid>> = reads.split('\n').map(sequence_from_string).collect();
    let alignment_type = match alignment_type {
        "1" => Type::Local,
        "2" => Type::GlobalForB,
        "3" => Type::Global,
        _ => panic!("Incorrect alignment type"),
    };

    let result = Template::new(
        template,
        reads.iter().map(|a| a.as_slice()).collect(),
        &alphabet,
        alignment_type,
    );
    result.generate_html()
}

#[tauri::command]
fn load_cif(path: &str, min_length: usize, warn: bool) -> Result<(String, String), String> {
    let result = open(path, StrictnessLevel::Loose);
    if let Ok(file) = result {
        let warnings = file.1.into_iter().map(|err| format!("{}", err)).join("\n");
        let pdb = file.0;
        let mut found_unknown = HashMap::new();
        let output = pdb
            .chains()
            .map(|c| {
                c.conformers()
                    .filter_map(|a| {
                        match AMINO_ACIDS
                            .iter()
                            .position(|err| *err == a.name())
                            .and_then(|v| AMINO_ACIDS_CHAR.get(v))
                        {
                            Some(s) => Some(s),
                            None => {
                                if warn && !IGNORE_LIST.contains(&a.name()) {
                                    found_unknown.insert(
                                        a.name(),
                                        1 + found_unknown.get(a.name()).unwrap_or(&0),
                                    );
                                };
                                None
                            }
                        }
                    })
                    .collect::<String>()
            })
            .filter(|a| a.len() >= min_length)
            .join("\n");
        let warnings = warnings + "\n" + &found_unknown.into_iter().map(|name|  {
            format!(
                "{}",
                PDBError::new(
                    ErrorLevel::GeneralWarning,
                    "Unrecognised residue",
                    format!(
                        "This name was not recognised as an Amino Acid or common solvent. It was found {} time{}.",
                        name.1,
                        if name.1 != 1 { "s" } else { "" }
                    ),
                    Context::show(name.0),
                )
            )
        }).join("\n");
        Ok((output, warnings))
    } else {
        Err(result
            .unwrap_err()
            .into_iter()
            .map(|a| format!("{}", a))
            .collect())
    }
}

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

#[tauri::command]
fn annotate_spectrum(
    index: usize,
    ppm: f64,
    mass: &str,
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
        "really_all" => Model::really_all(),
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
    // if let Some(threshold) = noise_threshold {
    //     spectrum.noise_filter(threshold);
    // }
    let use_charge = charge.map_or(spectrum.charge, Charge::new::<e>);
    let fragments = if mass == "monoisotopic" {
        generate_theoretical_fragments::<MonoIsotopic>(&peptide, use_charge, &model)
    } else if mass == "averageweight" {
        generate_theoretical_fragments::<AverageWeight>(&peptide, use_charge, &model)
    } else if mass == "hecklib" {
        generate_theoretical_fragments::<Hecklib>(&peptide, use_charge, &model)
    } else {
        panic!("Unknown mass type");
    };
    let annotated = spectrum.annotate(peptide, &fragments, &model);
    Ok((
        render::annotated_spectrum(&annotated, "spectrum", &fragments),
        render::fragment_table(&fragments),
        format!("{annotated:#?}\n{model:#?}"),
    ))
}

/// All amino acids. Includes Amber-specific naming conventions for (de-)protonated versions, CYS involved in
/// disulfide bonding and the like.
const AMINO_ACIDS: &[&str] = &[
    "ALA", "ARG", "ASH", "ASN", "ASP", "ASX", "CYS", "CYX", "GLH", "GLN", "GLU", "GLY", "HID",
    "HIE", "HIM", "HIP", "HIS", "ILE", "LEU", "LYN", "LYS", "MET", "PHE", "PRO", "SER", "THR",
    "TRP", "TYR", "VAL", "SEC", "PYL",
];

const AMINO_ACIDS_CHAR: &[char] = &[
    'A', 'R', 'N', 'N', 'D', 'B', 'C', 'C', 'Q', 'Q', 'E', 'G', 'H', 'H', 'H', 'H', 'H', 'I', 'L',
    'K', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'U', 'O',
];

const IGNORE_LIST: &[&str] = &["HOH", "WAT", "ADP", "DMS"]; // Common solvents I recognised

fn main() {
    tauri::Builder::default()
        .manage(Mutex::new(State {
            spectra: Vec::new(),
        }))
        .invoke_handler(tauri::generate_handler![
            align_sequences,
            load_cif,
            load_mgf,
            annotate_spectrum
        ])
        .run(tauri::generate_context!())
        .expect("error while running tauri application");
}

#[test]
fn deserialise_test() {
    dbg!(serde_json::to_string(&vec![
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
