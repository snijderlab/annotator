#![cfg_attr(
    all(not(debug_assertions), target_os = "windows"),
    windows_subsystem = "windows"
)]

use itertools::Itertools;
use mass_alignment::{template::Template, *};
use pdbtbx::*;
use std::collections::HashMap;

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
        let warnings = file.1.into_iter().map(|e| format!("{}", e)).join("\n");
        let pdb = file.0;
        let mut found_unknown = HashMap::new();
        let output = pdb
            .chains()
            .map(|c| {
                c.conformers()
                    .filter_map(|a| {
                        match AMINO_ACIDS
                            .iter()
                            .position(|e| *e == a.name())
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
        .invoke_handler(tauri::generate_handler![align_sequences, load_cif])
        .run(tauri::generate_context!())
        .expect("error while running tauri application");
}
