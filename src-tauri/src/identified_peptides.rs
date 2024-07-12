use itertools::Itertools;
use rustyms::{
    align::{align, matrix::BLOSUM62},
    error::*,
    identification::*,
    system::da,
    *,
};
use serde::{Deserialize, Serialize};

use crate::{html_builder, state::IdentifiedPeptideFile, ModifiableState};
use rayon::prelude::*;

#[tauri::command]
pub async fn load_identified_peptides_file<'a>(
    path: &'a str,
    state: ModifiableState<'a>,
) -> Result<usize, CustomError> {
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
        Some("csv") => {
            PeaksData::parse_file(path)
                .map(|peptides| {
                    state.lock().unwrap().identified_peptide_files.push(
                        IdentifiedPeptideFile::new(
                            path.to_string(),
                            peptides.filter_map(|p| p.ok()).map(|p| p.into()).collect(),
                        ),
                    );
                })
                .or_else(|_| {
                    NovorData::parse_file(path).map(|peptides| {
                        state.lock().unwrap().identified_peptide_files.push(
                            IdentifiedPeptideFile::new(
                                path.to_string(),
                                peptides.filter_map(|p| p.ok()).map(|p| p.into()).collect(),
                            ),
                        )
                    })
                })
                .map_err(|_| {
                    CustomError::error(
                        "Unknown file",
                        "Could not be recognised as either a Peaks or Novor file",
                        Context::None,
                    )
                })
        }
        Some("tsv") => {
            MSFraggerData::parse_file(path)
                .map(|peptides| {
                    state
                        .lock()
                        .unwrap()
                        .identified_peptide_files
                        .push(IdentifiedPeptideFile::new(
                            path.to_string(),
                            peptides.filter_map(|p| p.ok()).map(|p| p.into()).collect(),
                        ))
                })
                .or_else(|_| {
                    SageData::parse_file(path).map(|peptides| {
                        state.lock().unwrap().identified_peptide_files.push(
                            IdentifiedPeptideFile::new(
                                path.to_string(),
                                peptides.filter_map(|p| p.ok()).map(|p| p.into()).collect(),
                            ),
                        )
                    })
                })
                .map_err(|_| {
                    CustomError::error(
                        "Unknown file",
                        "Could not be recognised as either a MSFragger or Sage file",
                        Context::None,
                    )
                })
        }
        Some("psmtsv") => OpairData::parse_file(path).map(|peptides| {
            state
                .lock()
                .unwrap()
                .identified_peptide_files
                .push(IdentifiedPeptideFile::new(
                    path.to_string(),
                    peptides.filter_map(|p| p.ok()).map(|p| p.into()).collect(),
                ))
        }),
        Some("fasta") => FastaData::parse_file(path).map(|peptides| {
            state
                .lock()
                .unwrap()
                .identified_peptide_files
                .push(IdentifiedPeptideFile::new(
                    path.to_string(),
                    peptides.into_iter().map(|p| p.into()).collect(),
                ))
        }),
        _ => Err(CustomError::error(
            "Unknown extension",
            "Use CSV, TSV, PSMTSV, or Fasta, or any of these as a gzipped file (eg csv.gz).",
            Context::None,
        )),
    }?;
    Ok(state.lock().unwrap().identified_peptide_files.len())
}

#[tauri::command]
pub async fn close_identified_peptides_file(
    file: usize,
    state: ModifiableState<'_>,
) -> Result<(), CustomError> {
    let mut state = state.lock().map_err(|_| {
        CustomError::error(
            "Could not lock mutex",
            "You are likely doing to many things in parallel",
            Context::none(),
        )
    })?;
    if let Some(pos) = state
        .identified_peptide_files
        .iter()
        .position(|f| f.id == file)
    {
        state.identified_peptide_files.remove(pos);
        Ok(())
    } else {
        Err(CustomError::error("File does not exist", "This selected file could not be closed as it does not exist, did you already close it?", Context::none()))
    }
}

/// Get the id, file name, path, and numbe rof peptides of all open identified peptides files
#[tauri::command]
pub async fn get_identified_peptides_files(
    state: ModifiableState<'_>,
) -> Result<Vec<(usize, String, String, usize)>, ()> {
    Ok(state
        .lock()
        .map_err(|_| ())?
        .identified_peptide_files
        .iter()
        .map(|f| (f.id, f.file_name(), f.path.clone(), f.peptides.len()))
        .collect_vec())
}

#[tauri::command]
pub async fn search_peptide<'a>(
    text: &'a str,
    minimal_match_score: f64,
    minimal_peptide_score: f64,
    state: ModifiableState<'a>,
) -> Result<String, CustomError> {
    let state = state.lock().map_err(|_| {
        CustomError::error(
            "Cannot search",
            "The state is locked, are you trying to do many things at the same time?",
            Context::None,
        )
    })?;
    let search = LinearPeptide::<Linked>::pro_forma(text, Some(&state.database))?
        .simple()
        .ok_or_else(|| {
            CustomError::error(
                "Invalid search peptide",
                "A search peptide should be simple",
                Context::None,
            )
        })?;
    let data = state
        .identified_peptide_files
        .iter()
        .flat_map(|file| {
            file.peptides
                .iter()
                .enumerate()
                .map(|(index, p)| (file.id, index, p))
        })
        .filter(|(_, _, p)| p.score.map_or(true, |score| score >= minimal_peptide_score))
        .par_bridge()
        .map(|(id, index, peptide)| {
            (
                index,
                align::<4, VerySimple, Simple>(
                    &peptide.peptide,
                    &search,
                    BLOSUM62,
                    Tolerance::new_absolute(da(0.1)),
                    align::AlignType::GLOBAL_B,
                ),
                peptide,
                id,
            )
        })
        .filter(|(_, alignment, _, _)| alignment.normalised_score() >= minimal_match_score)
        .collect::<Vec<_>>()
        .into_iter()
        .k_largest_by(25, |a, b| {
            a.1.score().normalised.cmp(&b.1.score().normalised)
        })
        .map(|(index, alignment, peptide, id)| {
            let start = alignment.start_a();
            let end = alignment.start_a() + alignment.len_a();
            vec![
                // TODO: select the right file
                format!(
                    "<a onclick=\"load_peptide({id}, {index})\">F{}:{index}</a>",
                    id + 1
                ),
                format!(
                    "{}<span class='match'>{}</span>{}",
                    peptide.peptide.sub_peptide(..start).to_string(),
                    peptide.peptide.sub_peptide(start..end).to_string(),
                    peptide.peptide.sub_peptide(end..).to_string(),
                ),
                format!("{:.3}", alignment.normalised_score()),
                peptide
                    .score
                    .map(|score| format!("{:.3}", score))
                    .unwrap_or_default(),
            ]
        })
        .collect::<Vec<_>>();
    Ok(html_builder::HtmlElement::table(
        Some(&[
            "Index".to_string(),
            "Sequence".to_string(),
            "Match Score".to_string(),
            "Peptide Score".to_string(),
        ]),
        data,
    )
    .to_string())
}

#[derive(Serialize, Deserialize)]
pub struct Settings {
    peptide: String,
    charge: Option<usize>,
    mode: Option<String>,
    scan_index: Option<usize>,
}

impl Settings {
    fn from_peptide(peptide: &IdentifiedPeptide, scan: Option<usize>) -> Self {
        Self {
            peptide: peptide.peptide.to_string(),
            charge: peptide.metadata.charge().map(|v| v.value),
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
pub fn load_identified_peptide(
    file: usize,
    index: usize,
    state: ModifiableState,
) -> Option<Settings> {
    if let Ok(state) = state.lock() {
        state
            .identified_peptide_files
            .iter()
            .find(|f| f.id == file)
            .and_then(|file| {
                file.peptides.get(index).map(|peptide| {
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
            })
    } else {
        None
    }
}
