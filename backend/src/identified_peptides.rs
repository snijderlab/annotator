use crate::{html_builder, state::IdentifiedPeptideFile, ModifiableState};
use align::AlignScoring;
use itertools::Itertools;
use rayon::prelude::*;
use rustyms::{align::align, error::*, identification::*, *};
use serde::{Deserialize, Serialize};

/// Open a file and get all individual peptide errors.
/// # Errors
/// When the file could not be opened correctly.
#[tauri::command]
pub async fn load_identified_peptides_file<'a>(
    path: &'a str,
    state: ModifiableState<'a>,
) -> Result<Option<CustomError>, CustomError> {
    let state = state.lock().unwrap();
    let mut peptide_errors = Vec::new();
    let peptides = open_identified_peptides_file(path, state.database())?;
    state
        .identified_peptide_files_mut()
        .push(IdentifiedPeptideFile::new(
            path.to_string(),
            peptides
                .filter_map(|p| match p {
                    Ok(p) => Some(p),
                    Err(e) => {
                        peptide_errors.push(e);
                        None
                    }
                })
                .collect(),
        ));
    if peptide_errors.is_empty() {
        Ok(None)
    } else {
        Ok(Some(
            CustomError::warning(
                "Could not parse all peptides",
                format!(
                    "Out of all peptides {} gave rise to errors while parsing",
                    peptide_errors.len()
                ),
                Context::show(path),
            )
            .with_underlying_errors(peptide_errors),
        ))
    }
}

#[tauri::command]
pub async fn close_identified_peptides_file(
    file: usize,
    state: ModifiableState<'_>,
) -> Result<(), CustomError> {
    state.lock().map_err(|_| {
        CustomError::error(
            "Could not lock mutex",
            "You are likely doing too many things in parallel",
            Context::none(),
        )
    }).and_then(|state| {
            let pos = state
            .identified_peptide_files()
            .iter()
            .position(|f| f.id == file);
            if let Some(pos) = pos
            {
                state.identified_peptide_files_mut().remove(pos);
                Ok(())
            } else {
                Err(CustomError::error("File does not exist", "This selected file could not be closed as it does not exist, did you already close it?", Context::none()))
            }
        }
    )
}

/// Get the id, file name, path, and number of peptides of all open identified peptides files
#[tauri::command]
pub async fn get_identified_peptides_files(
    state: ModifiableState<'_>,
) -> Result<Vec<(usize, String, String, usize)>, ()> {
    Ok(state
        .lock()
        .map_err(|_| ())?
        .identified_peptide_files()
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
        .into_simple_linear()
        .ok_or_else(|| {
            CustomError::error(
                "Invalid search peptide",
                "A search peptide should be simple",
                Context::None,
            )
        })?;
    let data = state
        .identified_peptide_files()
        .iter()
        .flat_map(|file| {
            file.peptides
                .iter()
                .enumerate()
                .map(|(index, p)| (file.id, index, p))
        })
        .filter(|(_, _, p)| {
            p.score.map_or(true, |score| {
                score >= minimal_peptide_score && p.peptide().is_some()
            })
        })
        .par_bridge()
        .map(|(id, index, peptide)| {
            (
                index,
                align::<4, SemiAmbiguous, SimpleLinear>(
                    peptide.peptide().unwrap(),
                    &search,
                    AlignScoring::default(),
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
            let sequence = peptide.peptide().cloned().unwrap_or_default();
            vec![
                format!(
                    "<a onclick=\"load_peptide({id}, {index})\">F{}:{index}</a>",
                    id + 1
                ),
                format!(
                    "{}<span class='match'>{}</span>{}",
                    sequence.sub_peptide(..start).to_string(),
                    sequence.sub_peptide(start..end).to_string(),
                    sequence.sub_peptide(end..).to_string(),
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
pub struct IdentifiedPeptideSettings {
    peptide: String,
    charge: Option<usize>,
    mode: Option<String>,
    warning: Option<String>,
}

impl IdentifiedPeptideSettings {
    fn from_peptide(peptide: &IdentifiedPeptide, warning: Option<String>) -> Self {
        let mut str_peptide = String::new();
        if let Some(peptide) = peptide.peptide() {
            peptide.display(&mut str_peptide, true, false).unwrap()
        }
        Self {
            peptide: str_peptide,
            charge: peptide.charge().map(|v| v.value),
            mode: peptide
                .mode()
                .map(|mode| {
                    if mode.to_lowercase() == "hcd" || mode.to_lowercase() == "cid" {
                        "CidHcd"
                    } else {
                        mode
                    }
                })
                .map(|mode| mode.to_string()),
            warning,
        }
    }
}

#[tauri::command]
pub async fn load_identified_peptide(
    file: usize,
    index: usize,
    state: ModifiableState<'_>,
) -> Result<IdentifiedPeptideSettings, &'static str> {
    if let Ok(mut state) = state.lock() {
        let peptide = state
            .identified_peptide_files()
            .iter()
            .find(|f| f.id == file)
            .and_then(|file| file.peptides.get(index))
            .cloned();
        peptide
            .map(|peptide| {
                let scan_indices = peptide.scan_indices().unwrap_or_default();
                let native_ids = peptide.spectrum_native_ids().unwrap_or_default();
                let raw_file = peptide.raw_file().map(|p| p.to_owned());
                let mut message = None;

                // If there are scan indices unselect all selected spectra and select all
                if !scan_indices.is_empty() || !native_ids.is_empty() {
                    let mut name_matching = None;
                    let mut stem_matching = None;

                    let name = raw_file
                        .as_ref()
                        .and_then(|r| r.file_name())
                        .map(|n| n.to_string_lossy().to_lowercase())
                        .unwrap_or_default();
                    let stem = raw_file
                        .as_ref()
                        .and_then(|r| r.file_stem())
                        .map(|n| n.to_string_lossy().to_lowercase())
                        .unwrap_or_default();

                    // Search for a rawfile with the same name (stem+ext) or same stem, prefer the one with the same name
                    for (index, file) in state.spectra.iter_mut().enumerate() {
                        let path = file.details().path;
                        let path = std::path::Path::new(&path);
                        if name_matching.is_none()
                            && path
                                .file_name()
                                .map(|n| n.to_string_lossy().to_lowercase())
                                .unwrap_or_default()
                                == name
                        {
                            name_matching = Some(index);
                        }
                        if stem_matching.is_none()
                            && path
                                .file_stem()
                                .map(|n| n.to_string_lossy().to_lowercase())
                                .unwrap_or_default()
                                == stem
                        {
                            stem_matching = Some(index);
                        }
                        file.clear_selected();
                    }

                    if let Some(index) = name_matching.or(stem_matching) {
                        if !scan_indices.is_empty() {
                            for scan_index in scan_indices {
                                let _ = state.spectra[index].select_index(scan_index);
                            }
                        } else {
                            for native_id in native_ids {
                                let _ = state.spectra[index].select_native_id(native_id);
                            }
                        }
                    } else {
                        message = Some(format!("Could not find a raw file with name '{stem}' either load the correct raw file or manually load the spectra '{}' from the correct raw file", if !scan_indices.is_empty() {scan_indices.iter().join(",")} else {native_ids.iter().join(",")}))
                    }
                }

                IdentifiedPeptideSettings::from_peptide(&peptide, message)
            })
            .ok_or("The identified peptide could not be found")
    } else {
        Err("Could not lock mutex, are you running too many things at the same time?")
    }
}