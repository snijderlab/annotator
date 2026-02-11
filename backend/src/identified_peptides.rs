use context_error::{BasicKind, BoxedError, Context, CreateError, FullErrorContent, combine_error};
use itertools::Itertools;
use mzalign::{AlignScoring, AlignType, Alignment};
use mzcore::{prelude::*, sequence::Linked};
use mzident::*;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use crate::{ModifiableState, html_builder, psm_file::PSMFile, state::State};

/// Open a file and get all individual peptide errors.
/// # Errors
/// When the file could not be opened correctly.
#[tauri::command]
pub async fn load_identified_peptides_file<'a>(
    path: &'a str,
    state: ModifiableState<'a>,
) -> Result<Option<String>, String> {
    let mut state = state.lock().await;
    annotator_open_psm_file(std::path::Path::new(path), &mut state)
}

pub fn annotator_open_psm_file(
    path: &std::path::Path,
    state: &mut tokio::sync::MutexGuard<'_, State>,
) -> Result<Option<String>, String> {
    let mut peptide_errors = Vec::new();
    let peptides = open_psm_file(path, &state.ontologies, false).map_err(|e| e.to_html(false))?;
    state.psm_files_mut().push(PSMFile::new(
        path.to_string_lossy().to_string(),
        peptides
            .filter_map(|p| match p {
                Ok(p) => Some(p),
                Err(e) => {
                    combine_error(&mut peptide_errors, e);
                    None
                }
            })
            .collect(),
    ));
    if peptide_errors.is_empty() {
        Ok(None)
    } else {
        Ok(Some(
            BoxedError::new(
                BasicKind::Warning,
                "Could not parse all peptides",
                "All peptides with an error are ignored",
                Context::default().source(path.to_string_lossy()).to_owned(),
            )
            .add_underlying_errors(peptide_errors)
            .to_html(true),
        ))
    }
}

#[tauri::command]
pub async fn close_identified_peptides_file(
    file: usize,
    state: ModifiableState<'_>,
) -> Result<(), String> {
    let state = state.lock().await;
    let pos = state.psm_files().iter().position(|f| f.id == file);
    if let Some(pos) = pos {
        state.psm_files_mut().remove(pos);
        Ok(())
    } else {
        Err(BoxedError::new(BasicKind::Error,"File does not exist", "This selected file could not be closed as it does not exist, did you already close it?", Context::none()).to_html(false))
    }
}

/// Get the id, file name, path, and number of peptides of all open identified peptides files
#[tauri::command]
pub async fn get_identified_peptides_files(
    state: ModifiableState<'_>,
) -> Result<Vec<(usize, String, String, usize)>, ()> {
    Ok(state
        .lock()
        .await
        .psm_files()
        .iter()
        .map(|f| (f.id, f.file_name(), f.path.clone(), f.peptides.len()))
        .collect_vec())
}

#[tauri::command]
pub async fn search_peptide<'a>(
    text: &'a str,
    minimal_match_score: f64,
    minimal_peptide_score: f64,
    amount: usize,
    state: ModifiableState<'a>,
) -> Result<String, String> {
    if amount == 0 {
        return Ok(
            "Amount of peptides asked for is 0, so here are 0 peptides that match the query."
                .to_string(),
        );
    }
    let state = state.lock().await;
    let query = std::sync::Arc::new(
        Peptidoform::<Linked>::pro_forma(text, &state.ontologies)
            .map_err(|errs| {
                BoxedError::new(
                    BasicKind::Error,
                    "Invalid ProForma definition",
                    "The string could not be parsed as a ProForma definition",
                    Context::full_line(0, text),
                )
                .add_underlying_errors(errs)
            })
            .map_err(|e| e.to_html(false))?
            .0
            .into_linear()
            .ok_or_else(|| {
                BoxedError::new(
                    BasicKind::Error,
                    "Invalid search peptide",
                    "A search peptide should not be cross-linked",
                    Context::none(),
                )
                .to_html(false)
            })?,
    );
    let data = state
        .psm_files()
        .par_iter()
        .flat_map(|file| {
            file.index()
                .par_align_one_filtered(
                    query.clone(),
                    |p| p.score.is_none_or(|score| score.0 >= minimal_peptide_score),
                    AlignScoring::default(),
                    AlignType::GLOBAL_B,
                )
                .filter(|alignment| alignment.normalised_score() >= minimal_match_score)
        })
        .fold(
            || Vec::with_capacity(amount),
            |mut acc, x| {
                if acc.len() < amount {
                    let index = acc
                        .binary_search_by(|a: &Alignment<_, _>| a.cmp(&x).reverse())
                        .unwrap_or_else(|v| v);
                    acc.insert(index, x);
                } else if acc.last().is_some_and(|v| *v < x) {
                    acc.pop();
                    let index = acc
                        .binary_search_by(|a: &Alignment<_, _>| a.cmp(&x).reverse())
                        .unwrap_or_else(|v| v);
                    acc.insert(index, x);
                }
                acc
            },
        )
        .reduce(Vec::new, |mut acc, h2| {
            for x in h2 {
                if acc.len() < amount {
                    let index = acc
                        .binary_search_by(|a: &Alignment<_, _>| a.cmp(&x).reverse())
                        .unwrap_or_else(|v| v);
                    acc.insert(index, x);
                } else if acc.last().is_some_and(|v| *v < x) {
                    acc.pop();
                    let index = acc
                        .binary_search_by(|a: &Alignment<_, _>| a.cmp(&x).reverse())
                        .unwrap_or_else(|v| v);
                    acc.insert(index, x);
                }
            }
            acc
        })
        .into_iter()
        .map(|alignment| {
            let start = alignment.start_a();
            let end = alignment.start_a() + alignment.len_a();
            let sequence = alignment.seq_a().peptidoform();
            let mut highlighted_match = String::new();
            sequence
                .sub_peptidoform(..start)
                .unwrap()
                .display(&mut highlighted_match, true, false, true)
                .unwrap();
            highlighted_match.push_str("<span class='match'>");
            sequence
                .sub_peptidoform(start..end)
                .unwrap()
                .display(&mut highlighted_match, false, false, true)
                .unwrap();
            highlighted_match.push_str("</span>");
            sequence
                .sub_peptidoform(end..)
                .unwrap()
                .display(&mut highlighted_match, false, true, true)
                .unwrap();
            vec![
                format!(
                    "<a onclick=\"load_peptide({0}, {1})\">F{2}:{1}</a>",
                    alignment.seq_a().id,
                    alignment.seq_a().index,
                    alignment.seq_a().id + 1
                ),
                highlighted_match,
                format!("{:.3}", alignment.normalised_score()),
                alignment
                    .seq_a()
                    .score
                    .map(|score| format!("{score:.3}"))
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

#[derive(Deserialize, Serialize)]
pub struct PSMSettings {
    pub peptide: String,
    pub charge: Option<isize>,
    pub mode: Option<usize>,
    pub warning: Option<String>,
}

impl PSMSettings {
    fn from_peptide<C, A>(state: &State, peptide: &PSM<C, A>, warning: Option<String>) -> Self {
        let mut str_peptide = String::new();
        if let Some(peptide) = peptide.compound_peptidoform_ion() {
            peptide.display(&mut str_peptide, true).unwrap()
        }
        Self {
            peptide: str_peptide,
            charge: peptide.charge().map(|v| v.value),
            mode: peptide.mode().and_then(|mode| {
                crate::model::get_model_index(&state.custom_models, mode.as_ref())
            }),
            warning,
        }
    }
}

#[tauri::command]
pub async fn load_identified_peptide(
    file: usize,
    index: usize,
    state: ModifiableState<'_>,
) -> Result<PSMSettings, &'static str> {
    let mut state = state.lock().await;
    let peptide = state
        .psm_files()
        .iter()
        .find(|f| f.id == file)
        .and_then(|file| file.peptides.get(index))
        .cloned();
    peptide
        .map(|peptide| {
            let mut message = None;
            match peptide.scans() {
                SpectrumIds::None => (),
                SpectrumIds::FileNotKnown(scans) => {
                    if let Some(index) = state.spectra.first().map(|r| r.id()) {
                        for scan in scans {
                            let _ = match scan {
                                SpectrumId::Index(i) => state.spectra[index].select_index(i),
                                SpectrumId::Number(i) => state.spectra[index].select_index(i-1),
                                SpectrumId::Native(n) => state.spectra[index].select_native_id(n),
                                SpectrumId::RetentionTime(rt) => state.spectra[index].select_retention_time(rt),
                            };
                        }
                    } else {
                        message = Some("No raw files are loaded".to_string());
                    }
                },
                SpectrumIds::FileKnown(scans) => {
                    for (raw_file, scans) in scans {
                        let mut name_matching = None;
                        let mut stem_matching = None;

                        let name = raw_file.file_name().unwrap_or_default().to_string_lossy().to_lowercase();
                        let stem = raw_file.file_stem().unwrap_or_default().to_string_lossy().to_lowercase();

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
                            for scan in scans {
                                let _ = match scan {
                                    SpectrumId::Index(i) => state.spectra[index].select_index(i),
                                    SpectrumId::Number(i) => state.spectra[index].select_index(i-1),
                                    SpectrumId::Native(n) => state.spectra[index].select_native_id(n),
                                    SpectrumId::RetentionTime(rt) => state.spectra[index].select_retention_time(rt),
                                };
                            }
                        } else {
                            message = Some(format!("Could not find a raw file with name '{stem}' either load the correct raw file or manually load the spectra '{}' from the correct raw file", scans.iter().join(";")))
                        }
                    }
                },
            }

            PSMSettings::from_peptide(&state, &peptide, message)
        })
        .ok_or("The PSM could not be found")
}
