use itertools::Itertools;
use rustyms::{
    align::{align, matrix::BLOSUM62, Alignment},
    error::*,
    identification::*,
    system::da,
    *,
};
use serde::{Deserialize, Serialize};

use crate::{html_builder, ModifiableState};

#[tauri::command]
pub async fn load_identified_peptides<'a>(
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
            .map_err(|_| {
                CustomError::error(
                    "Unknown file",
                    "Could not be recognised as either a Peaks or Novor file",
                    Context::None,
                )
            }),
        Some("tsv") => SageData::parse_file(path).map(|peptides| {
            state.lock().unwrap().peptides =
                peptides.filter_map(|p| p.ok()).map(|p| p.into()).collect()
        }),
        Some("psmtsv") => OpairData::parse_file(path).map(|peptides| {
            state.lock().unwrap().peptides =
                peptides.filter_map(|p| p.ok()).map(|p| p.into()).collect()
        }),
        Some("fasta") => FastaData::parse_file(path).map(|peptides| {
            state.lock().unwrap().peptides = peptides.into_iter().map(|p| p.into()).collect()
        }),
        _ => Err(CustomError::error(
            "Unknown extension",
            "Use CSV, TSV, PSMTSV, or Fasta, or any of these as a gzipped file (eg csv.gz).",
            Context::None,
        )),
    }?;
    Ok(state.lock().unwrap().peptides.len())
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
        .peptides
        .iter()
        .filter(|p| p.score.map_or(true, |score| score >= minimal_peptide_score))
        .enumerate()
        .map(|(index, peptide)| {
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
            )
        })
        .sorted_unstable_by(|a, b| b.1.score().normalised.cmp(&a.1.score().normalised))
        .filter(|(_, alignment, _)| alignment.normalised_score() >= minimal_match_score)
        .take(25)
        .map(|(index, alignment, peptide)| {
            let start = alignment.start_a();
            let end = alignment.start_a() + alignment.len_a();
            vec![
                format!("<a onclick=\"document.getElementById('details-identified-peptide-index').value={0};document.getElementById('details-identified-peptide-index').dispatchEvent(new FocusEvent('focus'))\">{0}</a>", index.to_string()),
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
pub fn load_identified_peptide(index: usize, state: ModifiableState) -> Option<Settings> {
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