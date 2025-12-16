use std::sync::{Arc, OnceLock, atomic::AtomicUsize};

use mzalign::AlignIndex;
use mzcore::{prelude::MassMode, sequence::Linked};
use mzident::{MaybePeptidoform, PSM, PSMMetaData};
use ordered_float::OrderedFloat;

use crate::state::IndexSequence;

pub struct PSMFile {
    pub id: usize,
    pub path: String,
    pub peptides: Vec<PSM<Linked, MaybePeptidoform>>,
    pub index: OnceLock<AlignIndex<4, IndexSequence>>,
}

impl PSMFile {
    pub fn file_name(&self) -> String {
        std::path::Path::new(&self.path)
            .file_name()
            .unwrap()
            .to_string_lossy()
            .to_string()
    }

    pub fn new(path: String, peptides: Vec<PSM<Linked, MaybePeptidoform>>) -> Self {
        Self {
            id: IDENTIFIED_PEPTIDE_FILE_COUNTER.fetch_add(1, std::sync::atomic::Ordering::SeqCst),
            path,
            peptides,
            index: OnceLock::default(),
        }
    }

    pub fn index(&self) -> &AlignIndex<4, IndexSequence> {
        if let Some(index) = self.index.get() {
            index
        } else {
            let index = AlignIndex::new(
                self.peptides
                    .iter()
                    .enumerate()
                    .filter_map(|(index, identified)| {
                        identified
                            .compound_peptidoform_ion()
                            .and_then(|p| p.into_owned().singular_peptidoform())
                            .and_then(|p| p.into_linear())
                            .map(|p| IndexSequence {
                                sequence: Arc::new(p),
                                score: identified.score.map(OrderedFloat),
                                id: self.id,
                                index,
                            })
                    }),
                MassMode::Monoisotopic,
            );
            let _ = self.index.set(index);
            self.index.get().unwrap()
        }
    }
}

static IDENTIFIED_PEPTIDE_FILE_COUNTER: AtomicUsize = AtomicUsize::new(0);
