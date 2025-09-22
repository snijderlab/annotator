use std::{
    cell::{Ref, RefCell, RefMut},
    fs::File,
    ops::RangeInclusive,
    sync::{Arc, OnceLock, atomic::AtomicUsize},
};

use mzdata::{
    io::{MZReaderType, RandomAccessSpectrumIterator, SpectrumSource},
    prelude::SpectrumLike,
    spectrum::MultiLayerSpectrum,
};
use ordered_float::OrderedFloat;
use rustyms::{
    align::AlignIndex,
    identification::{IdentifiedPeptidoform, MaybePeptidoform, MetaData},
    ontology::CustomDatabase,
    prelude::*,
    sequence::{Linked, SimpleLinear},
    system::OrderedTime,
};
use serde::{Deserialize, Serialize};

use crate::raw_data::RawFile;

pub struct State {
    pub spectra: Vec<RawFile>,
    pub identified_peptide_files: RefCell<Vec<IdentifiedPeptidoformFile>>,
    pub custom_modifications: CustomDatabase,
    pub custom_modifications_error: Option<(String, Vec<String>)>,
    pub custom_models: Vec<(String, FragmentationModel)>,
    pub custom_models_error: Option<(String, Vec<String>)>,
}

impl State {
    pub fn database(&self) -> Option<&CustomDatabase> {
        (!self.custom_modifications.is_empty()).then_some(&self.custom_modifications)
    }
    pub fn identified_peptide_files(&self) -> Ref<Vec<IdentifiedPeptidoformFile>> {
        self.identified_peptide_files.borrow()
    }
    pub fn identified_peptide_files_mut(&self) -> RefMut<Vec<IdentifiedPeptidoformFile>> {
        self.identified_peptide_files.borrow_mut()
    }
    pub fn spectra_and_models(&mut self) -> (&mut [RawFile], &[(String, FragmentationModel)]) {
        (&mut self.spectra, &self.custom_models)
    }
}

#[derive(Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct IndexSequence {
    pub sequence: Arc<Peptidoform<SimpleLinear>>,
    pub score: Option<OrderedFloat<f64>>,
    pub id: usize,
    pub index: usize,
}

impl HasPeptidoformImpl for IndexSequence {
    type Complexity = SimpleLinear;
    fn peptidoform(&self) -> &Peptidoform<Self::Complexity> {
        &self.sequence
    }
}

pub struct IdentifiedPeptidoformFile {
    pub id: usize,
    pub path: String,
    pub peptides: Vec<IdentifiedPeptidoform<Linked, MaybePeptidoform>>,
    pub index: OnceLock<AlignIndex<4, IndexSequence>>,
}

impl IdentifiedPeptidoformFile {
    pub fn file_name(&self) -> String {
        std::path::Path::new(&self.path)
            .file_name()
            .unwrap()
            .to_string_lossy()
            .to_string()
    }

    pub fn new(
        path: String,
        peptides: Vec<IdentifiedPeptidoform<Linked, MaybePeptidoform>>,
    ) -> Self {
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
                            .and_then(|p| p.into_simple_linear())
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
