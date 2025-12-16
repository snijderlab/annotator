use std::{
    cell::{Ref, RefCell, RefMut},
    sync::Arc,
};

use mzannotate::prelude::*;
use mzcore::{ontology::Ontologies, prelude::*, sequence::Linear};
use ordered_float::OrderedFloat;

use crate::{psm_file::PSMFile, raw_file::RawFile};

pub struct State {
    pub spectra: Vec<RawFile>,
    pub identified_peptide_files: RefCell<Vec<PSMFile>>,
    pub annotated_spectrum: Option<AnnotatedSpectrum>,
    pub ontologies: Ontologies,
    pub custom_modifications_error: Option<(String, Vec<String>)>,
    pub custom_models: Vec<(String, FragmentationModel)>,
    pub custom_models_error: Option<(String, Vec<String>)>,
    pub auto_open_errors: Vec<String>,
}

impl State {
    pub fn psm_files(&self) -> Ref<'_, Vec<PSMFile>> {
        self.identified_peptide_files.borrow()
    }
    pub fn psm_files_mut(&self) -> RefMut<'_, Vec<PSMFile>> {
        self.identified_peptide_files.borrow_mut()
    }
}

#[derive(Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct IndexSequence {
    pub sequence: Arc<Peptidoform<Linear>>,
    pub score: Option<OrderedFloat<f64>>,
    pub id: usize,
    pub index: usize,
}

impl HasPeptidoformImpl for IndexSequence {
    type Complexity = Linear;
    fn peptidoform(&self) -> &Peptidoform<Self::Complexity> {
        &self.sequence
    }
}
