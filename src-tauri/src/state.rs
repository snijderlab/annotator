use std::{
    cell::{Ref, RefCell, RefMut},
    fs::File,
    sync::atomic::AtomicUsize,
};

use mzdata::io::MZReaderType;
use rustyms::{identification::IdentifiedPeptide, ontologies::CustomDatabase};

pub struct State {
    pub spectra: Option<MZReaderType<File>>,
    pub identified_peptide_files: RefCell<Vec<IdentifiedPeptideFile>>,
    pub database: CustomDatabase,
}

impl State {
    pub fn database(&self) -> Option<&CustomDatabase> {
        (!self.database.is_empty()).then_some(&self.database)
    }
    pub fn identified_peptide_files(&self) -> Ref<Vec<IdentifiedPeptideFile>> {
        self.identified_peptide_files.borrow()
    }
    pub fn identified_peptide_files_mut(&self) -> RefMut<Vec<IdentifiedPeptideFile>> {
        self.identified_peptide_files.borrow_mut()
    }
}

pub struct IdentifiedPeptideFile {
    pub id: usize,
    pub path: String,
    pub peptides: Vec<IdentifiedPeptide>,
}

impl IdentifiedPeptideFile {
    pub fn file_name(&self) -> String {
        std::path::Path::new(&self.path)
            .file_name()
            .unwrap()
            .to_string_lossy()
            .to_string()
    }

    pub fn new(path: String, peptides: Vec<IdentifiedPeptide>) -> Self {
        Self {
            id: IDENTIFIED_PEPTIDE_FILE_COUNTER.fetch_add(1, std::sync::atomic::Ordering::SeqCst),
            path,
            peptides,
        }
    }
}

static IDENTIFIED_PEPTIDE_FILE_COUNTER: AtomicUsize = AtomicUsize::new(0);
