use std::{
    cell::{Ref, RefCell, RefMut},
    fs::File,
    sync::atomic::AtomicUsize,
};

use mzdata::io::{MZReaderType, SpectrumSource};
use rustyms::{identification::IdentifiedPeptide, ontologies::CustomDatabase};
use serde::{Deserialize, Serialize};

pub struct State {
    pub spectra: Vec<RawFile>,
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

pub struct RawFile {
    pub id: usize,
    pub rawfile: MZReaderType<File>,
    pub selected_spectra: Vec<usize>,
    pub path: String,
}

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct RawFileDetails {
    id: usize,
    path: String,
    spectra: usize,
}

impl RawFile {
    pub fn new(path: &str, file: MZReaderType<File>) -> Self {
        RawFile {
            id: RAW_FILE_COUNTER.fetch_add(1, std::sync::atomic::Ordering::SeqCst),
            rawfile: file,
            selected_spectra: Vec::new(),
            path: path.to_string(),
        }
    }

    pub fn details(&self) -> RawFileDetails {
        RawFileDetails {
            id: self.id,
            path: self.path.clone(),
            spectra: self.rawfile.len(),
        }
    }
}

static RAW_FILE_COUNTER: AtomicUsize = AtomicUsize::new(0);
static IDENTIFIED_PEPTIDE_FILE_COUNTER: AtomicUsize = AtomicUsize::new(0);
