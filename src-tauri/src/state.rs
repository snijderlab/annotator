use std::sync::atomic::AtomicUsize;

use rustyms::{identification::IdentifiedPeptide, ontologies::CustomDatabase, RawSpectrum};

pub struct State {
    pub spectra: Vec<RawSpectrum>,
    pub identified_peptide_files: Vec<IdentifiedPeptideFile>,
    pub database: CustomDatabase,
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
