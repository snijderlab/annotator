use rustyms::{identification::IdentifiedPeptide, ontologies::CustomDatabase, RawSpectrum};

pub struct State {
    pub spectra: Vec<RawSpectrum>,
    pub peptides: Vec<IdentifiedPeptide>,
    pub database: CustomDatabase,
}
