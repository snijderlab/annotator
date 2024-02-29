use rustyms::{identification::IdentifiedPeptide, RawSpectrum};

pub struct State {
    pub spectra: Vec<RawSpectrum>,
    pub peptides: Vec<IdentifiedPeptide>,
}
