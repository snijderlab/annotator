use std::{
    cell::{Ref, RefCell, RefMut},
    fs::File,
    sync::atomic::AtomicUsize,
};

use mzdata::{
    io::{MZReaderType, SpectrumSource},
    prelude::SpectrumLike,
    spectrum::MultiLayerSpectrum,
};
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

pub enum RawFile {
    File {
        id: usize,
        rawfile: MZReaderType<File>,
        selected_spectra: Vec<usize>,
        path: String,
    },
    Single {
        id: usize,
        spectrum: MultiLayerSpectrum,
        selected: bool,
        title: String,
    },
}

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct RawFileDetails {
    pub id: usize,
    pub path: String,
    pub spectra: usize,
    pub single: bool,
    pub selected: bool,
}

impl RawFile {
    pub fn id(&self) -> usize {
        match self {
            RawFile::Single { id, .. } | RawFile::File { id, .. } => *id,
        }
    }

    pub fn clear_selected(&mut self) {
        match self {
            Self::File {
                selected_spectra, ..
            } => selected_spectra.clear(),
            Self::Single { selected, .. } => *selected = false,
        }
    }

    pub fn unselect_index(&mut self, index: usize) {
        match self {
            Self::File {
                selected_spectra, ..
            } => {
                if let Some(index) = selected_spectra.iter().position(|i| *i == index) {
                    selected_spectra.remove(index);
                }
            }
            Self::Single { selected, .. } => {
                if index == 0 {
                    *selected = false;
                }
            }
        }
    }

    pub fn select_index(&mut self, index: usize) -> Result<(), &'static str> {
        match self {
            Self::File {
                rawfile,
                selected_spectra,
                ..
            } => {
                if !selected_spectra.contains(&index) {
                    if rawfile.get_spectrum_by_index(index).is_some() {
                        selected_spectra.push(index);
                        selected_spectra.sort();
                        Ok(())
                    } else {
                        Err("Spectrum index does not exist")
                    }
                } else {
                    Ok(())
                }
            }
            Self::Single { selected, .. } => {
                if index == 0 {
                    *selected = true;
                    Ok(())
                } else {
                    Err("Outside of file range")
                }
            }
        }
    }

    pub fn select_native_id(&mut self, native_id: String) -> Result<(), &'static str> {
        match self {
            Self::File {
                rawfile,
                selected_spectra,
                ..
            } => rawfile
                .get_spectrum_by_id(&native_id)
                .map(|s| {
                    if !selected_spectra.contains(&s.index()) {
                        selected_spectra.push(s.index());
                        selected_spectra.sort();
                    }
                })
                .ok_or("Native ID does not exist"),
            Self::Single {
                selected, spectrum, ..
            } => {
                if spectrum.description.id == native_id {
                    *selected = true;
                    Ok(())
                } else {
                    Err("Native ID does not exist")
                }
            }
        }
    }

    pub fn new_file(path: &str, file: MZReaderType<File>) -> Self {
        RawFile::File {
            id: RAW_FILE_COUNTER.fetch_add(1, std::sync::atomic::Ordering::SeqCst),
            rawfile: file,
            selected_spectra: Vec::new(),
            path: path.to_string(),
        }
    }

    pub fn new_single(spectrum: MultiLayerSpectrum, title: String) -> Self {
        RawFile::Single {
            id: RAW_FILE_COUNTER.fetch_add(1, std::sync::atomic::Ordering::SeqCst),
            spectrum,
            selected: true,
            title,
        }
    }

    pub fn details(&self) -> RawFileDetails {
        match self {
            RawFile::File {
                id, rawfile, path, ..
            } => RawFileDetails {
                id: *id,
                path: path.clone(),
                spectra: rawfile.len(),
                single: false,
                selected: false,
            },
            RawFile::Single {
                id,
                title,
                selected,
                ..
            } => RawFileDetails {
                id: *id,
                path: title.clone(),
                spectra: 1,
                single: true,
                selected: *selected,
            },
        }
    }

    pub fn get_selected_spectra(&mut self) -> Box<dyn Iterator<Item = MultiLayerSpectrum> + '_> {
        match self {
            Self::File {
                rawfile,
                selected_spectra,
                ..
            } => Box::new(
                selected_spectra
                    .iter()
                    .filter_map(|index| rawfile.get_spectrum_by_index(*index)),
            ),
            Self::Single {
                spectrum, selected, ..
            } => Box::new(std::iter::once(spectrum.clone()).take(usize::from(*selected))),
        }
    }
}

static RAW_FILE_COUNTER: AtomicUsize = AtomicUsize::new(0);
static IDENTIFIED_PEPTIDE_FILE_COUNTER: AtomicUsize = AtomicUsize::new(0);
