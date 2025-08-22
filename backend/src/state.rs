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

#[derive(Clone, Debug, Deserialize, Serialize)]
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

    pub fn select_retention_time(
        &mut self,
        rt: RangeInclusive<OrderedTime>,
    ) -> Result<(), &'static str> {
        match self {
            Self::File {
                rawfile,
                selected_spectra,
                ..
            } => rawfile
                .start_from_time(rt.start().value)
                .map(|iter| {
                    iter.take_while(|s| s.start_time() <= rt.end().value)
                        .for_each(|s| {
                            if s.ms_level() == 2 && !selected_spectra.contains(&s.index()) {
                                selected_spectra.push(s.index());
                                selected_spectra.sort();
                            }
                        })
                })
                .map_err(|_| "Retention time outside of range"),
            Self::Single { .. } => {
                Err("Cannot select based on retention time for a single spectrum")
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
