use std::{fs::File, ops::RangeInclusive, sync::atomic::AtomicUsize};

use mzdata::{
    io::{MZReaderType, RandomAccessSpectrumIterator, SpectrumSource},
    prelude::SpectrumLike,
    spectrum::MultiLayerSpectrum,
};
use mzpeaks::PeakCollection;
use rustyms::{
    quantities::Tolerance,
    system::{MassOverCharge, OrderedTime, Time},
};
use serde::{Deserialize, Serialize};

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
                .start_from_time(rt.start().value) // TODO: likely mutates the file so fix
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

    pub fn get_ms1(
        &mut self,
        precursor: MassOverCharge,
        rt: RangeInclusive<OrderedTime>,
        tolerance: Tolerance<MassOverCharge>,
    ) -> Option<Vec<(Time, f32)>> {
        match self {
            Self::File { rawfile, .. } => Some(
                rawfile
                    .iter()
                    .filter(|s| s.ms_level() == 1)
                    .skip_while(|s| {
                        s.start_time() <= rt.start().get::<rustyms::system::time::min>()
                    })
                    .take_while(|s| s.start_time() <= rt.end().get::<rustyms::system::time::min>())
                    .filter_map(|mut s| {
                        let start = s.start_time();
                        // Some((Time::new::<rustyms::system::s>(start), 0.0))

                        s.pick_peaks(0.5);
                        s.peaks.and_then(|c| {
                            c._closest_peak(precursor.value, tolerance.into(), 0)
                                .map(|i| (Time::new::<rustyms::system::s>(start), c[i].intensity))
                        })
                    })
                    .collect(),
            ),
            _ => None,
        }
    }
}

static RAW_FILE_COUNTER: AtomicUsize = AtomicUsize::new(0);
