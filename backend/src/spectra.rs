use std::{io::ErrorKind, ops::RangeInclusive, path::Path, str::FromStr};

use context_error::{BasicKind, BoxedError, Context, CreateError, FullErrorContent};
use itertools::Itertools;
use mzannotate::prelude::*;
use mzcore::system::{Mass, OrderedTime, dalton};
use mzdata::{
    Param,
    io::{
        MZFileReader, SpectrumSource,
        proxi::{PROXIError, PROXIParam},
        usi::{Identifier, USIParseError},
    },
    meta::DissociationMethodTerm,
    params::{CURIE, ParamDescribed, Unit, Value, ValueRef},
    prelude::{IonMobilityMeasure, IonProperties, PrecursorSelection, SpectrumLike},
    spectrum::{
        ArrayType, BinaryDataArrayType, DataArray, MultiLayerSpectrum, SignalContinuity,
        SpectrumDescription, SpectrumSummary, bindata::BinaryCompressionType,
    },
};
use mzpeaks::{CentroidPeak, peak_set::PeakSetVec};
use mzsignal::{ArrayPairLike, PeakPicker};
use serde::{Deserialize, Serialize};

use crate::{
    ModifiableState,
    identified_peptides::IdentifiedPeptideSettings,
    metadata_render::OptionalString,
    model::get_mzdata_model,
    render::display_mass,
    state::{RawFile, RawFileDetails},
};

#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct SelectedSpectrumDetails {
    id: usize,
    short: String,
    description: String,
}

#[tauri::command]
pub async fn load_usi<'a>(
    usi: &'a str,
    state: ModifiableState<'a>,
) -> Result<IdentifiedPeptideSettings, String> {
    let mut usi = mzdata::io::usi::USI::from_str(usi)
        .map_err(|e| match e {
            USIParseError::MalformedIndex(index, e, full) => BoxedError::new(
                BasicKind::Error,
                "Invalid USI: Invalid index",
                format!("The index is not valid: {e}"),
                Context::line(
                    None,
                    &full,
                    full.find(&index).unwrap_or_default(),
                    index.len(),
                )
                .to_owned(),
            ),
            USIParseError::MissingDataset(full) => BoxedError::new(
                BasicKind::Error,
                "Invalid USI: Missing dataset",
                "The dataset section of the USI is missing",
                Context::show(full),
            ),
            USIParseError::MissingRun(full) => BoxedError::new(
                BasicKind::Error,
                "Invalid USI: Missing run",
                "The run section of the USI is missing",
                Context::show(full),
            ),
            USIParseError::UnknownIndexType(index, full) => BoxedError::new(
                BasicKind::Error,
                "Invalid USI: Unknown index type",
                "The index has to be one of 'index', 'scan', or 'nativeId'",
                Context::line(
                    None,
                    &full,
                    full.find(&index).unwrap_or_default(),
                    index.len(),
                )
                .to_owned(),
            ),
            USIParseError::UnknownProtocol(protocol, full) => BoxedError::new(
                BasicKind::Error,
                "Invalid USI: Unknown protocol",
                "The protocol has to be 'mzspec'",
                Context::line(
                    None,
                    &full,
                    full.find(&protocol).unwrap_or_default(),
                    protocol.len(),
                )
                .to_owned(),
            ),
        })
        .map_err(|err| err.to_html(false))?;

    // Remove the user added interpretation and add a fake entry.
    // The original peptide will be added back later after downloading.
    // This is done because the PROXI servers do not necessarily support
    // the full ProForma spec and this decoy approach allows the user to
    // still use the full spec in the Annotator.
    let peptide = usi.interpretation.take();
    usi.interpretation = Some("A/1".to_string());

    let (backend, spectra) = usi
        .download_spectrum_async(None, None)
        .await
        .map_err(|e| match e {
            PROXIError::IO(backend, error) => BoxedError::new(
                BasicKind::Error,
                format!("Error while retrieving USI from {backend}"),
                error.to_string(),
                Context::none(),
            ),
            PROXIError::Error {
                backend,
                title,
                detail,
                ..
            } => BoxedError::new(
                BasicKind::Error,
                format!("Error while retrieving USI from {backend}: {title}"),
                detail,
                Context::none(),
            ),
            PROXIError::PeakUnavailable(backend, _) => BoxedError::new(
                BasicKind::Error,
                format!("Error while retrieving USI from {backend}"),
                "The peak data is not available for this dataset",
                Context::none(),
            ),
            PROXIError::NotFound => BoxedError::new(
                BasicKind::Error,
                "Error while retrieving USI",
                "No PROXI server responded to the request",
                Context::none(),
            ),
        })
        .map_err(|err| err.to_html(false))?;

    let mut spectrum = spectra
        .into_iter()
        .find(|s| {
            s.status
                .is_none_or(|s| s == mzdata::io::proxi::Status::Readable)
        })
        .unwrap();

    if let Ok(mut state) = state.lock() {
        // Fix the stored sequence after the decoy shuffling
        if let Some(peptide) = peptide.clone() {
            let param = PROXIParam::new(
                CURIE::new(mzdata::params::ControlledVocabulary::Unknown, 0),
                "sequence",
                Value::String(peptide),
            );
            if let Some(index) = spectrum
                .attributes
                .iter()
                .position(|p| p.name.eq_ignore_ascii_case("sequence"))
            {
                spectrum.attributes[index] = param;
            } else {
                spectrum.attributes.push(param);
            }
        } else {
            spectrum
                .attributes
                .retain(|p| !p.name.eq_ignore_ascii_case("sequence"));
        }
        let mut spectrum: MultiLayerSpectrum = spectrum.into();
        spectrum.description.id = format!(
            "{}:{}:{}",
            usi.dataset,
            usi.run_name,
            usi.identifier.map_or("-".to_string(), |i| match i {
                Identifier::Scan(s) => format!("scan:{s}"),
                Identifier::Index(s) => format!("index:{s}"),
                Identifier::NativeID(s) => format!("nativeId:{}", s.iter().join(",")),
            })
        );
        let mode = spectrum.precursor().and_then(|p| {
            get_mzdata_model(p.activation.methods(), &state.custom_models)
                .map(|(i, _)| i)
                .ok()
                .filter(|i| *i != crate::model::NONE)
        });

        state.spectra.push(RawFile::new_single(
            spectrum,
            format!("{backend} {} {}", usi.dataset, usi.run_name),
        ));

        Ok(IdentifiedPeptideSettings {
            peptide: peptide.unwrap_or_default(),
            charge: None,
            mode,
            warning: None,
        })
    } else {
        Err(BoxedError::small(
            BasicKind::Error,
            "Could not lock mutext",
            "Are you doing too many things at once?",
        )
        .to_html(false))
    }
}

#[tauri::command]
pub async fn load_raw<'a>(
    path: &'a str,
    state: ModifiableState<'a>,
) -> Result<RawFileDetails, String> {
    match mzdata::io::MZReaderType::open_path(path) {
        Ok(mut file) => {
            let file = if file.len() == 1 {
                let spec = file.next().unwrap();
                RawFile::new_single(dbg!(spec), path.to_string())
            } else {
                RawFile::new_file(path, file)
            };
            let spectra = &mut state.lock().unwrap().spectra;
            let details = file.details();
            spectra.push(file);
            Ok(details)
        }
        Err(err) => {
            let error = err.to_string();

            if err.kind() == ErrorKind::Other
                && (error == "It was not possible to find a compatible framework version."
                    || error
                        == "Feature which requires certain version of the hosting layer binaries was used on a version which doesn't support it."
                    || error == "One of the dependent libraries is missing.")
            {
                Err("The .NET 8.0 Runtime is needed to open Thermo RAW files. <a target='_blank' href='https://dotnet.microsoft.com/en-us/download/dotnet/8.0'>Which can be downloaded here.</a> Additionally on windows you can use <code style='user-select:all'>winget install Microsoft.DotNet.Runtime.8</code> for a quick install.".to_string())
            } else {
                Err(error)
            }
        }
    }
}

#[tauri::command]
pub async fn load_clipboard<'a>(data: &'a str, state: ModifiableState<'a>) -> Result<(), String> {
    let lines = data.lines().collect_vec();
    if data.is_empty() {
        return Err("Empty clipboard".to_string());
    }
    let (title, spectrum) = match lines[0].trim() {
        "#	m/z	Res.	S/N	I	I %	FWHM" => load_bruker_clipboard(&lines)
            .map(MultiLayerSpectrum::from_spectrum_like)
            .map(|s| ("Bruker", s)),

        "m/z	Charge	Intensity	FragmentType	MassShift	Position" => load_stitch_clipboard(&lines)
            .map(MultiLayerSpectrum::from_spectrum_like)
            .map(|s| ("Stitch", s)),

        "Mass/Charge Intensity" => load_sciex_clipboard(&lines)
            .map(MultiLayerSpectrum::from_spectrum_like)
            .map(|s| ("Sciex", s)),

        "SPECTRUM - MS" => load_thermo_clipboard(&lines)
            .map(MultiLayerSpectrum::from_spectrum_like)
            .map(|s| ("Thermo", s)),

        _ => Err("Not a recognised format (Bruker/Stitch/Sciex/Thermo)".to_string()),
    }?;

    state
        .lock()
        .unwrap()
        .spectra
        .push(RawFile::new_single(spectrum, format!("Clipboard {title}")));

    Ok(())
}

fn load_bruker_clipboard(lines: &[&str]) -> Result<mzdata::spectrum::CentroidSpectrum, String> {
    let mut spectrum = mzdata::spectrum::CentroidSpectrum::default();
    spectrum.description.signal_continuity = SignalContinuity::Centroid;
    let mut index = 0;

    for (line_number, line) in lines.iter().enumerate().skip(1) {
        let line_number = line_number + 1; // Humans like 1 based counting...
        let cells = line.split('\t').collect_vec();

        if cells.len() != 8 {
            return Err(format!("Incorrect number of columns at line {line_number}"));
        }
        if let (Ok(mass_over_charge), Ok(intensity)) =
            (cells[1].parse::<f64>(), cells[4].parse::<f32>())
        {
            spectrum.peaks.peaks.push(CentroidPeak {
                mz: mass_over_charge,
                intensity,
                index,
            });
            index += 1;
        } else {
            return Err(format!(
                "Could not read numbers at line {line_number} '{line}' '{}' '{}'",
                cells[1], cells[4]
            ));
        }
    }
    Ok(spectrum)
}

fn load_stitch_clipboard(lines: &[&str]) -> Result<mzdata::spectrum::CentroidSpectrum, String> {
    let mut spectrum = mzdata::spectrum::CentroidSpectrum::default();
    spectrum.description.signal_continuity = SignalContinuity::Centroid;
    let mut index = 0;

    for (line_number, line) in lines.iter().enumerate().skip(1) {
        let line_number = line_number + 1; // Humans like 1 based counting...
        let cells = line.split('\t').collect_vec();
        if cells.len() != 6 {
            return Err(format!("Incorrect number of columns at line {line_number}"));
        }
        if let (Ok(mass_over_charge), Ok(intensity)) =
            (cells[0].parse::<f64>(), cells[2].parse::<f32>())
        {
            spectrum.peaks.peaks.push(CentroidPeak {
                mz: mass_over_charge,
                intensity,
                index,
            });
            index += 1;
        } else {
            return Err(format!(
                "Could not read numbers at line {line_number} '{line}' '{}' '{}'",
                cells[0], cells[2]
            ));
        }
    }
    Ok(spectrum)
}

fn load_sciex_clipboard(lines: &[&str]) -> Result<mzdata::spectrum::CentroidSpectrum, String> {
    let mut spectrum = mzdata::spectrum::CentroidSpectrum::default();
    spectrum.description.signal_continuity = SignalContinuity::Centroid;
    let mut index = 0;

    for (line_number, line) in lines.iter().enumerate().skip(1) {
        let line_number = line_number + 1; // Humans like 1 based counting...
        let cells = line.split(' ').collect_vec();
        if cells.len() != 2 {
            return Err(format!("Incorrect number of columns at line {line_number}"));
        }
        if let (Ok(mass_over_charge), Ok(intensity)) =
            (cells[0].parse::<f64>(), cells[1].parse::<f32>())
        {
            spectrum.peaks.peaks.push(CentroidPeak {
                mz: mass_over_charge,
                intensity,
                index,
            });
            index += 1;
        } else {
            return Err(format!(
                "Could not read numbers at line {line_number} '{line}' '{}' '{}'",
                cells[0], cells[1]
            ));
        }
    }
    Ok(spectrum)
}

fn load_thermo_clipboard(lines: &[&str]) -> Result<mzdata::spectrum::RawSpectrum, String> {
    let mut spectrum = mzdata::spectrum::RawSpectrum::default();
    spectrum.description.signal_continuity = SignalContinuity::Profile;
    let mut data_mz = DataArray::default();
    data_mz.name = ArrayType::MZArray;
    data_mz.dtype = BinaryDataArrayType::Float64;
    data_mz.compression = BinaryCompressionType::Decoded;
    let mut data_intensity = DataArray::default();
    data_intensity.name = ArrayType::IntensityArray;
    data_intensity.dtype = BinaryDataArrayType::Float32;
    data_intensity.compression = BinaryCompressionType::Decoded;

    for (line_number, line) in lines.iter().enumerate().skip(1) {
        if line.eq_ignore_ascii_case("spectrum - ms") || line.eq_ignore_ascii_case("mass intensity")
        {
            continue;
        }
        let line_number = line_number + 1; // Humans like 1 based counting...
        if line_number == 2 {
            spectrum.description.id = line.to_string();
        }
        let colon_cells = line.split(':').collect_vec();
        if colon_cells.len() == 2 {
            match (colon_cells[0], colon_cells[1]) {
                ("Scan #", scan) => {
                    match scan.trim().split_once("-") {
                        Some((first, end)) => {
                            spectrum.description.params.add_param(Param {
                                name: "Scan number start".to_string(),
                                value: Value::Int(first.parse().map_err(|err| {
                                    format!("Could not read first scan number: {err}")
                                })?),
                                accession: None,
                                controlled_vocabulary: None,
                                unit: Unit::Unknown,
                            });
                            spectrum.description.params.add_param(Param {
                                name: "Scan number end".to_string(),
                                value: Value::Int(end.parse().map_err(|err| {
                                    format!("Could not read end scan number: {err}")
                                })?),
                                accession: None,
                                controlled_vocabulary: None,
                                unit: Unit::Unknown,
                            });
                        }
                        None => spectrum.description.params.add_param(Param {
                            name: "Scan number".to_string(),
                            value: Value::Int(
                                scan.trim()
                                    .parse()
                                    .map_err(|err| format!("Could not read scan number: {err}"))?,
                            ),
                            accession: None,
                            controlled_vocabulary: None,
                            unit: Unit::Unknown,
                        }),
                    }
                }
                ("RT", rt) => match rt.trim().split_once("-") {
                    Some((first, end)) => {
                        spectrum.description.params.add_param(Param {
                            name: "Retention time start".to_string(),
                            value: Value::Float(first.parse().map_err(|err| {
                                format!("Could not read first retention time: {err}")
                            })?),
                            accession: None,
                            controlled_vocabulary: None,
                            unit: Unit::Minute,
                        });
                        spectrum.description.params.add_param(Param {
                            name: "Retention time end".to_string(),
                            value: Value::Float(end.parse().map_err(|err| {
                                format!("Could not read end retention time: {err}")
                            })?),
                            accession: None,
                            controlled_vocabulary: None,
                            unit: Unit::Minute,
                        });
                    }
                    None => spectrum.description.params.add_param(Param {
                        name: "Retention time".to_string(),
                        value: Value::Float(
                            rt.trim()
                                .parse()
                                .map_err(|err| format!("Could not read retention time: {err}"))?,
                        ),
                        accession: None,
                        controlled_vocabulary: None,
                        unit: Unit::Minute,
                    }),
                },
                _ => (),
            }
        }
        let cells = line.split('\t').collect_vec();
        if cells.len() == 1 {
            continue;
        }
        if cells.len() != 2 {
            return Err(format!("Incorrect number of columns at line {line_number}"));
        }
        if let (Ok(mass_over_charge), Ok(intensity)) =
            (cells[0].parse::<f64>(), cells[1].parse::<f32>())
        {
            let _ = data_mz.push(mass_over_charge);
            let _ = data_intensity.push(intensity);
        } else {
            return Err(format!(
                "Could not read numbers at line {line_number} '{line}' '{}' '{}'",
                cells[0], cells[1]
            ));
        }
    }
    spectrum.arrays.add(data_mz);
    spectrum.arrays.add(data_intensity);
    Ok(spectrum)
}

#[tauri::command]
pub fn select_spectrum_index(
    file_index: usize,
    index: usize,
    state: ModifiableState,
) -> Result<(), &'static str> {
    state
        .lock()
        .unwrap()
        .spectra
        .iter_mut()
        .find(|f| f.id() == file_index)
        .ok_or("File index not valid")
        .and_then(|file| file.select_index(index))
}

#[tauri::command]
pub fn select_spectrum_native_id(
    file_index: usize,
    native_id: String,
    state: ModifiableState,
) -> Result<(), &'static str> {
    state
        .lock()
        .unwrap()
        .spectra
        .iter_mut()
        .find(|f| f.id() == file_index)
        .ok_or("File index not valid")
        .and_then(|file| file.select_native_id(native_id))
}

#[tauri::command]
pub fn select_retention_time(
    file_index: usize,
    rt: RangeInclusive<OrderedTime>,
    state: ModifiableState,
) -> Result<(), &'static str> {
    state
        .lock()
        .unwrap()
        .spectra
        .iter_mut()
        .find(|f| f.id() == file_index)
        .ok_or("File index not valid")
        .and_then(|file| file.select_retention_time(rt))
}

#[tauri::command]
pub fn close_raw_file(file_index: usize, state: ModifiableState) {
    let _ = state.lock().map(|mut state| {
        state
            .spectra
            .iter_mut()
            .position(|f| f.id() == file_index)
            .map(|index| state.spectra.remove(index))
    });
}

#[tauri::command]
pub fn deselect_spectrum(file_index: usize, index: usize, state: ModifiableState) {
    let _ = state.lock().map(|mut state| {
        state
            .spectra
            .iter_mut()
            .find(|f| f.id() == file_index)
            .map(|file| file.unselect_index(index))
    });
}

#[tauri::command]
pub fn get_open_raw_files(state: ModifiableState) -> Vec<RawFileDetails> {
    state
        .lock()
        .unwrap()
        .spectra
        .iter()
        .map(|file| file.details())
        .collect_vec()
}

#[tauri::command]
pub fn get_selected_spectra(
    state: ModifiableState,
) -> Vec<(usize, bool, Vec<SelectedSpectrumDetails>)> {
    if let Ok(mut state) = state.lock() {
        let (spectra, models) = state.spectra_and_models();
        spectra
            .iter_mut()
            .map(|file| {
                let id = file.id();
                (
                    id,
                    matches!(file, RawFile::Single { .. }),
                    file.get_selected_spectra()
                        .map(|spectrum| SelectedSpectrumDetails {
                            id: spectrum.index(),
                            short: spectrum.description.id.to_string(),
                            description: spectrum_description(
                                spectrum.description(),
                                &spectrum.peaks().fetch_summaries(),
                                models,
                            ),
                        })
                        .collect_vec(),
                )
            })
            .collect_vec()
    } else {
        Vec::new()
    }
}

pub fn spectrum_description(
    description: &SpectrumDescription,
    summary: &SpectrumSummary,
    models: &[(String, FragmentationModel)],
) -> String {
    format!(
        "index: {} id: {}<br>time: {:.3} min signal mode: {:?} ms level: {} ion mobility: {}<br>mz range: {:.1} — {:.1} peak count: {} tic: {:.3e} base peak intensity: {:.3e} resultion: {}<br>{}<br>{}{}",
        description.index,
        description.id,
        description.acquisition.start_time() / 60.0,
        description.signal_continuity,
        description.ms_level,
        description
            .acquisition
            .scans
            .first()
            .and_then(|s| s.ion_mobility())
            .to_optional_string(),
        summary.mz_range.0,
        summary.mz_range.1,
        summary.count,
        summary.tic,
        summary.base_peak.intensity,
        description
            .acquisition
            .first_scan()
            .and_then(|s| s.resolution())
            .map_or("-".to_string(), |v| match v {
                ValueRef::Float(f) => format!("{f:.3e}"),
                ValueRef::Int(i) => format!("{i:.3e}"),
                _ => "-".to_string(),
            }),
        description
            .acquisition
            .first_scan()
            .and_then(|s| s.filter_string())
            .map_or("No filter string".to_string(), |v| format!("Filter: {v}")),
        description.precursor.first().map_or("No precursor".to_string(), |p| {
            let i = p.isolation_window();
            format!(
                "Precursor mass: {} charge: {} target: {} range: {} — {} method: {} energy: {:.1}",
                display_mass(Mass::new::<dalton>(p.neutral_mass()), None),
                p.charge().map_or("-".to_string(), |v| format!("{v:+.0}")),
                i.target,
                i.lower_bound,
                i.upper_bound,
                get_mzdata_model(p.activation.methods(), models).map_or_else(|n| n, |(_, n)| n),
                p.activation.energy,
            )
        }),
        if let Some(param) = &description.params().iter().find(|p| p.name == "sequence") {
            format!(
                "<br>Sequence: <span style='-webkit-user-select:all;user-select:all;'>{}</span>",
                param.value
            )
        } else {
            String::new()
        },
    )
}

pub fn create_selected_spectrum(
    state: &mut crate::State,
    filter: f32,
) -> Result<MultiLayerSpectrum, BoxedError<'_, BasicKind>> {
    let mut spectra = Vec::new();
    for file in state.spectra.iter_mut() {
        spectra.extend(file.get_selected_spectra());
    }
    let spectrum = if spectra.is_empty() {
        return Err(BoxedError::new(
            BasicKind::Error,
            "No selected spectra",
            "Select a spectrum from an open raw file, or open a raw file if none are opened yet",
            Context::none(),
        ));
    } else if spectra.len() == 1 {
        let mut spectrum = spectra.pop().unwrap();
        if filter != 0.0 && spectrum.arrays.is_some() {
            // TODO: When centroided data is loaded denoising is not applied
            spectrum.denoise(filter).map_err(|err| {
                BoxedError::new(
                    BasicKind::Error,
                    "Spectrum could not be denoised",
                    err.to_string(),
                    Context::none(),
                )
            })?;
        }

        if spectrum.signal_continuity() == SignalContinuity::Profile {
            spectrum
                .pick_peaks_with(&PeakPicker::default())
                .map_err(|err| {
                    BoxedError::new(
                        BasicKind::Error,
                        "Spectrum could not be peak picked",
                        err.to_string(),
                        Context::none(),
                    )
                })?;
            if let Some(p) = spectrum.peaks.as_mut() {
                p.peaks.retain(|p| p.intensity > 0.1)
            } // Stupid filter to remove very low peaks
            spectrum.description.signal_continuity = SignalContinuity::Centroid; // Not done by the above function
            dbg!(&spectrum);
        } else if spectrum.arrays.is_some() && spectrum.peaks.is_none() {
            // USI spectra are mostly loaded as the binary array maps instead of peaks regardless of the signal continuity level
            spectrum.peaks = spectrum.arrays.as_ref().map(|a| a.into());
        }
        spectrum
    } else {
        let description = spectra
            .first()
            .map(|s| s.description.clone())
            .unwrap_or_default();
        let averaged = mzdata::spectrum::average_spectra(&spectra, 0.001);
        let mut binding = averaged.intensity_array().to_owned();
        let picked = if filter != 0.0 {
            let denoised =
                mzsignal::denoise::denoise(&averaged.mz_array, &mut binding, filter).unwrap();
            mzsignal::peak_picker::pick_peaks(&averaged.mz_array, denoised).unwrap()
        } else {
            mzsignal::peak_picker::pick_peaks(&averaged.mz_array, &averaged.intensity_array)
                .unwrap()
        };

        MultiLayerSpectrum::new(
            description,
            None,
            Some(PeakSetVec::new(
                picked.into_iter().map(Into::into).collect(),
            )),
            None,
        )
    };

    Ok(spectrum)
}

#[tauri::command]
pub fn save_spectrum<'a>(
    state: ModifiableState,
    filter: f32,
    path: &Path,
    sequence: &'a str,
    model: &'a str,
    charge: Option<usize>,
) -> Result<(), String> {
    let state = state.lock().map_err(|e| {
        BoxedError::new(
            BasicKind::Error,
            "Could not write file",
            "Mutex locked, are you doing too much at the same time?",
            Context::show(e.to_string()),
        )
        .to_html(false)
    })?;
    let file = std::fs::OpenOptions::new()
        .create(true)
        .truncate(true)
        .write(true)
        .open(path)
        .map_err(|e| {
            BoxedError::new(
                BasicKind::Error,
                "Could not write file",
                "File could not be opened for writing",
                Context::show(e.to_string()),
            )
            .to_html(false)
        })?;
    let ext = path
        .extension()
        .map(|s| s.to_string_lossy().to_ascii_lowercase());
    let mut spectrum = state.annotated_spectrum.clone().ok_or_else(|| BoxedError::small(
                BasicKind::Error,
                "No spectrum present",
                "No spectrum was present so no spectrum can be saved. Annotate a spectrum or load one from a library to save a spectrum.",
            )
            .to_html(false))?;
    if !sequence.is_empty() {
        spectrum.description.params.push(Param::new_key_value(
            "sequence",
            mzdata::params::Value::String(sequence.to_string()),
        ));
    }
    if model != "all" && model != "none" && model != "custom" {
        if let Some(p) = spectrum.description.precursor.first_mut() {
            p.activation.methods_mut().append(&mut match model {
                "ethcd" => vec![
                    DissociationMethodTerm::ElectronTransferDissociation,
                    DissociationMethodTerm::CollisionInducedDissociation,
                ],
                "hot_eacid" => vec![
                    DissociationMethodTerm::ElectronActivatedDissociation,
                    DissociationMethodTerm::CollisionInducedDissociation,
                ],
                "ead" => vec![DissociationMethodTerm::ElectronActivatedDissociation],
                "cidhcd" => vec![DissociationMethodTerm::CollisionInducedDissociation],
                "etd" => vec![DissociationMethodTerm::ElectronTransferDissociation],
                "td_etd" => vec![DissociationMethodTerm::ElectronTransferDissociation],
                _ => Vec::new(),
            })
        }
    }
    if let Some(charge) = charge {
        if let Some(p) = spectrum.description.precursor.first_mut() {
            if let Some(i) = p.ions.first_mut() {
                i.charge = Some(charge as i32);
            }
        }
    }
    match ext.as_deref() {
        Some("mgf") => {
            let mut writer: mzdata::io::mgf::MGFWriterType<
                std::fs::File,
                mzannotate::prelude::AnnotatedPeak<mzannotate::prelude::Fragment>,
                mzpeaks::DeconvolutedPeak,
                mzdata::io::mgf::MZDataMGFStyle,
            > = mzdata::io::mgf::MGFWriterType::new(file);
            writer.write(&spectrum).map(|_| ()).map_err(|e| {
                BoxedError::new(
                    BasicKind::Error,
                    "Could not write file",
                    "Could not write MGF",
                    Context::show(e.to_string()),
                )
                .to_html(false)
            })
        }
        Some("mzml") => {
            let mut writer: mzdata::io::mzml::MzMLWriterType<
                std::fs::File,
                mzannotate::prelude::AnnotatedPeak<mzannotate::prelude::Fragment>,
                mzpeaks::DeconvolutedPeak,
            > = mzdata::io::mzml::MzMLWriterType::new(file);
            writer.write_spectrum(&spectrum).map_err(|e| {
                BoxedError::new(
                    BasicKind::Error,
                    "Could not write file",
                    "Could not write mzML",
                    Context::show(e.to_string()),
                )
                .to_html(false)
            })
        }
        Some("txt") => {
            let mut writer =
                mzannotate::mzspeclib::MzSpecLibTextWriter::new(std::io::BufWriter::new(file));
            writer.write_header().map_err(|e| {
                BoxedError::new(
                    BasicKind::Error,
                    "Could not write file",
                    "Could not write mzSpecLib",
                    Context::show(e.to_string()),
                )
                .to_html(false)
            })?;
            writer.write_spectrum(&spectrum).map_err(|e| {
                BoxedError::new(
                    BasicKind::Error,
                    "Could not write file",
                    "Could not write mzSpecLib",
                    Context::show(e.to_string()),
                )
                .to_html(false)
            })
        }
        _ => Err(BoxedError::new(
            BasicKind::Error,
            "Could not write file",
            "Invalid path, use mgf, mzml, or mzspeclib.txt as extension",
            Context::show(path.to_string_lossy()).to_owned(),
        )
        .to_html(false)),
    }
}
