use std::{io::ErrorKind, ops::RangeInclusive, str::FromStr};

use itertools::Itertools;
use mzdata::{
    io::{proxi::PROXIError, usi::USIParseError, MZFileReader, SpectrumSource},
    meta::DissociationMethodTerm,
    params::{ParamDescribed, Unit, Value},
    prelude::{IonProperties, PrecursorSelection, SpectrumLike},
    spectrum::{
        bindata::BinaryCompressionType, ArrayType, BinaryDataArrayType, DataArray,
        MultiLayerSpectrum, SignalContinuity,
    },
    Param,
};
use mzpeaks::CentroidPeak;
use rustyms::{
    error::{Context, CustomError},
    system::{dalton, Mass, OrderedTime},
};
use serde::{Deserialize, Serialize};

use crate::{
    identified_peptides::IdentifiedPeptideSettings,
    metadata_render::OptionalString,
    render::display_mass,
    state::{RawFile, RawFileDetails},
    ModifiableState,
};

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct SelectedSpectrumDetails {
    id: usize,
    short: String,
    description: String,
}

#[tauri::command]
pub async fn load_usi<'a>(
    usi: &'a str,
    state: ModifiableState<'a>,
) -> Result<IdentifiedPeptideSettings, CustomError> {
    let usi = mzdata::io::usi::USI::from_str(usi).map_err(|e| match e {
        USIParseError::MalformedIndex(index, e, full) => CustomError::error(
            "Invalid USI: Invalid index",
            format!("The index is not valid: {e}"),
            Context::line(
                None,
                &full,
                full.find(&index).unwrap_or_default(),
                index.len(),
            ),
        ),
        USIParseError::MissingDataset(full) => CustomError::error(
            "Invalid USI: Missing dataset",
            "The dataset section of the USI is missing",
            Context::show(full),
        ),
        USIParseError::MissingRun(full) => CustomError::error(
            "Invalid USI: Missing run",
            "The run section of the USI is missing",
            Context::show(full),
        ),
        USIParseError::UnknownIndexType(index, full) => CustomError::error(
            "Invalid USI: Unknown index type",
            "The index has to be one of 'index', 'scan', or 'nativeId'",
            Context::line(
                None,
                &full,
                full.find(&index).unwrap_or_default(),
                index.len(),
            ),
        ),
        USIParseError::UnknownProtocol(protocol, full) => CustomError::error(
            "Invalid USI: Unknown protocol",
            "The protocol has to be 'mzspec'",
            Context::line(
                None,
                &full,
                full.find(&protocol).unwrap_or_default(),
                protocol.len(),
            ),
        ),
    })?;

    let (backend, spectra) =
        usi.download_spectrum_async(None, None)
            .await
            .map_err(|e| match e {
                PROXIError::IO(backend, error) => CustomError::error(
                    format!("Error while retrieving USI from {backend}"),
                    error,
                    Context::None,
                ),
                PROXIError::Error {
                    backend,
                    title,
                    detail,
                    ..
                } => CustomError::error(
                    format!("Error while retrieving USI from {backend}: {title}"),
                    detail,
                    Context::None,
                ),
                PROXIError::PeakUnavailable(backend, _) => CustomError::error(
                    format!("Error while retrieving USI from {backend}"),
                    "The peak data is not available for this dataset",
                    Context::None,
                ),
                PROXIError::NotFound => CustomError::error(
                    "Error while retrieving USI",
                    "No PROXI server responded to the request",
                    Context::None,
                ),
            })?;

    let spectrum = spectra
        .into_iter()
        .find(|s| {
            s.status
                .is_none_or(|s| s == mzdata::io::proxi::Status::Readable)
        })
        .unwrap();

    if let Ok(mut state) = state.lock() {
        state.spectra.push(RawFile::new_single(
            spectrum.into(),
            format!("{backend} {} {}", usi.dataset, usi.run_name),
        ));

        Ok(IdentifiedPeptideSettings {
            peptide: usi.interpretation.unwrap_or_default(),
            charge: None,
            mode: None,
            warning: None,
        })
    } else {
        Err(CustomError::error(
            "Could not lock mutext",
            "Are you doing too many things at once?",
            Context::None,
        ))
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
                RawFile::new_single(file.next().unwrap(), path.to_string())
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
                    || error == "Feature which requires certain version of the hosting layer binaries was used on a version which doesn't support it." 
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
pub fn get_selected_spectra(state: ModifiableState) -> Vec<(usize, Vec<SelectedSpectrumDetails>)> {
    state
        .lock()
        .unwrap()
        .spectra
        .iter_mut()
        .map(|file| {
            let id = file.id();
            (id,
            file.get_selected_spectra().map(|spectrum | {
                let d = spectrum.description();
                let p = spectrum.precursor();
                let summary = spectrum.peaks().fetch_summaries();
                SelectedSpectrumDetails {
                    id: spectrum.index(),
                    short: d.id.to_string(),
                    description: format!(
                        "index: {} id: {}<br>time: {:.3} signal mode: {:?} ms level: {} ion mobility: {}<br>mz range: {:.1} — {:.1} peak count: {} tic: {:.3e} base peak intensity: {:.3e}<br>{}{}",
                        spectrum.index(),
                        d.id,
                        spectrum.start_time(),
                        spectrum.signal_continuity(),
                        spectrum.ms_level(),
                        spectrum.ion_mobility().to_optional_string(),
                        summary.mz_range.0,
                        summary.mz_range.1,
                        summary.count,
                        summary.tic,
                        summary.base_peak.intensity,
                        p.map_or("No precursor".to_string(), |p| {
                            let i = p.isolation_window();
                            format!(
                                "Precursor mass: {} charge: {} target: {} range: {} — {} method: {} energy: {:.1}",
                                display_mass(Mass::new::<dalton>(p.neutral_mass()), None),
                                p.charge().map_or("-".to_string(), |v| format!("{v:+.0}")),
                                i.target,
                                i.lower_bound,
                                i.upper_bound,
                                p.activation
                                    .method()
                                    .map_or("-", |m| match m {
                                        DissociationMethodTerm::BeamTypeCollisionInducedDissociation => "BeamCID",
                                        DissociationMethodTerm::BlackbodyInfraredRadiativeDissociation => "BlackbodyInfraredRadiativeDissociation",
                                        DissociationMethodTerm::CollisionInducedDissociation => "CID",
                                        DissociationMethodTerm::DissociationMethod => "Dissociation",
                                        DissociationMethodTerm::ElectronActivatedDissociation => "EAD",
                                        DissociationMethodTerm::ElectronCaptureDissociation => "ECD",
                                        DissociationMethodTerm::ElectronTransferDissociation => "ETD",
                                        DissociationMethodTerm::HigherEnergyBeamTypeCollisionInducedDissociation => "BeamHCD",
                                        DissociationMethodTerm::InfraredMultiphotonDissociation => "InfraredMultiphotonDissociation",
                                        DissociationMethodTerm::InSourceCollisionInducedDissociation => "isCID",
                                        DissociationMethodTerm::LIFT => "LIFT",
                                        DissociationMethodTerm::LowEnergyCollisionInducedDissociation => "lowCID",
                                        DissociationMethodTerm::NegativeElectronTransferDissociation => "nETD",
                                        DissociationMethodTerm::Photodissociation => "PD",
                                        DissociationMethodTerm::PlasmaDesorption => "PlasmaDesorption",
                                        DissociationMethodTerm::PostSourceDecay => "PostSourceDecay",
                                        DissociationMethodTerm::PulsedQDissociation => "PulsedDissociation",
                                        DissociationMethodTerm::SupplementalBeamTypeCollisionInducedDissociation => "sBeamCID",
                                        DissociationMethodTerm::SupplementalCollisionInducedDissociation => "sCID",
                                        DissociationMethodTerm::SurfaceInducedDissociation => "SurfaceInducedDissocation",
                                        DissociationMethodTerm::SustainedOffResonanceIrradiation => "SustainedOffResonanceIrradiation",
                                        DissociationMethodTerm::TrapTypeCollisionInducedDissociation => "trapCID",
                                        DissociationMethodTerm::UltravioletPhotodissociation => "UVPD",
                                    }),
                                p.activation.energy,
                            )
                        }),
                        if let Some(param) = &spectrum.params().iter().find(|p| p.name == "sequence") {
                            format!("<br>Sequence: <span style='-webkit-user-select:all;user-select:all;'>{}</span>", param.value)
                        } else {
                            String::new()
                        },
                    )
                }
            }).collect_vec())
        }).collect_vec()
}
