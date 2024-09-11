use std::{any::Any, io::ErrorKind};

use itertools::Itertools;
use mzdata::{
    io::MZFileReader,
    params::{ParamDescribed, Unit, Value},
    prelude::{IonProperties, PrecursorSelection, SpectrumLike},
    spectrum::{
        bindata::BinaryCompressionType, ActivationMethod, ArrayType, BinaryDataArrayType,
        DataArray, MultiLayerSpectrum, SignalContinuity,
    },
    Param,
};
use mzpeaks::CentroidPeak;
use rustyms::system::{dalton, Mass};
use serde::{Deserialize, Serialize};

use crate::{
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
pub async fn load_raw<'a>(
    path: &'a str,
    state: ModifiableState<'a>,
) -> Result<RawFileDetails, String> {
    match mzdata::io::MZReaderType::open_path(path) {
        Ok(file) => {
            let file = RawFile::new_file(path, file);
            let spectra = &mut state.lock().unwrap().spectra;
            let details = file.details();
            spectra.push(file);
            Ok(details)
        }
        Err(err) => {
            let error = err.to_string();

            if err.kind() == ErrorKind::Other
                && error == "It was not possible to find a compatible framework version."
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
    let spectrum = match lines[0].trim() {
        "#	m/z	Res.	S/N	I	I %	FWHM" => {
            load_bruker_clipboard(&lines).map(MultiLayerSpectrum::from_spectrum_like)
        }
        "m/z	Charge	Intensity	FragmentType	MassShift	Position" => {
            load_stitch_clipboard(&lines).map(MultiLayerSpectrum::from_spectrum_like)
        }
        "Mass/Charge Intensity" => {
            load_sciex_clipboard(&lines).map(MultiLayerSpectrum::from_spectrum_like)
        }
        "SPECTRUM - MS" => {
            load_thermo_clipboard(&lines).map(MultiLayerSpectrum::from_spectrum_like)
        }
        _ => Err("Not a recognised format (Bruker/Stitch/Sciex/Thermo)".to_string()),
    }?;

    state
        .lock()
        .unwrap()
        .spectra
        .push(RawFile::new_clipboard(spectrum));

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
        if line.to_lowercase() == "spectrum - ms" || line.to_ascii_lowercase() == "mass	intensity" {
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
pub fn unselect_spectrum(file_index: usize, index: usize, state: ModifiableState) {
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
                SelectedSpectrumDetails {
                    id,
                    short: d.id.to_string(),
                    description: format!(
                        "{}<br>time: {} signal mode: {:?} ms level: {} ion mobility: {}<br>{}{}",
                        d.id,
                        spectrum.start_time(),
                        spectrum.signal_continuity(),
                        spectrum.ms_level(),
                        spectrum.ion_mobility().to_optional_string(),
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
                                        ActivationMethod::BeamTypeCollisionInducedDissociation => "BeamCID",
                                        ActivationMethod::CollisionInducedDissociation => "CID",
                                        ActivationMethod::ElectronActivationDissociation => "EAD",
                                        ActivationMethod::ElectronCaptureDissociation => "ECD",
                                        ActivationMethod::ElectronTransferDissociation => "ETD",
                                        ActivationMethod::HighEnergyCollisionInducedDissociation => "HCD",
                                        ActivationMethod::InSourceCollisionInducedDissociation => "isCID",
                                        ActivationMethod::LowEnergyCollisionInducedDissociation => "lowCID",
                                        ActivationMethod::NegativeElectronTransferDissociation => "nETD",
                                        ActivationMethod::Other(_) => "other",
                                        ActivationMethod::Photodissociation => "PD",
                                        ActivationMethod::SupplementalBeamTypeCollisionInducedDissociation => "sBeamCID",
                                        ActivationMethod::SupplementalCollisionInducedDissociation => "sCID",
                                        ActivationMethod::TrapTypeCollisionInducedDissociation => "trapCID",
                                        ActivationMethod::UltravioletPhotodissociation => "UVPD",
                                    }),
                                p.activation.energy,
                            )
                        }),
                        if let Some(param) = SpectrumLike::params(&spectrum).iter().find(|p| p.name == "sequence") {
                            format!("<br>Sequence: <span style='user-select:all'>{}</span>", param.value)
                        } else {
                            String::new()
                        },
                    )
                }
            }).collect_vec())
        }).collect_vec()
}
