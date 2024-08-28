use itertools::Itertools;
use mzdata::{io::{MZFileReader, SpectrumSource}, prelude::{IonProperties, PrecursorSelection, SpectrumLike}, spectrum::ActivationMethod};
use rustyms::{spectrum::RawPeak, system::{dalton, mz, Mass, MassOverCharge}, RawSpectrum};
use serde::{Deserialize, Serialize};


use crate::{metadata_render::OptionalString, render::display_mass, state::{RawFile, RawFileDetails}, ModifiableState};

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct SelectedSpectrumDetails {
    id: usize, 
    short: String,
    description: String,
}

#[tauri::command]
pub async fn load_raw<'a>(path: &'a str, state: ModifiableState<'a>) -> Result<RawFileDetails, String> {
    match mzdata::io::MZReaderType::open_path(path) {
        Ok(file) => {
            let file = RawFile::new(path, file);
            let spectra = &mut state.lock().unwrap().spectra;
            let details = file.details();
            spectra.push(file);
            Ok(details)
        }
        Err(err) => Err(err.to_string()),
    }
}


#[tauri::command]
pub fn load_clipboard(data: &str, state: ModifiableState) -> Result<usize, String> {
    let lines = data.lines().collect_vec();
    if data.is_empty() {
        return Err("Empty clipboard".to_string());
    }
    let spectrum = match lines[0].trim() {
        "#	m/z	Res.	S/N	I	I %	FWHM" => load_bruker_clipboard(&lines),
        "m/z	Charge	Intensity	FragmentType	MassShift	Position" => load_stitch_clipboard(&lines),
        "Mass/Charge Intensity" => load_sciex_clipboard(&lines),
        _ => Err("Not a recognised format (Bruker/Stitch/Sciex)".to_string()),
    }?;

    let mut new_spectrum = RawSpectrum::default();
    new_spectrum.extend(spectrum);
    new_spectrum.title = "Clipboard".to_string();
    new_spectrum.charge = None;

    // state.lock().unwrap().spectra = vec![new_spectrum]; // TODO: figure out how to store clipboard spectra
    Ok(1)
}

fn load_bruker_clipboard(lines: &[&str]) -> Result<Vec<RawPeak>, String> {
    let mut spectrum = Vec::new();
    for (line_number, line) in lines.iter().enumerate().skip(1) {
        let line_number = line_number + 1; // Humans like 1 based counting...
        let cells = line.split('\t').collect_vec();
        if cells.len() != 8 {
            return Err(format!("Incorrect number of columns at line {line_number}"));
        }
        if let (Ok(mass_over_charge), Ok(intensity)) = (cells[1].parse(), cells[4].parse()) {
            spectrum.push(RawPeak {
                mz: MassOverCharge::new::<mz>(mass_over_charge),
                intensity,
            })
        } else {
            return Err(format!(
                "Could not read numbers at line {line_number} '{line}' '{}' '{}'",
                cells[1], cells[4]
            ));
        }
    }
    Ok(spectrum)
}

fn load_stitch_clipboard(lines: &[&str]) -> Result<Vec<RawPeak>, String> {
    let mut spectrum = Vec::new();
    for (line_number, line) in lines.iter().enumerate().skip(1) {
        let line_number = line_number + 1; // Humans like 1 based counting...
        let cells = line.split('\t').collect_vec();
        if cells.len() != 6 {
            return Err(format!("Incorrect number of columns at line {line_number}"));
        }
        if let (Ok(mass_over_charge), Ok(intensity)) = (cells[0].parse(), cells[2].parse()) {
            spectrum.push(RawPeak {
                mz: MassOverCharge::new::<mz>(mass_over_charge),
                intensity,
            })
        } else {
            return Err(format!(
                "Could not read numbers at line {line_number} '{line}' '{}' '{}'",
                cells[0], cells[2]
            ));
        }
    }
    Ok(spectrum)
}

fn load_sciex_clipboard(lines: &[&str]) -> Result<Vec<RawPeak>, String> {
    let mut spectrum = Vec::new();
    for (line_number, line) in lines.iter().enumerate().skip(1) {
        let line_number = line_number + 1; // Humans like 1 based counting...
        let cells = line.split(' ').collect_vec();
        if cells.len() != 2 {
            return Err(format!("Incorrect number of columns at line {line_number}"));
        }
        if let (Ok(mass_over_charge), Ok(intensity)) = (cells[0].parse(), cells[1].parse()) {
            spectrum.push(RawPeak {
                mz: MassOverCharge::new::<mz>(mass_over_charge),
                intensity,
            })
        } else {
            return Err(format!(
                "Could not read numbers at line {line_number} '{line}' '{}' '{}'",
                cells[0], cells[1]
            ));
        }
    }
    Ok(spectrum)
}

#[tauri::command]
pub fn select_spectrum_index(file_index: usize, index: usize, state: ModifiableState) -> Result<(), &'static str> {
    state
    .lock()
    .unwrap()
    .spectra
    .iter_mut()
    .find(|f| f.id == file_index)
    .ok_or("File index not valid")
    .and_then(|file| 
        if index >= file.rawfile.len() {
            Err("Outside of file range")
        } else if !file.selected_spectra.contains(&index) {
            file.selected_spectra.push(index); 
            file.selected_spectra.sort();
            Ok(())
        } else {
            Ok(())
        })
}

#[tauri::command]
pub fn select_spectrum_native_id(file_index: usize, native_id: String, state: ModifiableState) -> Result<(), &'static str> {
    state
    .lock()
    .unwrap()
    .spectra
    .iter_mut()
    .find(|f| f.id == file_index)
    .ok_or("File index not valid")
    .and_then(|file|  
        file.rawfile.get_spectrum_by_id(&native_id).map(|s| {
            if !file.selected_spectra.contains(&s.index()) {
                file.selected_spectra.push(s.index()); 
                file.selected_spectra.sort();
            }
        }).ok_or("Scan number does not exist"))
}

#[tauri::command]
pub fn close_raw_file(file_index: usize, state: ModifiableState) {
    let _ = state
    .lock()
    .map(|mut state| 
        state.spectra.iter_mut()
        .position(|f| f.id == file_index).map(|index| state.spectra.remove(index))
    );
}

#[tauri::command]
pub fn unselect_spectrum(file_index: usize, index: usize, state: ModifiableState) {
    let _ = state
    .lock()
    .map(|mut state| 
        state.spectra.iter_mut()
        .find(|f| f.id == file_index).map(|file| {
            if let Some(index) = file.selected_spectra.iter().position(|i| *i == index) {
            file.selected_spectra.remove(index);
            }
        })
    );
}

#[tauri::command]
pub fn get_open_raw_files(state: ModifiableState) -> Vec<RawFileDetails> {
    state
    .lock()
    .unwrap()
    .spectra.iter().map(|file|file.details()).collect_vec()
}

#[tauri::command]
pub fn get_selected_spectra(state: ModifiableState) -> Vec<(usize, Vec<SelectedSpectrumDetails>)> {
    state
        .lock()
        .unwrap()
        .spectra
        .iter_mut()
        .map(|file| (file.id, 
            file.selected_spectra.iter().map(|index| {
                let spectrum =file.rawfile.get_spectrum_by_index(*index).expect("Spectrum index not valid");
                let d = spectrum.description();
                let p = spectrum.precursor();
                SelectedSpectrumDetails {
                    id: *index,
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
                        if let Some(param) = spectrum.params().iter().find(|p| p.name == "sequence") {
                            format!("<br>Sequence: <span style='user-select:all'>{}</span>", param.value)
                        } else {
                            String::new()
                        },                        
                    )
                }
            }).collect_vec())
        ).collect_vec()
}