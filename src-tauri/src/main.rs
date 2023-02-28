#![cfg_attr(
    all(not(debug_assertions), target_os = "windows"),
    windows_subsystem = "windows"
)]

use itertools::Itertools;
use mass_alignment::{template::Template, *};
use pdbtbx::*;
use rustyms::{Mass, MassOverCharge, Zero};
use std::collections::{HashMap, HashSet};
use std::fmt::Write;

// Learn more about Tauri commands at https://tauri.app/v1/guides/features/command
#[tauri::command]
fn align_sequences(template: &str, reads: &str, alignment_type: &str) -> String {
    let alphabet = Alphabet::default();
    let template = sequence_from_string(template);
    let reads: Vec<Vec<AminoAcid>> = reads.split('\n').map(sequence_from_string).collect();
    let alignment_type = match alignment_type {
        "1" => Type::Local,
        "2" => Type::GlobalForB,
        "3" => Type::Global,
        _ => panic!("Incorrect alignment type"),
    };

    let result = Template::new(
        template,
        reads.iter().map(|a| a.as_slice()).collect(),
        &alphabet,
        alignment_type,
    );
    result.generate_html()
}

#[tauri::command]
fn load_cif(path: &str, min_length: usize, warn: bool) -> Result<(String, String), String> {
    let result = open(path, StrictnessLevel::Loose);
    if let Ok(file) = result {
        let warnings = file.1.into_iter().map(|err| format!("{}", err)).join("\n");
        let pdb = file.0;
        let mut found_unknown = HashMap::new();
        let output = pdb
            .chains()
            .map(|c| {
                c.conformers()
                    .filter_map(|a| {
                        match AMINO_ACIDS
                            .iter()
                            .position(|err| *err == a.name())
                            .and_then(|v| AMINO_ACIDS_CHAR.get(v))
                        {
                            Some(s) => Some(s),
                            None => {
                                if warn && !IGNORE_LIST.contains(&a.name()) {
                                    found_unknown.insert(
                                        a.name(),
                                        1 + found_unknown.get(a.name()).unwrap_or(&0),
                                    );
                                };
                                None
                            }
                        }
                    })
                    .collect::<String>()
            })
            .filter(|a| a.len() >= min_length)
            .join("\n");
        let warnings = warnings + "\n" + &found_unknown.into_iter().map(|name|  {
            format!(
                "{}",
                PDBError::new(
                    ErrorLevel::GeneralWarning,
                    "Unrecognised residue",
                    format!(
                        "This name was not recognised as an Amino Acid or common solvent. It was found {} time{}.",
                        name.1,
                        if name.1 != 1 { "s" } else { "" }
                    ),
                    Context::show(name.0),
                )
            )
        }).join("\n");
        Ok((output, warnings))
    } else {
        Err(result
            .unwrap_err()
            .into_iter()
            .map(|a| format!("{}", a))
            .collect())
    }
}

use rustyms::{
    e, fragment::FragmentType, generate_theoretical_fragments, spectrum::*, AverageWeight, Charge,
    Model,
};
static mut spectra: Vec<RawSpectrum> = Vec::new();

// Learn more about Tauri commands at https://tauri.app/v1/guides/features/command
#[tauri::command]
fn load_mgf(path: &str) -> Result<String, String> {
    if let Ok(v) = rustyms::mgf::open(path) {
        unsafe {
            spectra = v;
        }
        Ok(format!("Loaded {} spectra", unsafe { spectra.len() }))
    } else {
        Err("Could not load mgf file".to_string())
    }
}

#[tauri::command]
fn annotate_spectrum(index: usize, peptide: &str) -> Result<(String, String), String> {
    if index >= unsafe { spectra.len() } {
        return Err("Nonexistent spectrum index".to_string());
    }
    let model = Model::all();
    let peptide = rustyms::peptide::Peptide::pro_forma(peptide)?;
    let spectrum = unsafe { &spectra[index] };
    let fragments =
        generate_theoretical_fragments::<AverageWeight>(&peptide, Charge::new::<e>(1.0), &model);
    let annotated = spectrum.annotate(peptide, fragments, &model);
    Ok((
        render_annotated_spectrum(&annotated, "spectrum"),
        format!("{annotated:?}"),
    ))
}

fn render_annotated_spectrum(spectrum: &AnnotatedSpectrum, id: &str) -> String {
    let mut output = "<div class='spectrum'>".to_string();
    let (limits, overview) = get_overview(spectrum);

    write!(output, "<div class='manual-zoom'>").unwrap();
    write!(output, "<label for='{id}-mz-min'>Mz Min</label>").unwrap();
    write!(
        output,
        "<input id='{id}-mz-min' class='mz-min' type='number' value='0'/>"
    )
    .unwrap();
    write!(output, "<label for='{id}-mz-max'>Mz Max</label>").unwrap();
    write!(
        output,
        "<input id='{id}-mz-max' class='mz-max' type='number' value='{}'/>",
        limits.0.value
    )
    .unwrap();
    write!(
        output,
        "<label for='{id}-intensity-max'>Intensity Max</label>"
    )
    .unwrap();
    write!(
        output,
        "<input id='{id}-intensity-max' class='intensity-max' type='number' value='{}'/>",
        limits.2
    )
    .unwrap();
    write!(output, "</div>").unwrap();
    write!(output, "<div class='wrapper unassigned first'>").unwrap();
    create_ion_legend(&mut output, &format!("{id}-1"));
    render_peptide(&mut output, spectrum, overview);
    render_spectrum(&mut output, spectrum, limits, "first");
    write!(output, "</div>").unwrap();
    // Second spectrum
    write!(output, "</div>").unwrap();
    output
}

fn get_overview(
    spectrum: &AnnotatedSpectrum,
) -> ((MassOverCharge, f64, f64), Vec<HashSet<FragmentType>>) {
    let mut output = vec![HashSet::new(); spectrum.peptide.sequence.len()];
    let mut max_mz: MassOverCharge = MassOverCharge::zero();
    let mut max_intensity: f64 = 0.0;
    let mut max_intensity_unassigned: f64 = 0.0;
    for peak in &spectrum.spectrum {
        max_mz = max_mz.max(peak.experimental_mz);
        max_intensity_unassigned = max_intensity_unassigned.max(peak.intensity);
        if let Some(f) = &peak.annotation {
            max_intensity = max_intensity.max(peak.intensity);
            output[f.sequence_index].insert(f.ion);
        }
    }
    (
        (
            max_mz * 1.01,
            max_intensity.max(max_intensity_unassigned) * 1.01,
            max_intensity_unassigned * 1.01,
        ),
        output,
    )
}

fn create_ion_legend(output: &mut String, id: &str) {
    write!(
        output,
        "<div class='legend'>
    <span class='title'>Ion legend</span>
    <div class='ion-series'>
        <div class='top'>
            <span class='ion w' tabindex='0'>w</span>
            <span class='ion x' tabindex='0'>x</span>
            <span class='ion y' tabindex='0'>y</span>
            <span class='ion z' tabindex='0'>z</span>
        </div><div class='bottom'>
            <span class='ion a' tabindex='0'>a</span>
            <span class='ion b' tabindex='0'>b</span>
            <span class='ion c' tabindex='0'>c</span>
            <span class='ion d' tabindex='0'>d</span>
        </div>
    </div>
    <span class='other'>Other</span>
    <input id='{id}_unassigned' type='checkbox' checked class='unassigned'/>
    <label for='{id}_unassigned' class='unassigned' tabindex='0'>Unassigned</label>
    <label class='label'>
        Ion
        <sup>Charge</sup>
        <sub style='margin-left:-6ch;margin-right:.5rem;'>Position</sub>
        Show for top:
        <input id='{id}_label' type='range' min='0' max='100' value='100'/>
        <input id='{id}_label_value' type='number' min='0' max='100' value='100'/>
        %
    </label>
</div>"
    )
    .unwrap();
}

fn render_peptide(
    output: &mut String,
    spectrum: &AnnotatedSpectrum,
    overview: Vec<HashSet<FragmentType>>,
) {
    write!(output, "<div class='peptide'>").unwrap();
    if spectrum.peptide.n_term.is_some() {
        write!(output, "<span class='modification term'></span>").unwrap();
    }
    for (index, (pos, ions)) in spectrum.peptide.sequence.iter().zip(overview).enumerate() {
        let mut classes = String::new();
        if pos.1.is_some() {
            write!(classes, " class='modification'").unwrap();
        }
        write!(
            output,
            "<span data-pos='{}'{classes} tabindex='0'>{}",
            index + 1,
            pos.0.char()
        )
        .unwrap();
        for ion in ions {
            write!(output, "<span class='corner {ion}'></span>").unwrap();
        }
        write!(output, "</span>").unwrap();
    }
    if spectrum.peptide.c_term.is_some() {
        write!(output, "<span class='modification term'></span>").unwrap();
    }
    write!(output, "</div>").unwrap();
}

fn render_spectrum(
    output: &mut String,
    spectrum: &AnnotatedSpectrum,
    limits: (MassOverCharge, f64, f64),
    selection: &str,
) {
    write!(
        output,
        "<div class='canvas-wrapper label' aria-hidden='true'>"
    )
    .unwrap();
    write!(output, "<div class='y-axis'><span class='n0'>0</span>").unwrap();
    write!(output, "<span>{}</span>", limits.2 / 4.0).unwrap();
    write!(output, "<span>{}</span>", limits.2 / 2.0).unwrap();
    write!(output, "<span>{}</span>", 3.0 * limits.2 / 4.0).unwrap();
    write!(output, "<span class='last'>{}</span>", limits.2).unwrap();
    write!(output, "</div>").unwrap();
    write!(output, "<div class='canvas' style='--min-mz:0;--max-mz:{};--max-intensity:{};' data-initial-max-mz='{}' data-initial-max-intensity='{}' data-initial-max-intensity-assigned='{}'>", limits.0.value, limits.2, limits.0.value, limits.2, limits.1).unwrap();
    write!(
        output,
        "<span class='selection {selection}' hidden='true'></span>"
    )
    .unwrap();
    write!(output, "<div class='zoom-out' tabindex='0'>Zoom Out</div>").unwrap();

    for peak in &spectrum.spectrum {
        match &peak.annotation {
            None => {
                write!(
                    output,
                    "<span class='peak unassigned' style='--mz:{};--intensity:{};'></span>",
                    peak.experimental_mz.value, peak.intensity
                )
                .unwrap();
            }
            Some(f) => {
                write!(
                    output,
                    "<span class='peak {} label' style='--mz:{};--intensity:{};' data-pos='{}'>",
                    f.ion,
                    peak.experimental_mz.value,
                    peak.intensity,
                    f.sequence_index + 1
                )
                .unwrap();
                let ch = format!("{:+}", peak.charge.value);
                write!(
                    output,
                    "<span>{}<sup>{}</sup><sub style='margin-left:-{}ch'>{}</sub></span></span>",
                    f.ion,
                    ch,
                    ch.len(),
                    f.sequence_index
                )
                .unwrap();
            }
        }
    }
    write!(output, "</div>").unwrap();
    write!(output, "<div class='x-axis'><span class='n0'>0</span>").unwrap();
    write!(output, "<span>{}</span>", limits.0.value / 4.0).unwrap();
    write!(output, "<span>{}</span>", limits.0.value / 2.0).unwrap();
    write!(output, "<span>{}</span>", 3.0 * limits.0.value / 4.0).unwrap();
    write!(output, "<span class='last'>{}</span>", limits.0.value).unwrap();
    write!(output, "</div></div>").unwrap();
}

/// All amino acids. Includes Amber-specific naming conventions for (de-)protonated versions, CYS involved in
/// disulfide bonding and the like.
const AMINO_ACIDS: &[&str] = &[
    "ALA", "ARG", "ASH", "ASN", "ASP", "ASX", "CYS", "CYX", "GLH", "GLN", "GLU", "GLY", "HID",
    "HIE", "HIM", "HIP", "HIS", "ILE", "LEU", "LYN", "LYS", "MET", "PHE", "PRO", "SER", "THR",
    "TRP", "TYR", "VAL", "SEC", "PYL",
];

const AMINO_ACIDS_CHAR: &[char] = &[
    'A', 'R', 'N', 'N', 'D', 'B', 'C', 'C', 'Q', 'Q', 'E', 'G', 'H', 'H', 'H', 'H', 'H', 'I', 'L',
    'K', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'U', 'O',
];

const IGNORE_LIST: &[&str] = &["HOH", "WAT", "ADP", "DMS"]; // Common solvents I recognised

fn main() {
    tauri::Builder::default()
        .invoke_handler(tauri::generate_handler![
            align_sequences,
            load_cif,
            load_mgf,
            annotate_spectrum
        ])
        .run(tauri::generate_context!())
        .expect("error while running tauri application");
}
