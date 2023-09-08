use std::{collections::HashSet, fmt::Write};

use itertools::Itertools;
use rustyms::{
    AnnotatedSpectrum, Charge, ComplexPeptide, Fragment, FragmentType, LinearPeptide, Mass,
    MassOverCharge, MolecularFormula, Zero,
};

// TODO: finish multimeric:
// - [ ] Fix highlighting (now class based)
// - [ ] Add highlighting on a peptide level (see all ions for one peptide)
// - [ ] Fix the output if only one peptide is supplied
// - [ ] Fix general stat table
// - [ ] Check what is wrong with the annotation for chimeric model CidHcd
// - [ ] Do some more thinking on the multi annotation display for single peaks
// - [ ] Show peptides horizontally besides each other if there is space (justified to both page ends?)
// - [ ] Show peptide number somewhere? (the number should be clear any ways, but maybe this is nice, or also nice to clearly separate both)
// Some more general things:
// - [ ] Make sure the change of sequence positions covered is correct
// - [ ] Fix colour of precursor peaks
// - [ ] Make all ion series automatically ignore the very last pos (with all AA, so not feasible to be generated)

use crate::html_builder::{HtmlElement, HtmlTag};

pub fn fragment_table(fragments: &[Fragment]) -> String {
    let mut output = "<table><thead><tr><th>Sequence Index</th><th>SeriesNumber</th><th>Ion</th><th>mz</th><th>Charge</th><th>Neutral loss</th></tr></thead><tbody>".to_string();
    for fragment in fragments {
        write!(
            &mut output,
            "<tr><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td></tr>",
            fragment
                .ion
                .position()
                .map(|p| p.sequence_index.to_string())
                .unwrap_or("-".to_string()),
            fragment
                .ion
                .position()
                .map(|p| p.series_number.to_string())
                .unwrap_or("-".to_string()),
            fragment.ion,
            fragment.mz().map_or(f64::NAN, |v| v.value),
            fragment.charge.value,
            fragment
                .neutral_loss
                .map(|n| n.to_string())
                .unwrap_or("-".to_string()),
        )
        .unwrap();
    }
    write!(&mut output, "</tbody></table>").unwrap();
    output
}

pub fn annotated_spectrum(
    spectrum: &AnnotatedSpectrum,
    id: &str,
    fragments: &[Fragment],
) -> String {
    let mut output = String::new();
    let (limits, overview) = get_overview(spectrum);
    let (graph_data, graph_boundaries) = spectrum_graph_boundaries(spectrum, fragments);

    spectrum_top_buttons(&mut output, id, &limits, &graph_boundaries).unwrap();

    create_ion_legend(&mut output, &format!("{id}-1"));
    render_peptide(&mut output, spectrum, overview);
    render_spectrum(&mut output, spectrum, &graph_boundaries, limits, "first");
    // Spectrum graph
    spectrum_graph(&mut output, &graph_boundaries, &graph_data, limits.0.value);
    write!(output, "</div></div>").unwrap();
    // General stats
    general_stats(&mut output, spectrum, fragments);
    // Spectrum table
    collapsible(
        &mut output,
        &format!("{id}-table-1"),
        "Peaks table".to_string(),
        spectrum_table(spectrum, &format!("{id}-table-1")),
    );

    //write!(output, "</div>").unwrap();
    output
}

fn spectrum_top_buttons(
    output: &mut String,
    id: &str,
    limits: &(MassOverCharge, f64, f64),
    spectrum_graph_boundaries: &(f64, f64, f64, f64, f64, f64, f64, f64, f64, f64),
) -> core::fmt::Result {
    write!(output, "<div class='settings manual-zoom'><p>Spectrum</p>")?;
    write!(output, "<label for='{id}-mz-min'>Mz Min</label>")?;
    write!(
        output,
        "<input id='{id}-mz-min' class='mz-min' type='number' value='0'/>"
    )?;
    write!(output, "<label for='{id}-mz-max'>Mz Max</label>")?;
    write!(
        output,
        "<input id='{id}-mz-max' class='mz-max' type='number' value='{}'/>",
        limits.0.value
    )?;
    write!(
        output,
        "<label for='{id}-intensity-max'>Intensity Max</label>"
    )?;
    write!(
        output,
        "<input id='{id}-intensity-max' class='intensity-max' type='number' value='{}'/>",
        limits.2
    )?;
    write!(output, "</div>")?;
    Ok(())
}

type Boundaries = (f64, f64, f64, f64, f64, f64, f64, f64, f64, f64);
type SpectrumGraphData = Vec<(
    Vec<Fragment>,
    (f64, Fragment),
    (f64, Fragment),
    MassOverCharge,
    Mass,
    f64,
)>;

fn spectrum_graph_boundaries(
    spectrum: &AnnotatedSpectrum,
    fragments: &[Fragment],
) -> (SpectrumGraphData, Boundaries) {
    let data: SpectrumGraphData = spectrum
        .spectrum
        .iter()
        .map(|point| {
            let distance = fragments
                .iter()
                .filter(|frag| frag.ion.to_string() == "c" || frag.ion.to_string() == "z")
                .fold(
                    (
                        (
                            f64::MAX,
                            Fragment::new(
                                MolecularFormula::default(),
                                Charge::zero(),
                                0,
                                FragmentType::precursor,
                                String::new(),
                            ),
                        ),
                        (
                            f64::MAX,
                            Fragment::new(
                                MolecularFormula::default(),
                                Charge::zero(),
                                0,
                                FragmentType::precursor,
                                String::new(),
                            ),
                        ),
                    ),
                    |acc, frag: &Fragment| {
                        if let Some(mz) = frag.mz() {
                            let rel = ((mz - point.experimental_mz) / mz
                                * MassOverCharge::new::<rustyms::mz>(1e6))
                            .value;
                            let abs = (mz - point.experimental_mz).value;
                            (
                                if acc.0 .0.abs() < rel.abs() {
                                    acc.0
                                } else {
                                    (rel, frag.clone())
                                },
                                if acc.1 .0.abs() < abs.abs() {
                                    acc.1
                                } else {
                                    (abs, frag.clone())
                                },
                            )
                        } else {
                            acc
                        }
                    },
                );
            (
                point.annotation.clone(),
                distance.0,                           // rel (ppm)
                distance.1,                           // abs (Da)
                point.experimental_mz,                // mz
                point.experimental_mz * point.charge, // mass
                point.intensity,                      // intensity
            )
        })
        .collect();
    let bounds = data.iter().fold(
        (
            f64::MIN,
            f64::MAX,
            f64::MIN,
            f64::MAX,
            f64::MIN,
            f64::MAX,
            f64::MIN,
            f64::MAX,
            f64::MIN,
            f64::MAX,
        ),
        |acc, point| {
            (
                acc.0.max(point.1 .0), // rel
                acc.1.min(point.1 .0),
                acc.2.max(point.2 .0), // abs
                acc.3.min(point.2 .0),
                acc.4.max(point.3.value), // mz
                acc.5.min(point.3.value),
                acc.6.max(point.4.value), // mass
                acc.7.min(point.4.value),
                acc.8.max(point.5), // intensity
                acc.9.min(point.5),
            )
        },
    );
    (data, bounds)
}

fn spectrum_graph(
    output: &mut String,
    boundaries: &Boundaries,
    data: &SpectrumGraphData,
    x_max: f64,
) {
    write!(output, "<div class='spectrum-graph-y-axis'>").unwrap();
    write!(output, "<span class='max'>{:.2}</span>", boundaries.2).unwrap();
    write!(
        output,
        "<span class='title abs'>Absolute distance to closest c/z ion (Da)</span><span class='title rel'>Relative distance to closest c/z ion (ppm)</span>"
    )
    .unwrap();
    write!(output, "<span class='min'>{:.2}</span>", boundaries.3).unwrap();
    density_graph::<256>(
        output,
        &data.iter().map(|p| p.1 .0).collect::<Vec<_>>(),
        "rel",
    );
    density_graph::<256>(
        output,
        &data.iter().map(|p| p.2 .0).collect::<Vec<_>>(),
        "abs",
    );
    write!(output, "</div>").unwrap();
    write!(output, "<div class='spectrum-graph canvas'>",).unwrap();
    write!(output, "<div class='x-axis'>").unwrap();
    write!(output, "<span class='min'>0</span>").unwrap();
    write!(output, "<span class='max'>{:.2}</span>", x_max).unwrap();
    write!(output, "</div>").unwrap();
    write!(
        output,
        "<span class='ruler'><span id='ruler-value'>XX</span></span>"
    )
    .unwrap();

    for point in data {
        write!(
            output,
            "<span class='point {}' style='--rel:{};--abs:{};--mz:{};--intensity:{}' data-ppm='{}' data-abs='{}' data-mz='{}' data-intensity='{}' data-reference-fragment-rel='{}' data-reference-fragment-abs='{}'></span>",
            get_classes(&point.0),
            point.1.0,
            point.2.0,
            point.3.value,
            point.5,
            point.1.0,
            point.2.0,
            point.3.value,
            point.5,
            point.1.1,
            point.2.1,
        )
        .unwrap();
    }
    write!(output, "</div>").unwrap();
}

/// Get all applicable classes for a set of annotations.
/// These are:
///   * The ion type(s) (eg 'precursor', 'y')
///   * The peptide(s) (eg 'p0', 'p1')
///   * The position(s) (eg 'p0-2', 'p2-123')
fn get_classes(annotations: &[Fragment]) -> String {
    let mut output = Vec::new();
    let mut shared_ion = annotations.get(0).map(|a| a.ion.to_string());
    for annotation in annotations {
        output.push(annotation.ion.to_string());
        output.push(format!("p{}", annotation.peptide_index));
        if let Some(pos) = annotation.ion.position() {
            output.push(format!(
                "p{}-{}",
                annotation.peptide_index, pos.sequence_index
            ));
        }
        if let Some(ion) = &shared_ion {
            if *ion != annotation.ion.to_string() {
                shared_ion = None;
            }
        }
    }
    if shared_ion.is_none() {
        output.push("multi".to_string())
    }
    output = output.into_iter().unique().collect();
    if annotations.is_empty() {
        "unassigned".to_string()
    } else {
        format!("label {}", output.join(" "))
    }
}

fn get_label(annotations: &[Fragment]) -> String {
    if annotations.is_empty() {
        String::new()
    } else {
        let mut shared_charge = Some(annotations[0].charge);
        let mut shared_ion = Some(annotations[0].ion.to_string());
        let mut shared_pos = Some(annotations[0].ion.position());
        let mut shared_peptide = Some(annotations[0].peptide_index);
        for a in annotations {
            if let Some(charge) = shared_charge {
                if charge != a.charge {
                    shared_charge = None;
                }
            }
            if let Some(ion) = &shared_ion {
                if *ion != a.ion.to_string() {
                    shared_ion = None;
                }
            }
            if let Some(pos) = shared_pos {
                if pos != a.ion.position() {
                    shared_pos = None;
                }
            }
            if let Some(peptide) = shared_peptide {
                if peptide != a.peptide_index {
                    shared_peptide = None;
                }
            }
        }

        if shared_charge.is_none()
            && shared_ion.is_none()
            && shared_pos.is_none()
            && shared_peptide.is_none()
        {
            "*".to_string()
        } else {
            let charge_str = shared_charge
                .map(|charge| format!("{:+}", charge.value))
                .unwrap_or("*".to_string());
            let ion_str = shared_ion.unwrap_or("*".to_string());
            let pos_str = shared_pos
                .map(|pos| {
                    pos.map(|p| p.series_number.to_string())
                        .unwrap_or(String::new())
                })
                .unwrap_or("*".to_string());
            let peptide_str = shared_peptide
                .map(|pep| (pep + 1).to_string())
                .unwrap_or("*".to_string());

            let multi = if annotations.len() > 1 {
                let mut multi = String::new();
                for annotation in annotations {
                    let ch = format!("{:+}", annotation.charge.value);
                    write!(
                        multi,
                        "<span>{}<sup>{}</sup><sub style='margin-left:-{}ch'>{}p{}</sub></span>",
                        annotation.ion,
                        ch,
                        ch.len(),
                        annotation
                            .ion
                            .position()
                            .map(|p| p.series_number.to_string())
                            .unwrap_or(String::new()),
                        annotation.peptide_index + 1,
                    )
                    .unwrap();
                }
                format!("<span class='multi'>{multi}</span>")
            } else {
                String::new()
            };

            format!(
                "<span>{}<sup>{}</sup><sub style='margin-left:-{}ch'>{}p{}</sub></span>{}",
                ion_str,
                charge_str,
                charge_str.len(),
                pos_str,
                peptide_str,
                multi
            )
        }
    }
}

fn spectrum_graph_header(
    output: &mut String,
    boundaries: &(f64, f64, f64, f64, f64, f64, f64, f64, f64, f64),
) -> std::fmt::Result {
    write!(output, "<input type='radio' name='y-axis' id='absolute' value='absolute' checked/><label for='absolute'>Absolute</label>")?;
    write!(output, "<input type='radio' name='y-axis' id='relative' value='relative'/><label for='relative'>Relative</label>")?;
    write!(output, "<input type='checkbox' name='intensity' id='intensity' value='intensity'/><label for='intensity'>Intensity</label>")?;
    write!(output, "<div class='manual-zoom'>")?;
    write!(output, "<label for='y-min'>Y Min</label>")?;
    write!(
        output,
        "<input id='y-min' class='y-min' type='number' value='{}'/>",
        boundaries.3
    )?;
    write!(output, "<label for='y-max'>Y Max</label>")?;
    write!(
        output,
        "<input id='y-max' class='y-max' type='number' value='{}'/>",
        boundaries.2
    )?;
    write!(output, "</div>")?;
    Ok(())
}

type PositionCoverage = Vec<Vec<HashSet<FragmentType>>>;

fn get_overview(spectrum: &AnnotatedSpectrum) -> ((MassOverCharge, f64, f64), PositionCoverage) {
    let mut output: PositionCoverage = spectrum
        .peptide
        .peptides()
        .iter()
        .map(|p| vec![HashSet::new(); p.sequence.len()])
        .collect();
    let mut max_mz: MassOverCharge = MassOverCharge::zero();
    let mut max_intensity: f64 = 0.0;
    let mut max_intensity_unassigned: f64 = 0.0;
    for peak in &spectrum.spectrum {
        max_mz = max_mz.max(peak.experimental_mz);
        max_intensity_unassigned = max_intensity_unassigned.max(peak.intensity);
        if !peak.annotation.is_empty() {
            max_intensity = max_intensity.max(peak.intensity);
            peak.annotation.iter().for_each(|frag| {
                frag.ion
                    .position()
                    .map(|i| output[frag.peptide_index][i.sequence_index].insert(frag.ion));
            });
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
            <span class='ion c-term' tabindex='0'>C-term</span>
        </div><div class='bottom'>
            <span class='ion n-term' tabindex='0'>N-term</span>
            <span class='ion a' tabindex='0'>a</span>
            <span class='ion b' tabindex='0'>b</span>
            <span class='ion c' tabindex='0'>c</span>
            <span class='ion d' tabindex='0'>d</span>
            <span class='ion v' tabindex='0'>v</span>
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
    <input id='{id}_mass_label' type='checkbox' class='mass-label'/>
    <label for='{id}_mass_label' class='mass-label' tabindex='0'>Show top masses</label>
</div>"
    )
    .unwrap();
}

fn render_peptide(output: &mut String, spectrum: &AnnotatedSpectrum, overview: PositionCoverage) {
    write!(output, "<div class='complex-peptide'>").unwrap();
    match &spectrum.peptide {
        ComplexPeptide::Singular(peptide) => {
            render_linear_peptide(output, peptide, &overview[0], 0)
        }
        ComplexPeptide::Multimeric(peptides) => {
            for (index, peptide) in peptides.iter().enumerate() {
                render_linear_peptide(output, peptide, &overview[index], index);
            }
        }
    }
    write!(output, "</div>").unwrap();
}

fn render_linear_peptide(
    output: &mut String,
    peptide: &LinearPeptide,
    overview: &[HashSet<FragmentType>],
    peptide_index: usize,
) {
    write!(output, "<div class='peptide'>").unwrap();
    if peptide.n_term.is_some() {
        write!(output, "<span class='modification term'></span>").unwrap();
    }
    for (index, (pos, ions)) in peptide.sequence.iter().zip(overview).enumerate() {
        let mut classes = String::new();
        if !pos.modifications.is_empty() {
            write!(classes, " class='modification'").unwrap();
        }
        write!(
            output,
            "<span data-pos='{}-{}'{classes} tabindex='0' title='N terminal position: {}, C terminal position: {}'>{}",
            peptide_index+1,
            index + 1,
            index + 1,
            peptide.sequence.len() - index,
            pos.aminoacid.char()
        )
        .unwrap();
        for ion in ions {
            write!(
                output,
                "<span class='corner {}'></span>",
                ion.to_string().trim_end_matches('·')
            )
            .unwrap();
        }
        write!(output, "</span>").unwrap();
    }
    if peptide.c_term.is_some() {
        write!(output, "<span class='modification term'></span>").unwrap();
    }
    write!(output, "</div>").unwrap();
}

fn render_spectrum(
    output: &mut String,
    spectrum: &AnnotatedSpectrum,
    boundaries: &Boundaries,
    limits: (MassOverCharge, f64, f64),
    selection: &str,
) {
    write!(
        output,
        "<div class='canvas-wrapper label' aria-hidden='true' style='--min-mz:0;--max-mz:{0};--max-intensity:{1};--y-max:{2};--y-min:{3};--abs-max-initial:{2};--abs-min-initial:{3};--rel-max-initial:{4};--rel-min-initial:{5};' data-initial-max-mz='{0}' data-initial-max-intensity='{1}' data-initial-max-intensity-assigned='{6}'>",
        limits.0.value, limits.2, boundaries.2, boundaries.3, boundaries.0, boundaries.1, limits.1
    )
    .unwrap();
    write!(
        output,
        "<div class='y-axis'><span class='tick'><span class='n0'>0</span></span>"
    )
    .unwrap();
    write!(
        output,
        "<span class='tick'><span>{:.2}</span></span>",
        limits.2 / 4.0
    )
    .unwrap();
    write!(
        output,
        "<span class='tick'><span>{:.2}</span></span>",
        limits.2 / 2.0
    )
    .unwrap();
    write!(
        output,
        "<span class='tick'><span>{:.2}</span></span>",
        3.0 * limits.2 / 4.0
    )
    .unwrap();
    write!(
        output,
        "<span class='tick'><span class='last'>{:.2}</span></span>",
        limits.2
    )
    .unwrap();
    write!(output, "</div>").unwrap();
    write!(output, "<div class='canvas'>").unwrap();
    write!(
        output,
        "<span class='selection {selection}' hidden='true'></span>"
    )
    .unwrap();
    write!(
        output,
        "<div class='zoom-out' tabindex='0'>Reset Zoom</div>"
    )
    .unwrap();

    for peak in &spectrum.spectrum {
        write!(
            output,
            "<span class='peak {}' style='--mz:{};--intensity:{};' data-label='{}'>{}</span>",
            get_classes(&peak.annotation),
            peak.experimental_mz.value,
            peak.intensity,
            (peak.experimental_mz.value * 10.0).round() / 10.0,
            get_label(&peak.annotation),
        )
        .unwrap();
    }
    write!(output, "</div>").unwrap();
    write!(
        output,
        "<div class='x-axis'><span class='tick'><span class='n0'>0</span></span>"
    )
    .unwrap();
    write!(
        output,
        "<span class='tick'><span>{:.2}</span></span>",
        limits.0.value / 4.0
    )
    .unwrap();
    write!(
        output,
        "<span class='tick'><span>{:.2}</span></span>",
        limits.0.value / 2.0
    )
    .unwrap();
    write!(
        output,
        "<span class='tick'><span>{:.2}</span></span>",
        3.0 * limits.0.value / 4.0
    )
    .unwrap();
    write!(
        output,
        "<span class='tick'><span class='last'>{:.2}</span></span>",
        limits.0.value
    )
    .unwrap();
    write!(output, "</div>").unwrap();
}

fn spectrum_table(spectrum: &AnnotatedSpectrum, table_id: &str) -> String {
    let mut output = String::new();
    write!(
        output,
        "<label class='background'><input type='checkbox'/>Show background peaks</label>
        <table id='{table_id}' class='wide-table'>
            <tr>
                <th>Position</th>
                <th>Ion type</th>
                <th>Loss</th>
                <th>Intensity</th>
                <th>mz Theoretical</th>
                <th>mz Error (Th)</th>
                <th>mz Error (ppm)</th>
                <th>Charge</th>
                <th>Series Number</th>
            </tr>"
    )
    .unwrap();
    for peak in &spectrum.spectrum {
        if peak.annotation.is_empty() {
            write!(
                output,
                "<tr class='unassigned'>
                <td>-</td>
                <td>-</td>
                <td>-</td>
                <td>{:.2}</td>
                <td>{:.2}</td>
                <td>-</td>
                <td>-</td>
                <td>{:+}</td>
                <td>-</td>
            </tr>",
                peak.intensity, peak.experimental_mz.value, peak.charge.value
            )
            .unwrap();
        } else {
            for annotation in &peak.annotation {
                write!(
                    output,
                    "<tr>
                    <td>{}</td>
                    <td>{}</td>
                    <td>{}</td>
                    <td>{:.2}</td>
                    <td>{:.2}</td>
                    <td>{:.5}</td>
                    <td>{:.2}</td>
                    <td>{:+}</td>
                    <td>{}</td>
                </tr>",
                    annotation
                        .ion
                        .position()
                        .map_or("*".to_string(), |i| (i.sequence_index + 1).to_string()),
                    annotation.ion,
                    annotation
                        .neutral_loss
                        .map_or(String::new(), |v| v.to_string()),
                    peak.intensity,
                    peak.experimental_mz.value,
                    annotation
                        .mz()
                        .map_or(f64::NAN, |v| (v - peak.experimental_mz).abs().value),
                    annotation
                        .mz()
                        .map_or(f64::NAN, |v| ((v - peak.experimental_mz).abs() / v * 1e6)
                            .value),
                    annotation.charge.value,
                    annotation
                        .ion
                        .position()
                        .map_or("*".to_string(), |i| i.series_number.to_string()),
                )
                .unwrap();
            }
        }
    }
    write!(output, "</table>").unwrap();
    output
}

fn general_stats(output: &mut String, spectrum: &AnnotatedSpectrum, fragments: &[Fragment]) {
    let mut mass_row = String::new();
    let mut fragments_row = String::new();
    let mut peaks_row = String::new();
    let mut intensity_row = String::new();
    let mut positions_row = String::new();

    for (peptide_index, peptide) in spectrum.peptide.peptides().iter().enumerate() {
        let precursor = if let Some(formula) = peptide.formula() {
            if let (Some(mono), Some(avg)) = (formula.monoisotopic_mass(), formula.average_weight())
            {
                format!("{:.3} Da | avg: {:.3} Da", mono.value, avg.value)
            } else {
                "No defined mass for precursor".to_string()
            }
        } else {
            "No defined molecular formula for precursor".to_string()
        };
        let (num_annotated, intensity_annotated) = spectrum
            .spectrum
            .iter()
            .filter(|p| {
                p.annotation
                    .iter()
                    .any(|a| a.peptide_index == peptide_index)
            })
            .fold((0, 0.0), |(n, intensity), p| {
                (n + 1, intensity + p.intensity)
            });
        let total_intensity: f64 = spectrum.spectrum.iter().map(|p| p.intensity).sum();
        let percentage_fragments_found = num_annotated as f64 / fragments.len() as f64 * 100.0;
        let percentage_peaks_annotated =
            num_annotated as f64 / spectrum.spectrum.len() as f64 * 100.0;
        let percentage_intensity_annotated = intensity_annotated / total_intensity * 100.0;
        let percentage_positions_covered = spectrum
            .spectrum
            .iter()
            .flat_map(|p| {
                p.annotation
                    .iter()
                    .filter(|a| a.peptide_index == peptide_index)
                    .filter_map(|a| a.ion.position())
            })
            .map(|pos| pos.sequence_index)
            .unique()
            .count() as f64
            / peptide.len() as f64
            * 100.0;
        write!(mass_row, "<td>{precursor}</td>").unwrap();
        write!(
            fragments_row,
            "<td>{percentage_fragments_found:.2} % ({num_annotated}/{})</td>",
            fragments.len()
        )
        .unwrap();
        write!(
            peaks_row,
            "<td>{percentage_peaks_annotated:.2} % ({num_annotated}/{})</td>",
            spectrum.spectrum.len()
        )
        .unwrap();
        write!(
            intensity_row,
            "<td>{percentage_intensity_annotated:.2} %</td>"
        )
        .unwrap();
        write!(
            positions_row,
            "<td>{percentage_positions_covered:.2} %</td>"
        )
        .unwrap();
    }

    if spectrum.peptide.peptides().len() > 1 {
        write!(output, "<tr><td>Stat</td>").unwrap();
        for i in 0..spectrum.peptide.peptides().len() {
            write!(output, "<td>{}</td>", i + 1).unwrap();
        }
        write!(output, "</tr>").unwrap();
    }
    write!(
        output,
        "<table class='general-stats'>
    <tr><td>Precursor Mass (M)</td>{mass_row}</tr>
    <tr><td>Fragments found</td>{fragments_row}</tr>
    <tr><td>Peaks annotated</td>{peaks_row}</tr>
    <tr><td>Intensity annotated</td>{intensity_row}</tr>
    <tr><td>Sequence positions covered</td>{positions_row}</tr>
    </table>"
    )
    .unwrap();
}

fn collapsible(output: &mut String, id: &str, title: String, content: String) {
    write!(
        output,
        "{}{}{}",
        HtmlElement::new(HtmlTag::input)
            .header("type", "checkbox")
            .id(format!("collapsible-{id}")),
        HtmlElement::new(HtmlTag::label)
            .header("for", format!("collapsible-{id}"))
            .content(title),
        HtmlElement::new(HtmlTag::div)
            .class("collapsible")
            .id(id)
            .content(HtmlElement::new(HtmlTag::div).class("clear-fix"))
            .content(content),
    )
    .unwrap();
}

fn density_estimation<const STEPS: usize>(data: &[f64]) -> [f64; STEPS] {
    let mut densities = [0.0; STEPS];
    if data.is_empty() {
        return densities;
    }
    let mut data = data.to_vec();
    data.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());
    let len = data.len() as f64;
    let min_value = data.iter().copied().reduce(f64::min).unwrap_or(f64::MIN);
    let max_value = data.iter().copied().reduce(f64::max).unwrap_or(f64::MAX);
    let mean: f64 = data.iter().sum::<f64>() / len;
    let stdev: f64 = (data.iter().map(|p| (mean - p).powi(2)).sum::<f64>() / len).sqrt();
    let half = len / 2.0;
    let iqr: f64 = data.iter().rev().take(half.floor() as usize).sum::<f64>() / half
        - data.iter().take(half.floor() as usize).sum::<f64>() / half;
    let h = 0.9 * stdev.min(iqr / 1.34) * len.powf(-0.2);

    // let kde = |x: f64| {
    //     1.0 / ((2.0 * std::f64::consts::PI).sqrt() * len * h)
    //         * data
    //             .iter()
    //             .map(|i| (-1.0 / (2.0 * h * h) * (x - i).powi(2)).exp())
    //             .sum::<f64>()
    // };

    let gaussian_kernel =
        |x: f64| 1.0 / (2.0 * std::f64::consts::PI).sqrt() * (-1.0 / 2.0 * x.powi(2)).exp();
    let kde = |x: f64| {
        1.0 / (len * h)
            * data
                .iter()
                .map(|i| gaussian_kernel((x - i) / h))
                .sum::<f64>()
    };

    for (i, density) in densities.iter_mut().enumerate() {
        *density = kde(min_value + (max_value - min_value) / (STEPS - 1) as f64 * i as f64);
    }

    densities
}

fn density_graph<const STEPS: usize>(output: &mut String, data: &[f64], class: &str) {
    let densities = density_estimation::<STEPS>(data);
    let max_density = densities
        .iter()
        .copied()
        .reduce(f64::max)
        .unwrap_or(f64::MAX);
    let mut path = String::new();
    for (i, density) in densities.iter().enumerate() {
        write!(
            &mut path,
            "{}{} {}",
            if i != 0 { " L " } else { "" },
            (max_density - density) / max_density * 100.0,
            i,
        )
        .unwrap();
    }
    write!(output, "<svg viewBox='-1 0 100 {}' preserveAspectRatio='none'><g class='density {class}'><path class='line' d='M {}'></path><path class='volume' d='M 100 0 L {} L {} 0 Z'></path></g></svg>",
STEPS-1,
path,
path, STEPS -1
).unwrap();
}

#[cfg(test)]
mod tests {
    use super::density_estimation;

    #[test]
    fn density_flat() {
        let d = dbg!(density_estimation::<10>(&[1.0, 2.0, 3.0, 4.0, 5.0]));
        assert!((d.iter().sum::<f64>() - 1.0).abs() < 0.005); // half a percent deviation from 1.0 density for the whole distribution
        assert_eq!(
            d.iter().copied().reduce(f64::min),
            d.iter().copied().reduce(f64::max)
        );
    }

    #[test]
    fn density_peak() {
        let d = dbg!(density_estimation::<10>(&[
            1.0, 2.0, 3.0, 4.0, 5.0, 2.0, 3.0, 4.0, 3.0
        ]));
        assert!((d.iter().sum::<f64>() - 1.0).abs() < 0.005); // half a percent deviation from 1.0 density for the whole distribution
        todo!()
    }

    #[test]
    fn density_sample() {
        let d = dbg!(density_estimation::<10>(&[-2.1, -1.3, -0.4, 1.9, 5.1, 6.2]));
        assert!((d.iter().sum::<f64>() - 1.0).abs() < 0.005); // half a percent deviation from 1.0 density for the whole distribution
        todo!()
    }
}
