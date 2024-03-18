use std::{cmp::Ordering, collections::HashSet, fmt::Write};

use itertools::{all, Itertools};
use rustyms::{
    fragment::*,
    model::Location,
    spectrum::{AnnotatedPeak, PeakSpectrum},
    system::*,
    AnnotatedSpectrum, ComplexPeptide, Element, LinearPeptide, MassMode, Model, Modification,
    MolecularFormula, MultiChemical,
};

use crate::html_builder::{HtmlElement, HtmlTag};

pub fn fragment_table(fragments: &[Fragment], multiple_peptides: bool) -> String {
    let mut output = format!("<table><thead><tr><th>Sequence Index</th><th>Series Number</th><th>Ion</th><th>mz</th><th>Charge</th><th>Neutral loss</th>{}</tr></thead><tbody>", if multiple_peptides {
        "<td>Peptide</td>"
    } else {
        ""
    });
    for fragment in fragments {
        let (sequence_index, series_number) = if let Some(pos) = fragment.ion.position() {
            (
                pos.sequence_index.to_string(),
                pos.series_number.to_string(),
            )
        } else if let Some(pos) = fragment.ion.glycan_position() {
            (pos.attachment(), pos.series_number.to_string())
        } else if let FragmentType::InternalGlycan(breakages) = &fragment.ion {
            (
                breakages
                    .get(0)
                    .map(|b| b.position().attachment())
                    .unwrap_or("-".to_string()),
                breakages
                    .iter()
                    .map(std::string::ToString::to_string)
                    .join(""),
            )
        } else {
            // precursor
            ("-".to_string(), "-".to_string())
        };
        write!(
            &mut output,
            "<tr><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td>{}</tr>",
            sequence_index,
            series_number,
            fragment.ion.label(),
            fragment.mz(MassMode::Monoisotopic).value,
            fragment.charge.value,
            fragment
                .neutral_loss
                .as_ref()
                .map_or("-".to_string(), |n| n.to_string()),
            if multiple_peptides {
                format!("<td>{}</td>", fragment.peptide_index + 1)
            } else {
                String::new()
            }
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
    model: &Model,
) -> (String, Limits) {
    let mut output = String::new();
    let (limits, overview) = get_overview(spectrum);
    let (graph_data, graph_boundaries, ions) =
        spectrum_graph_boundaries(spectrum, fragments, model);
    let multiple_peptides = !matches!(spectrum.peptide, ComplexPeptide::Singular(_));
    let multiple_glycans = spectrum.peptide.peptides().iter().any(|p| {
        p.sequence
            .iter()
            .filter(|seq| {
                seq.modifications
                    .iter()
                    .any(|m| matches!(m, Modification::GlycanStructure(_)))
            })
            .count()
            > 1
    });

    render_peptide(&mut output, spectrum, overview, multiple_peptides);
    render_spectrum(
        &mut output,
        spectrum,
        fragments,
        &graph_boundaries,
        &limits,
        "first",
        multiple_peptides,
        multiple_glycans,
    );
    // Spectrum graph
    spectrum_graph(
        &mut output,
        &graph_boundaries,
        &graph_data,
        &ions,
        limits.mz.value,
    );
    write!(output, "</div></div>").unwrap();
    // General stats
    general_stats(&mut output, spectrum, fragments);

    //write!(output, "</div>").unwrap();
    (output, limits)
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
    model: &Model,
) -> (SpectrumGraphData, Boundaries, String) {
    // let n = if model.c.0 != Location::None {
    //     FragmentType::c(_)
    // } else if model.b.0 != Location::None {
    //     FragmentType::b(_)
    // } else {
    //     FragmentType::a(_)
    // };
    // let c = if model.z.0 != Location::None {
    //     FragmentType::z(_)
    // } else if model.y.0 != Location::None {
    //     FragmentType::y(_)
    // } else {
    //     FragmentType::x(_)
    // };
    fn filter(model: &Model, ion: &FragmentType) -> bool {
        let n = if model.c.0 != Location::None {
            matches!(ion, FragmentType::c(_))
        } else if model.b.0 != Location::None {
            matches!(ion, FragmentType::b(_))
        } else {
            matches!(ion, FragmentType::a(_))
        };
        let c = if model.z.0 != Location::None {
            matches!(ion, FragmentType::z(_))
        } else if model.y.0 != Location::None {
            matches!(ion, FragmentType::y(_))
        } else {
            matches!(ion, FragmentType::x(_))
        };
        n || c
    }
    let data: SpectrumGraphData = spectrum
        .spectrum()
        .map(|point| {
            let distance = fragments
                .iter()
                .filter(|frag| filter(model, &frag.ion))
                .fold(
                    (
                        (
                            f64::MAX,
                            Fragment::new(
                                MolecularFormula::default(),
                                Charge::new::<e>(0.0),
                                0,
                                FragmentType::precursor,
                                String::new(),
                            ),
                        ),
                        (
                            f64::MAX,
                            Fragment::new(
                                MolecularFormula::default(),
                                Charge::new::<e>(0.0),
                                0,
                                FragmentType::precursor,
                                String::new(),
                            ),
                        ),
                    ),
                    |acc, frag: &Fragment| {
                        let mass_over_charge = frag.mz(MassMode::Monoisotopic);
                        let rel = ((mass_over_charge - point.experimental_mz) / mass_over_charge
                            * MassOverCharge::new::<rustyms::system::mz>(1e6))
                        .value;
                        let abs = (mass_over_charge - point.experimental_mz).value;
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
                    },
                );
            (
                point.annotation.clone(),
                distance.0,                           // rel (ppm)
                distance.1,                           // abs (Da)
                point.experimental_mz,                // mz
                point.experimental_mz * point.charge, // mass
                point.intensity.0,                    // intensity
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
    let n = if model.c.0 != Location::None {
        "c"
    } else if model.b.0 != Location::None {
        "b"
    } else {
        "a"
    };
    let c = if model.z.0 != Location::None {
        "z"
    } else if model.y.0 != Location::None {
        "y"
    } else {
        "x"
    };
    (data, bounds, format!("{n}/{c}"))
}

fn spectrum_graph(
    output: &mut String,
    boundaries: &Boundaries,
    data: &SpectrumGraphData,
    ions: &str,
    x_max: f64,
) {
    write!(output, "<div class='spectrum-graph-y-axis'>").unwrap();
    write!(output, "<span class='max'>{:.2}</span>", boundaries.2).unwrap();
    write!(
        output,
        "<span class='title abs'>Absolute distance to closest {ions} ion (Da)</span><span class='title rel'>Relative distance to closest {ions} ion (ppm)</span>"
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
    let mut first_peptide_index = None;
    for annotation in annotations {
        output.push(annotation.ion.label().to_string());
        output.push(format!("p{}", annotation.peptide_index));
        if let Some(num) = first_peptide_index {
            if num != annotation.peptide_index {
                output.push("mp".to_string());
            }
        } else {
            first_peptide_index = Some(annotation.peptide_index);
        }
        if let Some(pos) = annotation.ion.position() {
            output.push(format!(
                "p{}-{}",
                annotation.peptide_index, pos.sequence_index
            ));
        }
        if annotation.ion.glycan_position().is_some() {
            output.push("glycan".to_string());
        }
        if matches!(annotation.ion, FragmentType::InternalGlycan(_)) {
            output.push("glycan".to_string());
            output.push("internal_glycan".to_string());
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
        output.join(" ")
    }
}

fn get_label(annotations: &[Fragment], multiple_peptides: bool, multiple_glycans: bool) -> String {
    if annotations.is_empty() {
        String::new()
    } else {
        let mut shared_charge = Some(annotations[0].charge);
        let mut shared_ion = Some(annotations[0].ion.label());
        let mut shared_pos = Some(annotations[0].ion.position_label());
        let mut shared_peptide = Some(annotations[0].peptide_index);
        let mut shared_glycan = Some(annotations[0].ion.glycan_position().map(|g| g.attachment()));
        let mut shared_loss = Some(annotations[0].neutral_loss.clone());
        for a in annotations {
            if let Some(charge) = shared_charge {
                if charge != a.charge {
                    shared_charge = None;
                }
            }
            if let Some(ion) = shared_ion {
                if ion != a.ion.label() {
                    shared_ion = None;
                }
            }
            if let Some(pos) = &shared_pos {
                if *pos != a.ion.position_label() {
                    shared_pos = None;
                }
            }
            if let Some(peptide) = shared_peptide {
                if peptide != a.peptide_index {
                    shared_peptide = None;
                }
            }
            if let Some(glycan) = &shared_glycan {
                if *glycan != a.ion.glycan_position().map(|g| g.attachment()) {
                    shared_glycan = None;
                }
            }
            if let Some(loss) = &shared_loss {
                if loss != &a.neutral_loss {
                    shared_loss = None;
                }
            }
        }

        if shared_charge.is_none()
            && shared_ion.is_none()
            && shared_pos.is_none()
            && shared_peptide.is_none()
            && shared_glycan.is_none()
            && shared_loss.is_none()
        {
            "*".to_string()
        } else {
            let charge_str = shared_charge
                .map(|charge| format!("{:+}", charge.value))
                .unwrap_or("*".to_string());
            let ion_str = shared_ion.unwrap_or("*");
            let pos_str = shared_pos
                .map(|pos| pos.unwrap_or_default())
                .unwrap_or("*".to_string());
            let peptide_str = shared_peptide
                .map(|pep| (pep + 1).to_string())
                .unwrap_or("*".to_string());
            let glycan_str = shared_glycan
                .unwrap_or(Some("*".to_string()))
                .unwrap_or_default();
            let loss_str = shared_loss
                .map(|o| o.map(|n| n.hill_notation_html()))
                .unwrap_or(Some("*".to_string()))
                .unwrap_or_default();

            let multi = if annotations.len() > 1 {
                let mut multi = String::new();
                for annotation in annotations {
                    if let FragmentType::InternalGlycan(breakages) = &annotation.ion {
                        let ch = format!("{:+}", annotation.charge.value);
                        write!(
                            multi,
                            "<span>{}<sup>{:+}</sup><sub style='margin-left:-{}ch;min-width:{2}ch;display:inline-block;'>{}{}</sub>{}</span>",
                            breakages
                                .iter()
                                .filter(|b| !matches!(b, GlycanBreakPos::End(_)))
                                .map(|b| format!(
                                    "{}<sub>{}</sub>",
                                    b.label(),
                                    b.position().label()
                                ))
                                .join(""),
                            ch,
                            ch.len(),
                            if multiple_glycans {
                                breakages[0].position().attachment()
                            } else {
                                String::new()
                            },
                            if multiple_peptides {
                                format!("p{}", annotation.peptide_index + 1)
                            } else {
                                String::new()
                            },
                            annotation.neutral_loss.as_ref().map(|n| n.hill_notation_html()).unwrap_or_default(),
                        )
                        .unwrap();
                    } else {
                        let ch = format!("{:+}", annotation.charge.value);
                        write!(
                            multi,
                            "<span>{}<sup>{}</sup><sub style='margin-left:-{}ch;min-width:{2}ch;display:inline-block;'>{}{}{}</sub>{}</span>",
                            annotation.ion.label(),
                            ch,
                            ch.len(),
                            annotation.ion.position_label().unwrap_or(String::new()),
                            if multiple_glycans {
                                annotation.ion.glycan_position().map(|g| g.attachment()).unwrap_or(String::new())
                            } else {
                                String::new()
                            },
                            if multiple_peptides {
                                format!("p{}", annotation.peptide_index + 1)
                            } else {
                                String::new()
                            },
                            annotation.neutral_loss.as_ref().map(|n| n.hill_notation_html()).unwrap_or_default(),
                        )
                        .unwrap();
                    }
                }
                format!("<span class='multi'>{multi}</span>")
            } else {
                String::new()
            };
            let single_internal_glycan =
                matches!(annotations[0].ion, FragmentType::InternalGlycan(_))
                    && annotations.len() == 1;

            if single_internal_glycan {
                if let FragmentType::InternalGlycan(breakages) = &annotations[0].ion {
                    format!(
                        "<span>{}<sup>{}</sup><sub style='margin-left:-{}ch;min-width:{2}ch;display:inline-block;'>{}{}</sub>{}</span>",
                        breakages
                            .iter()
                            .filter(|b| !matches!(b, GlycanBreakPos::End(_)))
                            .map(|b| format!("{}<sub>{}</sub>", b.label(), b.position().label()))
                            .join(""),
                        charge_str,
                        charge_str.len(),
                        if multiple_glycans {
                            breakages[0].position().attachment()
                        } else {
                            String::new()
                        },
                        if multiple_peptides {
                            format!("p{}", peptide_str)
                        } else {
                            String::new()
                        },
                        annotations[0]
                            .neutral_loss.as_ref()
                            .map(|n| n.hill_notation_html())
                            .unwrap_or_default(),
                    )
                } else {
                    unreachable!();
                }
            } else {
                format!(
                    "<span>{}<sup>{}</sup><sub style='margin-left:-{}ch;min-width:{2}ch;display:inline-block;'>{}{}{}</sub>{}</span>{}",
                    ion_str,
                    charge_str,
                    charge_str.len(),
                    pos_str,
                    if multiple_glycans {
                        glycan_str
                    } else {
                        String::new()
                    },
                    if multiple_peptides {
                        format!("p{}", peptide_str)
                    } else {
                        String::new()
                    },
                    loss_str,
                    multi
                )
            }
        }
    }
}

type PositionCoverage = Vec<Vec<HashSet<FragmentType>>>;

pub struct Limits {
    pub mz: MassOverCharge,
    pub intensity: f64,
    pub intensity_unassigned: f64,
}

fn get_overview(spectrum: &AnnotatedSpectrum) -> (Limits, PositionCoverage) {
    let mut output: PositionCoverage = spectrum
        .peptide
        .peptides()
        .iter()
        .map(|p| vec![HashSet::new(); p.sequence.len()])
        .collect();
    let mut max_mz: MassOverCharge = MassOverCharge::new::<mz>(0.0);
    let mut max_intensity: f64 = 0.0;
    let mut max_intensity_unassigned: f64 = 0.0;
    for peak in spectrum.spectrum() {
        max_mz = max_mz.max(peak.experimental_mz);
        max_intensity_unassigned = max_intensity_unassigned.max(peak.intensity.0);
        if !peak.annotation.is_empty() {
            max_intensity = max_intensity.max(peak.intensity.0);
            peak.annotation.iter().for_each(|frag| {
                frag.ion
                    .position()
                    .map(|i| output[frag.peptide_index][i.sequence_index].insert(frag.ion.clone()));
            });
        }
    }
    (
        Limits {
            mz: max_mz * 1.01,
            intensity: max_intensity * 1.01,
            intensity_unassigned: max_intensity_unassigned * 1.01,
        },
        output,
    )
}

fn render_peptide(
    output: &mut String,
    spectrum: &AnnotatedSpectrum,
    overview: PositionCoverage,
    multiple_peptides: bool,
) {
    write!(output, "<div class='complex-peptide'>").unwrap();
    match &spectrum.peptide {
        ComplexPeptide::Singular(peptide) => {
            render_linear_peptide(output, peptide, &overview[0], 0, multiple_peptides)
        }
        ComplexPeptide::Multimeric(peptides) => {
            for (index, peptide) in peptides.iter().enumerate() {
                render_linear_peptide(output, peptide, &overview[index], index, multiple_peptides);
            }
        }
        _ => panic!("Unhandled complex peptide type"),
    }
    write!(output, "</div>").unwrap();
}

fn render_linear_peptide(
    output: &mut String,
    peptide: &LinearPeptide,
    overview: &[HashSet<FragmentType>],
    peptide_index: usize,
    multiple_peptides: bool,
) {
    write!(output, "<div class='peptide'>").unwrap();
    if multiple_peptides {
        write!(output, "<span class='name'>{}</span>", peptide_index + 1).unwrap();
    }
    if peptide.n_term.is_some() {
        write!(output, "<span class='modification term'></span>").unwrap();
    }
    for (index, (pos, ions)) in peptide.sequence.iter().zip(overview).enumerate() {
        let mut classes = String::new();
        if !pos.modifications.is_empty() {
            write!(classes, " modification").unwrap();
        }
        if !pos.possible_modifications.is_empty() {
            write!(classes, " possible-modification").unwrap();
        }
        if !classes.is_empty() {
            classes = format!(" class='{}'", classes.trim());
        }
        write!(
            output,
            "<span data-pos='{}-{}'{classes} tabindex='0' title='N terminal position: {}, C terminal position: {}'>{}",
            peptide_index,
            index,
            index + 1,
            peptide.sequence.len() - index,
            pos.aminoacid.char()
        )
        .unwrap();
        for ion in ions {
            write!(
                output,
                "<span class='corner {}'></span>",
                ion.label().trim_end_matches('·')
            )
            .unwrap();
        }
        write!(output, "</span>").unwrap();
    }
    if peptide.c_term.is_some() {
        write!(output, "<span class='modification term'></span>").unwrap();
    }
    if let Some(charge_carriers) = &peptide.charge_carriers {
        // TODO: When rustyms is updated use the display impl here
        write!(output, "<span class='charge-carriers'>/",).unwrap();
        write!(
            output,
            "{}",
            charge_carriers
                .charge_carriers
                .iter()
                .map(|c| c.0)
                .sum::<isize>()
        )
        .unwrap();
        if !charge_carriers.charge_carriers.iter().all(|c| {
            c.1 == MolecularFormula::new(&[(Element::H, None, 1), (Element::Electron, None, -1)])
                .unwrap()
        }) {
            write!(output, "[").unwrap();
            let mut first = true;
            for (amount, formula) in &charge_carriers.charge_carriers {
                if first {
                    first = false;
                } else {
                    write!(output, ",").unwrap();
                }
                let electron_index = formula
                    .elements()
                    .iter()
                    .position(|el| el.0 == Element::Electron);
                let charge = electron_index.map(|ei| match -formula.elements()[ei].2 {
                    1 => "+".to_string(),
                    -1 => "-".to_string(),
                    n => n.to_string(),
                });
                if let (Some(electron_index), Some(charge), 2) =
                    (electron_index, &charge, formula.elements().len())
                {
                    let element_index = 1 - electron_index;
                    write!(
                        output,
                        "{}{}{}",
                        amount,
                        formula.elements()[element_index].0,
                        charge
                    )
                    .unwrap();
                } else {
                    write!(
                        output,
                        "{}({}){}",
                        amount,
                        formula,
                        charge.unwrap_or_default()
                    )
                    .unwrap();
                }
            }
            write!(output, "]").unwrap();
        }
        write!(output, "</span>").unwrap();
    }
    write!(output, "</div>").unwrap();
}

fn render_spectrum(
    output: &mut String,
    spectrum: &AnnotatedSpectrum,
    fragments: &[Fragment],
    boundaries: &Boundaries,
    limits: &Limits,
    selection: &str,
    multiple_peptides: bool,
    multiple_glycans: bool,
) {
    write!(
        output,
        "<div class='canvas-wrapper label' aria-hidden='true' style='--min-mz:0;--max-mz:{0};--max-intensity:{1};--y-max:{2};--y-min:{3};--abs-max-initial:{2};--abs-min-initial:{3};--rel-max-initial:{4};--rel-min-initial:{5};' data-initial-max-mz='{0}' data-initial-max-intensity='{1}' data-initial-max-intensity-assigned='{6}'>",
        limits.mz.value, limits.intensity_unassigned, boundaries.2, boundaries.3, boundaries.0, boundaries.1, limits.intensity
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
        limits.intensity_unassigned / 4.0
    )
    .unwrap();
    write!(
        output,
        "<span class='tick'><span>{:.2}</span></span>",
        limits.intensity_unassigned / 2.0
    )
    .unwrap();
    write!(
        output,
        "<span class='tick'><span>{:.2}</span></span>",
        3.0 * limits.intensity_unassigned / 4.0
    )
    .unwrap();
    write!(
        output,
        "<span class='tick'><span class='last'>{:.2}</span></span>",
        limits.intensity_unassigned
    )
    .unwrap();
    write!(output, "</div>").unwrap();
    write!(output, "<div class='canvas canvas-spectrum'>").unwrap();
    write!(
        output,
        "<span class='selection {selection}' hidden='true'></span>"
    )
    .unwrap();

    for peak in spectrum.spectrum() {
        write!(
            output,
            "<span class='peak {}' style='--mz:{};--intensity:{};' data-label='{}' {}>{}</span>",
            get_classes(&peak.annotation),
            peak.experimental_mz.value,
            peak.intensity,
            (peak.experimental_mz.value * 10.0).round() / 10.0,
            if peak.intensity.0 / limits.intensity >= 0.1 {
                "data-show-label='true'"
            } else {
                ""
            },
            get_label(&peak.annotation, multiple_peptides, multiple_glycans),
        )
        .unwrap();
    }
    for peak in fragments {
        write!(
            output,
            "<span class='theoretical peak {}' style='--mz:{};' data-label='{}'>{}</span>",
            get_classes(&[peak.clone()]),
            peak.mz(MassMode::Monoisotopic).value,
            (peak.mz(MassMode::Monoisotopic).value * 10.0).round() / 10.0,
            get_label(&[peak.clone()], multiple_peptides, multiple_glycans),
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
        limits.mz.value / 4.0
    )
    .unwrap();
    write!(
        output,
        "<span class='tick'><span>{:.2}</span></span>",
        limits.mz.value / 2.0
    )
    .unwrap();
    write!(
        output,
        "<span class='tick'><span>{:.2}</span></span>",
        3.0 * limits.mz.value / 4.0
    )
    .unwrap();
    write!(
        output,
        "<span class='tick'><span class='last'>{:.2}</span></span>",
        limits.mz.value
    )
    .unwrap();
    write!(output, "</div>").unwrap();
}

pub fn spectrum_table(
    spectrum: &AnnotatedSpectrum,
    fragments: &[Fragment],
    multiple_peptides: bool,
) -> String {
    fn generate_text(annotation: &Fragment) -> (String, String, String) {
        if let Some(pos) = annotation.ion.position() {
            (
                (pos.sequence_index + 1).to_string(),
                pos.series_number.to_string(),
                annotation.ion.label().to_string(),
            )
        } else if let Some(pos) = annotation.ion.glycan_position() {
            (
                pos.attachment(),
                format!("{}{}", pos.series_number, pos.branch_names()),
                annotation.ion.label().to_string(),
            )
        } else if let FragmentType::InternalGlycan(breakages) = &annotation.ion {
            (
                breakages
                    .first()
                    .map(|b| b.position().attachment())
                    .unwrap_or("-".to_string()),
                "-".to_string(),
                breakages
                    .iter()
                    .filter(|b| !matches!(b, GlycanBreakPos::End(_)))
                    .map(|b| format!("{}<sub>{}</sub>", b.label(), b.position().label()))
                    .join(""),
            )
        } else {
            // precursor
            (
                "-".to_string(),
                "-".to_string(),
                annotation.ion.label().to_string(),
            )
        }
    }
    let mut output = String::new();
    write!(
        output,
        "<label class='show-unassigned'><input type='checkbox' switch/>Show background peaks</label>
        <label class='show-matched'><input type='checkbox' switch checked/>Show annotated peaks</label>
        <label class='show-missing-fragments'><input type='checkbox' switch/>Show missing fragments</label>
        <table id='spectrum-table' class='wide-table'>
            <thead><tr>
                {}
                <th>Position</th>
                <th>Ion type</th>
                <th>Loss</th>
                <th>Intensity</th>
                <th>mz Theoretical</th>
                <th>mz Error (Th)</th>
                <th>mz Error (ppm)</th>
                <th>Charge</th>
                <th>Series Number</th>
                <th>Additional label</th>
            </tr></thead><tdata>",
        if multiple_peptides {
            "<th>Peptide</th>"
        } else {
            ""
        }
    )
    .unwrap();
    // class followed by all data
    let mut data = Vec::new();
    for peak in spectrum.spectrum() {
        if peak.annotation.is_empty() {
            data.push((
                peak.experimental_mz.value,
                [
                    "unassigned".to_string(),
                    "-".to_string(),
                    "-".to_string(),
                    "-".to_string(),
                    "-".to_string(),
                    format!("{:.2}", peak.intensity),
                    format!("{:.2}", peak.experimental_mz.value),
                    "-".to_string(),
                    "-".to_string(),
                    "-".to_string(),
                    "-".to_string(),
                    "-".to_string(),
                ],
            ));
        } else {
            for annotation in &peak.annotation {
                let (sequence_index, series_number, label) = generate_text(annotation);
                data.push((
                    peak.experimental_mz.value,
                    [
                        "matched".to_string(),
                        if multiple_peptides {
                            format!("{}", annotation.peptide_index + 1)
                        } else {
                            String::new()
                        },
                        sequence_index.to_string(),
                        label.to_string(),
                        annotation
                            .neutral_loss
                            .as_ref()
                            .map_or(String::new(), |v| v.to_string()),
                        format!("{:.2}", peak.intensity),
                        format!("{:.2}", peak.experimental_mz.value),
                        format!(
                            "{:.5}",
                            (annotation.mz(MassMode::Monoisotopic) - peak.experimental_mz)
                                .abs()
                                .value
                        ),
                        format!(
                            "{:.2}",
                            ((annotation.mz(MassMode::Monoisotopic) - peak.experimental_mz).abs()
                                / annotation.mz(MassMode::Monoisotopic)
                                * 1e6)
                                .value
                        ),
                        format!("{:+}", annotation.charge.value),
                        series_number,
                        annotation.label.clone(),
                    ],
                ));
            }
        }
    }
    for fragment in fragments {
        if !spectrum
            .spectrum()
            .any(|p| p.annotation.iter().any(|a| a == fragment))
        {
            let (sequence_index, series_number, label) = generate_text(fragment);
            data.push((
                fragment.mz(MassMode::Monoisotopic).value,
                [
                    "fragment".to_string(),
                    if multiple_peptides {
                        format!("{}", fragment.peptide_index + 1)
                    } else {
                        String::new()
                    },
                    sequence_index.to_string(),
                    label.to_string(),
                    fragment
                        .neutral_loss
                        .as_ref()
                        .map_or(String::new(), |v| v.to_string()),
                    "-".to_string(),
                    format!("{:.2}", fragment.mz(MassMode::Monoisotopic).value),
                    "-".to_string(),
                    "-".to_string(),
                    format!("{:+}", fragment.charge.value),
                    series_number,
                    fragment.label.clone(),
                ],
            ))
        }
    }
    data.sort_unstable_by(|a, b| a.0.total_cmp(&b.0));
    for row in data {
        write!(output, "<tr class='{}'>", row.1[0]).unwrap();
        for cell in &row.1[if multiple_peptides { 1 } else { 2 }..] {
            write!(output, "<td>{}</td>", cell).unwrap();
        }
        write!(output, "</tr>").unwrap();
    }
    write!(output, "</tdata></table>").unwrap();
    output
}

#[allow(non_snake_case)]
#[derive(Debug, Default, Clone)]
struct IonStats<T> {
    a: T,
    b: T,
    c: T,
    d: T,
    v: T,
    w: T,
    x: T,
    y: T,
    z: T,
    precursor: T,
    oxonium: T,
    Y: T,
    peaks: T,
}

impl<T> IonStats<T>
where
    T: std::ops::Add<T, Output = T> + std::ops::AddAssign<T> + Copy + Default + PartialOrd<T>,
{
    fn add(
        mut self,
        annotation: &[Fragment],
        peptide_index: Option<usize>,
        value: T,
        allow_double_counting: bool,
    ) -> Self {
        let mut seen_a = false;
        let mut seen_b = false;
        let mut seen_c = false;
        let mut seen_d = false;
        let mut seen_v = false;
        let mut seen_w = false;
        let mut seen_x = false;
        let mut seen_y = false;
        let mut seen_z = false;
        let mut seen_precursor = false;
        let mut seen_oxonium = false;
        #[allow(non_snake_case)]
        let mut seen_Y = false;
        for fragment in annotation {
            if peptide_index.map_or(false, |index| index != fragment.peptide_index) {
                continue; // Skip any annotation not for this peptide
            }
            match fragment.ion {
                FragmentType::a(..) if !seen_a || allow_double_counting => {
                    self.a += value;
                    seen_a = true;
                }
                FragmentType::b(..) if !seen_b || allow_double_counting => {
                    self.b += value;
                    seen_b = true;
                }
                FragmentType::c(..) if !seen_c || allow_double_counting => {
                    self.c += value;
                    seen_c = true;
                }
                FragmentType::d(..) if !seen_d || allow_double_counting => {
                    self.d += value;
                    seen_d = true;
                }
                FragmentType::v(..) if !seen_v || allow_double_counting => {
                    self.v += value;
                    seen_v = true;
                }
                FragmentType::w(..) if !seen_w || allow_double_counting => {
                    self.w += value;
                    seen_w = true;
                }
                FragmentType::x(..) if !seen_x || allow_double_counting => {
                    self.x += value;
                    seen_x = true;
                }
                FragmentType::y(..) if !seen_y || allow_double_counting => {
                    self.y += value;
                    seen_y = true;
                }
                FragmentType::z(..) | FragmentType::z·(..) if !seen_z || allow_double_counting => {
                    self.z += value;
                    seen_z = true;
                }
                FragmentType::B(..) | FragmentType::InternalGlycan(..)
                    if !seen_oxonium || allow_double_counting =>
                {
                    self.oxonium += value;
                    seen_oxonium = true;
                }
                FragmentType::Y(..) if !seen_Y || allow_double_counting => {
                    self.Y += value;
                    seen_Y = true;
                }
                FragmentType::precursor if !seen_precursor || allow_double_counting => {
                    self.precursor += value;
                    seen_precursor = true;
                }
                _ => (),
            }
        }
        self.peaks += value;
        self
    }

    /// Call the callback on all nonzero lines (not the total line) in the `IonStats`.
    /// Get the name, value and the singular value provided.
    fn map<U: Copy, R>(&self, callback: impl Fn(&str, T, U) -> R, single: U) -> Vec<R> {
        let mut result = Vec::new();
        if self.a > T::default() {
            result.push(callback("a", self.a, single));
        }
        if self.b > T::default() {
            result.push(callback("b", self.b, single));
        }
        if self.c > T::default() {
            result.push(callback("c", self.c, single));
        }
        if self.d > T::default() {
            result.push(callback("d", self.d, single));
        }
        if self.v > T::default() {
            result.push(callback("v", self.v, single));
        }
        if self.w > T::default() {
            result.push(callback("w", self.w, single));
        }
        if self.x > T::default() {
            result.push(callback("x", self.x, single));
        }
        if self.y > T::default() {
            result.push(callback("y", self.y, single));
        }
        if self.z > T::default() {
            result.push(callback("z", self.z, single));
        }
        if self.precursor > T::default() {
            result.push(callback("precursor", self.precursor, single));
        }
        if self.oxonium > T::default() {
            result.push(callback("oxonium", self.oxonium, single));
        }
        if self.Y > T::default() {
            result.push(callback("Y", self.Y, single));
        }
        result
    }
    /// Call the callback on all nonzero lines (not the total line) in the `IonStats`.
    /// Get the name, value and the singular value provided.
    fn zip<U, R>(&self, callback: impl Fn(&str, T, U) -> R, other: IonStats<U>) -> Vec<R> {
        let mut result = Vec::new();
        if self.a > T::default() {
            result.push(callback("a", self.a, other.a));
        }
        if self.b > T::default() {
            result.push(callback("b", self.b, other.b));
        }
        if self.c > T::default() {
            result.push(callback("c", self.c, other.c));
        }
        if self.d > T::default() {
            result.push(callback("d", self.d, other.d));
        }
        if self.v > T::default() {
            result.push(callback("v", self.v, other.v));
        }
        if self.w > T::default() {
            result.push(callback("w", self.w, other.w));
        }
        if self.x > T::default() {
            result.push(callback("x", self.x, other.x));
        }
        if self.y > T::default() {
            result.push(callback("y", self.y, other.y));
        }
        if self.z > T::default() {
            result.push(callback("z", self.z, other.z));
        }
        if self.precursor > T::default() {
            result.push(callback("precursor", self.precursor, other.precursor));
        }
        if self.oxonium > T::default() {
            result.push(callback("oxonium", self.oxonium, other.oxonium));
        }
        if self.Y > T::default() {
            result.push(callback("Y", self.Y, other.Y));
        }
        result
    }
}

fn general_stats(output: &mut String, spectrum: &AnnotatedSpectrum, fragments: &[Fragment]) {
    fn format(found: usize, total: usize) -> String {
        format!(
            "{:.2}% ({}/{})",
            found as f64 / total as f64 * 100.0,
            found,
            total
        )
    }
    fn format_f64(found: f64, total: f64) -> String {
        format!(
            "{:.2}% ({:.2e}/{:.2e})",
            found / total * 100.0,
            found,
            total
        )
    }

    let mut mass_row = String::new();
    let mut fragments_row = String::new();
    let mut fragments_details_row = String::new();
    let mut peaks_row = String::new();
    let mut intensity_row = String::new();
    let mut intensity_details_row = String::new();
    let mut positions_row = String::new();

    let total_intensity: f64 = spectrum.spectrum().map(|p| p.intensity.0).sum();
    for (peptide_index, peptide) in spectrum.peptide.peptides().iter().enumerate() {
        let precursor = peptide
            .formulas()
            .iter()
            .map(|f| {
                format!(
                    "{} | avg: {}",
                    display_mass(f.monoisotopic_mass()),
                    display_mass(f.average_weight())
                )
            })
            .join(", ");
        let (num_annotated, intensity_annotated) = spectrum
            .spectrum()
            .filter(|p| {
                p.annotation
                    .iter()
                    .any(|a| a.peptide_index == peptide_index)
            })
            .fold(
                (IonStats::default(), IonStats::default()),
                |(n, intensity), p| {
                    (
                        n.add(&p.annotation, Some(peptide_index), 1, true),
                        intensity.add(&p.annotation, Some(peptide_index), p.intensity.0, false),
                    )
                },
            );
        let total_fragments = fragments
            .iter()
            .filter(|f| f.peptide_index == peptide_index)
            .fold(IonStats::default(), |n, p| {
                n.add(&[p.clone()], Some(peptide_index), 1, true)
            });

        let percentage_positions_covered = spectrum
            .spectrum()
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
            "<td>{}</td>",
            format(num_annotated.peaks, total_fragments.peaks)
        )
        .unwrap();

        write!(
            fragments_details_row,
            "<td><table>{}</table></td>",
            total_fragments
                .zip(
                    |name, total, annotated| {
                        format!(
                            "<tr><td>{}</td><td>{}</td></tr>",
                            name,
                            format(annotated, total)
                        )
                    },
                    num_annotated.clone(),
                )
                .join("")
        )
        .unwrap();

        write!(
            peaks_row,
            "<td>{}</td>",
            format(num_annotated.peaks, spectrum.spectrum().len())
        )
        .unwrap();
        write!(
            intensity_row,
            "<td>{}</td>",
            format_f64(intensity_annotated.peaks, total_intensity),
        )
        .unwrap();
        write!(
            intensity_details_row,
            "<td><table>{}</table></td>",
            intensity_annotated
                .map(
                    |name, annotated, total| {
                        format!(
                            "<tr><td>{}</td><td>{}</td></tr>",
                            name,
                            format_f64(annotated, total)
                        )
                    },
                    total_intensity,
                )
                .join("")
        )
        .unwrap();
        write!(positions_row, "<td>{percentage_positions_covered:.2}%</td>").unwrap();
    }

    write!(output, "<label><input type='checkbox' switch id='general-stats-show-details'>Show detailed fragment breakdown</label><table class='general-stats'>").unwrap();
    if spectrum.peptide.peptides().len() > 1 {
        write!(output, "<tr><td>General Statistics</td>").unwrap();
        for i in 0..spectrum.peptide.peptides().len() {
            write!(output, "<td>Peptide {}</td>", i + 1).unwrap();
        }
        write!(output, "<td>Combined</td></tr>").unwrap();
        // Add a combined stats column
        let (num_annotated, intensity_annotated) = spectrum
            .spectrum()
            .filter(|p| !p.annotation.is_empty())
            .fold(
                (IonStats::default(), IonStats::default()),
                |(n, intensity), p| {
                    (
                        n.add(&p.annotation, None, 1, true),
                        intensity.add(&p.annotation, None, p.intensity.0, false),
                    )
                },
            );
        let total_fragments = fragments.iter().fold(IonStats::default(), |n, f| {
            n.add(&[f.clone()], None, 1, true)
        });

        write!(mass_row, "<td>-</td>").unwrap();
        write!(
            fragments_row,
            "<td>{}</td>",
            format(num_annotated.peaks, total_fragments.peaks)
        )
        .unwrap();
        write!(
            fragments_details_row,
            "<td><table>{}</table></td>",
            total_fragments
                .zip(
                    |name, total, annotated| {
                        format!(
                            "<tr><td>{}</td><td>{}</td></tr>",
                            name,
                            format(annotated, total)
                        )
                    },
                    num_annotated.clone(),
                )
                .join("")
        )
        .unwrap();
        write!(
            peaks_row,
            "<td>{}</td>",
            format(num_annotated.peaks, spectrum.spectrum().len())
        )
        .unwrap();
        write!(
            intensity_row,
            "<td>{}</td>",
            format_f64(intensity_annotated.peaks, total_intensity)
        )
        .unwrap();
        write!(
            intensity_details_row,
            "<td><table>{}</table></td>",
            intensity_annotated
                .map(
                    |name, annotated, total| {
                        format!(
                            "<tr><td>{}</td><td>{}</td></tr>",
                            name,
                            format_f64(annotated, total)
                        )
                    },
                    total_intensity,
                )
                .join("")
        )
        .unwrap();
        write!(positions_row, "<td>-</td>").unwrap();
    }
    write!(
        output,
        "<tr><td>Precursor Mass (M)</td>{mass_row}</tr>
        <tr><td>Fragments found</td>{fragments_row}</tr>
        <tr class='fragments-detail'><td>Fragments detailed</td>{fragments_details_row}</tr>
        <tr><td>Peaks annotated</td>{peaks_row}</tr>
        <tr><td>Intensity annotated</td>{intensity_row}</tr>
        <tr class='fragments-detail'><td>Intensity detailed</td>{intensity_details_row}</tr>
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

pub fn display_mass(value: Mass) -> HtmlElement {
    let (num, suf, full) = engineering_notation(value.value, 3);
    let suf = suf.map_or(String::new(), |suf| suf.to_string());
    HtmlElement::new(HtmlTag::span)
        .class("mass")
        .header(
            "title",
            if suf.is_empty() {
                format!("{} Da", value.value)
            } else {
                format!("{} {}Da\n{} Da", full, suf, value.value)
            },
        )
        .content(format!("{} {}Da", num, suf))
}

/// Display the given value in engineering notation eg `1000` -> `10 k`, with the given number of decimal points and returns the suffix separately.
/// A value of `0.0` will result in the lowest possible suffix `0.0 q`.
fn engineering_notation(value: f64, precision: usize) -> (String, Option<char>, f64) {
    const BIG_SUFFIXES: &[char] = &[' ', 'k', 'M', 'G', 'T', 'P', 'E', 'Z', 'Y', 'R', 'Q'];
    const SMALL_SUFFIXES: &[char] = &[' ', 'm', 'μ', 'n', 'p', 'f', 'a', 'z', 'y', 'r', 'q'];
    let base = if value == 0.0 {
        0
    } else {
        ((value.abs().log10() / 3.0).floor() as isize).clamp(
            -(SMALL_SUFFIXES.len() as isize - 1),
            BIG_SUFFIXES.len() as isize - 1,
        )
    };
    match base.cmp(&0) {
        Ordering::Less => (
            format!(
                "{:.precision$}",
                value * 10_usize.pow(base.unsigned_abs() as u32 * 3) as f64,
            ),
            Some(SMALL_SUFFIXES[base.unsigned_abs()]),
            value * 10_usize.pow(base.unsigned_abs() as u32 * 3) as f64,
        ),
        Ordering::Equal => (format!("{value:.precision$}"), None, value),
        Ordering::Greater => (
            format!(
                "{:.precision$}",
                value / 10_usize.pow(base as u32 * 3) as f64,
            ),
            Some(BIG_SUFFIXES[base as usize]),
            value / 10_usize.pow(base as u32 * 3) as f64,
        ),
    }
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
