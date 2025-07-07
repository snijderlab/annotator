use std::fmt::Write;

use crate::{
    metadata_render::OptionalString,
    render::{display_formula, display_neutral_loss, label::display_sequence_index},
};
use itertools::Itertools;
use rustyms::{
    annotation::AnnotatedSpectrum, fragment::FragmentType, prelude::*, spectrum::PeakSpectrum,
};

pub fn spectrum_table(
    spectrum: &AnnotatedSpectrum,
    fragments: &[Fragment],
    multiple_peptidoform_ions: bool,
    multiple_peptidoforms: bool,
) -> String {
    fn generate_text(
        annotation: &Fragment,
        compound_peptidoform: &CompoundPeptidoformIon,
    ) -> (String, String, String) {
        let format_label = |label: (Option<String>, std::borrow::Cow<'_, str>)| {
            if let Some(sup) = label.0 {
                format!("<sup>{sup}</sup>{}", label.1)
            } else {
                label.1.to_string()
            }
        };
        let (sequence_index, series_number, label) = if let Some(pos) = annotation.ion.position() {
            (
                format!(
                    "{}{}",
                    compound_peptidoform.peptidoform_ions()
                        [annotation.peptidoform_ion_index.unwrap_or_default()]
                    .peptidoforms()[annotation.peptidoform_index.unwrap_or_default()]
                        [pos.sequence_index]
                        .aminoacid
                        .pro_forma_definition(),
                    display_sequence_index(pos.sequence_index)
                ),
                pos.series_number.to_string(),
                format_label(annotation.ion.label()),
            )
        } else if let Some(pos) = annotation.ion.glycan_position() {
            (
                pos.attachment(),
                format!("{}{}", pos.series_number, pos.branch_names()),
                annotation.ion.kind().to_string(),
            )
        } else if let FragmentType::B { b, y, .. } = &annotation.ion {
            (
                b.attachment.map_or("-".to_string(), |_| b.attachment()),
                format!("B<sub>{}</sub>", b.label())
                    + &y.iter()
                        .map(|b| format!("Y<sub>{}</sub>", b.label()))
                        .join(""),
                "B".to_string(),
            )
        } else if let FragmentType::Y(bonds) = &annotation.ion {
            (
                bonds
                    .first()
                    .map(|b| b.attachment())
                    .unwrap_or("-".to_string()),
                bonds.iter().map(|b| b.label()).join(""),
                "Y".to_string(),
            )
        } else if let FragmentType::Immonium(pos, aa) = &annotation.ion {
            (
                pos.map_or("-".to_string(), |p| {
                    display_sequence_index(p.sequence_index)
                }),
                pos.map_or("-".to_string(), |p| p.series_number.to_string()),
                format!(
                    "immonium {}",
                    aa.aminoacid
                        .one_letter_code()
                        .map(|c| c.to_string())
                        .or_else(|| aa.aminoacid.three_letter_code().map(|c| c.to_string()))
                        .unwrap_or_else(|| aa.aminoacid.to_string())
                ),
            )
        } else {
            // precursor
            (
                "-".to_string(),
                "-".to_string(),
                format_label(annotation.ion.label()),
            )
        };
        (
            sequence_index,
            series_number,
            format!(
                "<span title='mzPAF: {}'>{label}</span>",
                annotation.to_mzPAF()
            ),
        )
    }
    let mut output = String::new();
    write!(
        output,
        "<p id='export-sequence' title='The inputted peptide written as fully compatible with the ProForma spec. So removes things like custom modifications.'>Universal ProForma definition: <span>{}</span></p>
        <label class='show-unassigned'><input type='checkbox' switch/>Show background peaks</label>
        <label class='show-matched'><input type='checkbox' switch checked/>Show annotated peaks</label>
        <label class='show-missing-fragments'><input type='checkbox' switch/>Show missing fragments</label>
        <table id='spectrum-table' class='wide-table'>
            <thead><tr>
                {}
                {}
                <th>Position</th>
                <th>Ion type</th>
                <th>Loss</th>
                <th>Intensity</th>
                <th title='Experimental mz for background and matched peaks, but theoretical mz for missing peaks'>mz</th>
                <th>Formula</th>
                <th>mz Error (Th)</th>
                <th>mz Error (ppm)</th>
                <th>Charge</th>
                <th>Series Number</th>
                <th>Additional label</th>
            </tr></thead><tdata>",
        spectrum.peptide,
        if multiple_peptidoform_ions {
            "<th>Peptidoform</th>"
        } else {
            ""
        },
        if multiple_peptidoforms {
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
                    "<span title='?'>-</span>".to_string(),
                    "-".to_string(),
                    format!("<span title='{0}'>{0:.2}</span>", peak.intensity),
                    format!(
                        "<span title='{0}'>{0:.2}</span>",
                        peak.experimental_mz.value
                    ),
                    "-".to_string(),
                    "-".to_string(),
                    "-".to_string(),
                    "-".to_string(),
                    "-".to_string(),
                    "-".to_string(),
                ],
            ));
        } else {
            for annotation in &peak.annotation {
                let (sequence_index, series_number, label) =
                    generate_text(annotation, &spectrum.peptide);
                data.push((
                    peak.experimental_mz.value,
                    [
                        "matched".to_string(),
                        if multiple_peptidoform_ions {
                            annotation
                                .peptidoform_ion_index
                                .map_or("?".to_string(), |i| (i + 1).to_string())
                        } else {
                            String::new()
                        },
                        if multiple_peptidoforms {
                            annotation
                                .peptidoform_index
                                .map_or("?".to_string(), |i| (i + 1).to_string())
                        } else {
                            String::new()
                        },
                        sequence_index.to_string(),
                        label.to_string(),
                        annotation
                            .neutral_loss
                            .iter()
                            .map(|n| display_neutral_loss(n, true))
                            .join(""),
                        format!("<span title='{0}'>{0:.2}</span>", peak.intensity),
                        format!(
                            "<span title='{0}'>{0:.2}</span>",
                            peak.experimental_mz.value
                        ),
                        annotation
                            .formula
                            .as_ref()
                            .map(|f| display_formula(f, true))
                            .to_optional_string(),
                        annotation
                            .mz(MassMode::Monoisotopic)
                            .map_or("-".to_string(), |f| {
                                format!("{:.5}", (f - peak.experimental_mz).abs().value)
                            }),
                        annotation
                            .mz(MassMode::Monoisotopic)
                            .map_or("-".to_string(), |f| {
                                format!("{:.2}", f.ppm(peak.experimental_mz).value * 1e6)
                            }),
                        format!("{:+}", annotation.charge.value),
                        series_number,
                        annotation.formula.as_ref().map_or(String::new(), |f| {
                            f.labels().iter().map(|l| l.to_string()).join(",")
                        }),
                    ],
                ));
            }
        }
    }
    for fragment in fragments {
        if !spectrum.spectrum().any(|p| p.annotation.contains(fragment)) {
            let (sequence_index, series_number, label) =
                generate_text(&fragment.clone(), &spectrum.peptide);
            data.push((
                fragment
                    .mz(MassMode::Monoisotopic)
                    .map_or(f64::NAN, |v| v.value),
                [
                    "fragment".to_string(),
                    if multiple_peptidoform_ions {
                        fragment
                            .peptidoform_index
                            .map_or("?".to_string(), |i| (i + 1).to_string())
                    } else {
                        String::new()
                    },
                    if multiple_peptidoforms {
                        fragment
                            .peptidoform_index
                            .map_or("?".to_string(), |i| (i + 1).to_string())
                    } else {
                        String::new()
                    },
                    sequence_index.to_string(),
                    label.to_string(),
                    fragment
                        .neutral_loss
                        .iter()
                        .map(|n| display_neutral_loss(n, true))
                        .join(""),
                    "-".to_string(),
                    fragment
                        .mz(MassMode::Monoisotopic)
                        .map_or("-".to_string(), |v| {
                            format!("<span title='{0}'>{0:.2}</span>", v.value)
                        }),
                    fragment
                        .formula
                        .as_ref()
                        .map(|f| display_formula(f, true))
                        .to_optional_string(),
                    "-".to_string(),
                    "-".to_string(),
                    format!("{:+}", fragment.charge.value),
                    series_number,
                    fragment.formula.as_ref().map_or(String::new(), |f| {
                        f.labels().iter().map(|l| l.to_string()).join(",")
                    }),
                ],
            ))
        }
    }
    data.sort_unstable_by(|a, b| a.0.total_cmp(&b.0));
    for row in data {
        write!(output, "<tr class='{}'>", row.1[0]).unwrap();
        for cell in &row.1[if multiple_peptidoforms {
            if multiple_peptidoform_ions { 1 } else { 2 }
        } else {
            3
        }..]
        {
            write!(output, "<td>{cell}</td>").unwrap();
        }
        write!(output, "</tr>").unwrap();
    }
    write!(output, "</tdata></table>").unwrap();
    output
}
