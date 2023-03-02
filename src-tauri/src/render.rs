use std::{collections::HashSet, fmt::Write};

use rustyms::{fragment::FragmentType, spectrum::AnnotatedSpectrum, MassOverCharge, Zero};

use crate::html_builder::{HtmlElement, HtmlTag};

pub fn annotated_spectrum(spectrum: &AnnotatedSpectrum, id: &str) -> String {
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
    write!(output, "<div class='wrapper unassigned'>").unwrap();
    create_ion_legend(&mut output, &format!("{id}-1"));
    render_peptide(&mut output, spectrum, overview);
    render_spectrum(&mut output, spectrum, limits, "first");
    write!(output, "</div>").unwrap();
    // Second spectrum
    collapsible(
        &mut output,
        &format!("{id}-table-1"),
        "Peaks table".to_string(),
        spectrum_table(spectrum, &format!("{id}-table-1")),
    );
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
            f.ion.sequence_index().map(|i| output[i].insert(f.ion));
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
            write!(
                output,
                "<span class='corner {}'></span>",
                ion.to_string().trim_end_matches('·')
            )
            .unwrap();
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
                    f.ion
                        .sequence_index()
                        .map_or("*".to_string(), |i| (i + 1).to_string()),
                )
                .unwrap();
                let ch = format!("{:+}", peak.charge.value);
                write!(
                    output,
                    "<span>{}<sup>{}</sup><sub style='margin-left:-{}ch'>{}</sub></span></span>",
                    f.ion,
                    ch,
                    ch.len(),
                    f.ion
                        .sequence_index()
                        .map_or("*".to_string(), |i| (i + 1).to_string())
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

fn spectrum_table(spectrum: &AnnotatedSpectrum, table_id: &str) -> String {
    let mut output = String::new();
    write!(
        output,
        "<label class='background'><input type='checkbox'/>Show background peaks</label>
        <table id='{table_id}' class='wide-table'>
            <tr>
                <th>Position</th>
                <th>Ion type</th>
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
        match &peak.annotation {
            None => write!(
                output,
                "<tr class='unassigned'>
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
            ),
            Some(f) => write!(
                output,
                "<tr>
                <td>{}</td>
                <td>{}</td>
                <td>{:.2}</td>
                <td>{:.2}</td>
                <td>{:.5}</td>
                <td>{:.2}</td>
                <td>{:+}</td>
                <td>{}</td>
            </tr>",
                f.ion
                    .sequence_index()
                    .map_or("*".to_string(), |i| (i + 1).to_string()),
                f.ion,
                peak.intensity,
                peak.experimental_mz.value,
                (f.mz() - peak.experimental_mz).abs().value,
                ((f.mz() - peak.experimental_mz).abs() / f.mz() * 1e6).value,
                f.charge.value,
                0
            ),
        }
        .unwrap()
    }
    write!(output, "</table>").unwrap();
    output
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
