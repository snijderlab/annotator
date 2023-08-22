use std::{collections::HashSet, fmt::Write};

use rustyms::{AnnotatedSpectrum, Charge, Fragment, FragmentType, Mass, MassOverCharge, Zero};

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
            fragment.mz().value,
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
    let mut output = "<div class='spectrum' onload='SpectrumSetUp()'>".to_string();
    let (limits, overview) = get_overview(spectrum);
    let (graph_data, graph_boundaries) = spectrum_graph_boundaries(spectrum, fragments);

    spectrum_top_buttons(&mut output, id, &limits, &graph_boundaries).unwrap();

    write!(output, "<div class='wrapper unassigned'>").unwrap();
    create_ion_legend(&mut output, &format!("{id}-1"));
    render_peptide(&mut output, spectrum, overview);
    render_spectrum(&mut output, spectrum, &graph_boundaries, limits, "first");
    // Spectrum graph
    spectrum_graph(&mut output, &graph_boundaries, &graph_data, limits.0.value);
    write!(output, "</div></div>").unwrap();
    // Spectrum table
    collapsible(
        &mut output,
        &format!("{id}-table-1"),
        "Peaks table".to_string(),
        spectrum_table(spectrum, &format!("{id}-table-1")),
    );

    write!(output, "</div>").unwrap();
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
    write!(
        output,
        "<div class='settings spectrum-graph-setup'><p>Spectrum graph</p>"
    )?;
    spectrum_graph_header(output, spectrum_graph_boundaries)?;
    write!(output, "</div>")?;
    write!(
        output,
        "<div class='settings render-setup'><p>Render setup</p>"
    )?;
    write!(output, "<label for='{id}-width'>Width</label>")?;
    write!(
        output,
        "<input id='{id}-width' class='width' type='text' value='100%'/>",
    )?;
    write!(output, "<label for='{id}-height'>Height</label>")?;
    write!(
        output,
        "<input id='{id}-height' class='height' type='text' value='250px'/>",
    )?;
    write!(
        output,
        "<label for='{id}-fs-peptide'>Peptide font size</label>"
    )?;
    write!(
        output,
        "<input id='{id}-fs-peptide' class='fs-peptide' type='text' value='1.25rem'/>",
    )?;
    write!(
        output,
        "<label for='{id}-peptide-stroke'>Peptide stroke width</label>"
    )?;
    write!(
        output,
        "<input id='{id}-peptide-stroke' class='stroke-peptide' type='text' value='2px'/>",
    )?;
    write!(
        output,
        "<label for='{id}-fs-spectrum'>Spectrum font size</label>"
    )?;
    write!(
        output,
        "<input id='{id}-fs-spectrum' class='fs-spectrum' type='text' value='1rem'/>",
    )?;
    write!(
        output,
        "<label for='{id}-spectrum-stroke'>Spectrum stroke width</label>"
    )?;
    write!(
        output,
        "<input id='{id}-spectrum-stroke' class='stroke-spectrum' type='text' value='2px'/>",
    )?;
    write!(
        output,
        "<input id='{id}-compact' class='compact' type='checkbox'/>",
    )?;
    write!(output, "<label for='{id}-compact'>Compact peptide</label>")?;
    write!(output, "</div>")?;
    Ok(())
}

type Boundaries = (f64, f64, f64, f64, f64, f64, f64, f64, f64, f64);
type SpectrumGraphData = Vec<(
    Option<Fragment>,
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
                            Fragment::new(Mass::zero(), Charge::zero(), FragmentType::precursor),
                        ),
                        (
                            f64::MAX,
                            Fragment::new(Mass::zero(), Charge::zero(), FragmentType::precursor),
                        ),
                    ),
                    |acc, frag: &Fragment| {
                        let rel = ((frag.mz() - point.experimental_mz) / frag.mz()
                            * MassOverCharge::new::<rustyms::mz>(1e6))
                        .value;
                        let abs = (frag.mz() - point.experimental_mz).value;
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
                point
                    .annotation
                    .as_ref()
                    .map(|a| Some(a.clone()))
                    .unwrap_or(None),
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
            "<span class='point {}' style='--rel:{};--abs:{};--mz:{};--intensity:{}' data-ppm='{}' data-abs='{}' data-mz='{}' data-intensity='{}' data-reference-fragment-rel='{}' data-reference-fragment-abs='{}' data-pos='{}'></span>",
            point.0.as_ref().map(|a| a.ion.to_string()).unwrap_or("unassigned".to_string()),
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
            point.0.as_ref().and_then(|a| a.ion.position().map(|p|p.sequence_index + 1)).unwrap_or(point.2.1.ion.position().map(|p| p.sequence_index + 1).unwrap_or(0))
        )
        .unwrap();
    }
    write!(output, "</div>").unwrap();
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
            f.ion
                .position()
                .map(|i| output[i.sequence_index].insert(f.ion));
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
            "<span data-pos='{}'{classes} tabindex='0' title='N terminal position: {}, C terminal position: {}'>{}",
            index + 1,
            index + 1,
            spectrum.peptide.sequence.len() - index,
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
    write!(output, "<div class='zoom-out' tabindex='0'>Zoom Out</div>").unwrap();

    for peak in &spectrum.spectrum {
        match &peak.annotation {
            None => {
                write!(
                    output,
                    "<span class='peak unassigned' style='--mz:{};--intensity:{};' data-label='{}'></span>",
                    peak.experimental_mz.value, peak.intensity, (peak.experimental_mz.value * 10.0).round() / 10.0
                )
                .unwrap();
            }
            Some(f) => {
                write!(
                    output,
                    "<span class='peak {} label' style='--mz:{};--intensity:{};' data-pos='{}' data-label='{}'>",
                    f.ion,
                    peak.experimental_mz.value,
                    peak.intensity,
                    f.ion
                        .position()
                        .map_or("*".to_string(), |i| (i.sequence_index + 1).to_string()),
                    (peak.experimental_mz.value * 10.0).round() / 10.0,
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
                        .position()
                        .map_or("*".to_string(), |i| i.series_number.to_string())
                )
                .unwrap();
            }
        }
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
                    .position()
                    .map_or("*".to_string(), |i| (i.sequence_index + 1).to_string()),
                f.ion,
                peak.intensity,
                peak.experimental_mz.value,
                (f.mz() - peak.experimental_mz).abs().value,
                ((f.mz() - peak.experimental_mz).abs() / f.mz() * 1e6).value,
                f.charge.value,
                f.ion
                    .position()
                    .map_or("*".to_string(), |i| i.series_number.to_string()),
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
