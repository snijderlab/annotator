use std::{collections::HashMap, fmt::Write};

use crate::{Theme, html_builder::HtmlTag};
use itertools::Itertools;
use mzannotate::fragment::*;
use mzcore::{
    prelude::*,
    sequence::{
        CrossLinkName, FlankingSequence, GnoComposition, Linked, Modification,
        SimpleModificationInner,
    },
};
use ordered_float::OrderedFloat;

use super::render_full_glycan;

/// Render a peptidoform with the given ions present. It returns a lookup with the unique ids for all peptides.
pub fn render_peptide(
    output: &mut String,
    compound_peptidoform: &CompoundPeptidoformIon,
    overview: Option<super::PositionCoverage>,
    local_confidence: Option<Vec<Vec<Vec<f64>>>>,
    flanking_sequences: Option<(&FlankingSequence, &FlankingSequence)>,
    glycan_footnotes: &mut Vec<String>,
    theme: Theme,
) -> Vec<(usize, usize)> {
    let mut unique_peptide_lookup = Vec::new();
    let multiple_peptidoforms = compound_peptidoform.peptidoform_ions().len() > 1;
    let multiple_peptides = compound_peptidoform
        .peptidoform_ions()
        .iter()
        .any(|p| p.peptidoforms().len() > 1);
    let max_position_intensity = overview
        .as_ref()
        .and_then(|o| {
            o.iter()
                .flat_map(|o| o.iter())
                .flat_map(|o| o.iter())
                .flat_map(|o| o.values())
                .max()
                .copied()
        })
        .unwrap_or_default();
    let max_combined_position_intensity = overview
        .as_ref()
        .and_then(|o| {
            o.iter()
                .flat_map(|o| o.iter())
                .flat_map(|o| o.iter())
                .map(|o| o.values().sum::<OrderedFloat<f32>>())
                .max()
        })
        .unwrap_or_default();
    write!(
        output,
        "<div class='complex-peptide{}{}' style='--max-position-intensity:{max_position_intensity};--max-combined-position-intensity:{max_combined_position_intensity};'>",
        overview.as_ref().map_or("", |_| " with-fragments"),
        local_confidence
            .as_ref()
            .map_or("", |_| " with-local-confidence")
    )
    .unwrap();
    let mut cross_link_lookup = Vec::new();
    for (peptidoform_ion_index, peptidoform_ion) in
        compound_peptidoform.peptidoform_ions().iter().enumerate()
    {
        for (peptidoform_index, peptidoform) in peptidoform_ion.peptidoforms().iter().enumerate() {
            render_linear_peptidoform(
                output,
                peptidoform,
                overview
                    .as_ref()
                    .map(|o| o[peptidoform_ion_index][peptidoform_index].as_slice()),
                local_confidence
                    .as_ref()
                    .map(|o| o[peptidoform_ion_index][peptidoform_index].as_slice()),
                peptidoform_ion_index,
                peptidoform_index,
                unique_peptide_lookup.len(),
                multiple_peptidoforms,
                multiple_peptides,
                flanking_sequences.filter(|_| !multiple_peptidoforms && !multiple_peptides),
                &mut cross_link_lookup,
                glycan_footnotes,
                theme,
            );
            unique_peptide_lookup.push((peptidoform_ion_index, peptidoform_index));
        }
    }
    write!(output, "</div>").unwrap();
    unique_peptide_lookup
}

fn render_linear_peptidoform(
    output: &mut String,
    peptidoform: &Peptidoform<Linked>,
    overview: Option<&[HashMap<FragmentType, OrderedFloat<f32>>]>,
    local_confidence: Option<&[f64]>,
    peptidoform_ion_index: usize,
    peptidoform_index: usize,
    unique_peptidoform_index: usize,
    multiple_peptidoform_ions: bool,
    multiple_peptidoforms: bool,
    flanking_sequences: Option<(&FlankingSequence, &FlankingSequence)>,
    cross_link_lookup: &mut Vec<CrossLinkName>,
    glycan_footnotes: &mut Vec<String>,
    theme: Theme,
) {
    write!(
        output,
        "<div class='peptide pu{unique_peptidoform_index} p{peptidoform_ion_index}'>"
    )
    .unwrap();
    if multiple_peptidoforms || multiple_peptidoform_ions {
        write!(
            output,
            "{}",
            HtmlTag::span
                .new()
                .class("name")
                .data([(
                    "pos",
                    format!("{peptidoform_ion_index}-{peptidoform_index}")
                )])
                .content(if multiple_peptidoform_ions && multiple_peptidoforms {
                    format!("{}.{}", peptidoform_ion_index + 1, peptidoform_index + 1)
                } else if multiple_peptidoform_ions {
                    (peptidoform_ion_index + 1).to_string()
                } else if multiple_peptidoforms {
                    (peptidoform_index + 1).to_string()
                } else {
                    String::new()
                })
        )
        .unwrap();
    }
    if let Some((n_flanking, _)) = flanking_sequences {
        let (symbol, title) = match n_flanking {
            FlankingSequence::Unknown => ("?".to_string(), "Flanking sequence unknown"),
            FlankingSequence::Terminal => ("-".to_string(), "Flanks the terminus"),
            FlankingSequence::AminoAcid(aa) => (aa.to_string(), "Flanking sequence"),
            FlankingSequence::Sequence(seq) => (seq.to_string(), "Flanking sequence"),
        };
        write!(
            output,
            "<span class='flanking' title='{title}'>{symbol}</span>",
        )
        .unwrap();
    }
    if !peptidoform.get_n_term().is_empty() {
        let mut possible_modifications = Vec::new();
        let mut cross_link = Vec::new();
        let mut modifications = Vec::new();
        let mut glycans = String::new();
        for m in peptidoform.get_n_term() {
            match m {
                Modification::Ambiguous {
                    group,
                    modification,
                    localisation_score,
                    ..
                } => possible_modifications.push(format!(
                    "{modification}\x23{group}{}",
                    localisation_score.map_or(String::new(), |s| format!("({s:3})"))
                )),
                Modification::CrossLink { peptide, name, .. } => {
                    let xl_index = cross_link_lookup
                        .iter()
                        .position(|xl| *xl == *name)
                        .unwrap_or_else(|| {
                            cross_link_lookup.push(name.clone());
                            cross_link_lookup.len() - 1
                        });
                    cross_link.push((
                        format!("c{xl_index}"),
                        format!("xl.{}p{peptide}", xl_index + 1),
                        format!("x{}", xl_index + 1),
                        format!("{name}"),
                    ));
                }
                Modification::Simple(m) => match &**m {
                    SimpleModificationInner::Gno {
                        composition: GnoComposition::Topology(structure),
                        id,
                        ..
                    } => {
                        let svg = render_full_glycan(
                            structure,
                            false,
                            true,
                            theme,
                            glycan_footnotes,
                            false,
                            false,
                            0,
                            0,
                        );
                        glycans.push_str(&svg);
                        modifications.push(id.to_string());
                    }
                    SimpleModificationInner::GlycanStructure(structure) => {
                        let svg = render_full_glycan(
                            structure,
                            false,
                            true,
                            theme,
                            glycan_footnotes,
                            false,
                            false,
                            0,
                            0,
                        );
                        glycans.push_str(&svg);
                    }
                    other => modifications.push(other.to_string()),
                },
            }
        }
        write!(
            output,
            "<span class='{} term' data-cross-links='{}' data-cross-links-compact='{}' title='N-term {}{}{}{}{}{}{}{}'>{}</span>", 
            (!possible_modifications.is_empty()).then_some("possible-modification").into_iter()
                .chain((!cross_link.is_empty()).then_some("cross-link").into_iter())
                .chain((!modifications.is_empty()).then_some("modification").into_iter())
                .chain(cross_link.iter().map(|c| c.0.as_str())).join(" "),
                cross_link.iter().map(|c| c.1.as_str()).join(" "),
                cross_link.iter().map(|c| c.2.as_str()).join(" "),
                if modifications.is_empty() {""} else {"Modification: "}, 
                modifications.join(", "),
                if !modifications.is_empty() && !possible_modifications.is_empty() {", "} else {""},
                if possible_modifications.is_empty() {""} else{"Modification of unknown position: "}, 
                possible_modifications.join(", "),
                if (!modifications.is_empty() || !possible_modifications.is_empty()) && !cross_link.is_empty() {", "} else {""},
                if cross_link.is_empty() {""} else {"Cross-link: "}, 
                cross_link.iter().map(|c| c.3.as_str()).join(", "),
                if !glycans.is_empty() {
                    format!("<div class='glycans'>{glycans}</div>")
                } else {String::new()},
        )
        .unwrap();
    }
    for (index, ((pos, ions), confidence)) in peptidoform
        .sequence()
        .iter()
        .zip(
            overview
                .map(|o| o.iter().map(Some).collect_vec())
                .unwrap_or(vec![None; peptidoform.len()]),
        )
        .zip(
            local_confidence
                .map(|o| o.iter().map(Some).collect_vec())
                .unwrap_or(vec![None; peptidoform.len()]),
        )
        .enumerate()
    {
        let mut classes = String::new();
        let mut xl_indices = Vec::new();
        let mut xl_names = Vec::new();
        let mut xl_peptides = Vec::new();
        let mut modification = false;
        let mut modifications = String::new();
        let mut glycans = String::new();
        let mut modifications_of_unknown_position = String::new();
        for m in &pos.modifications {
            match m {
                Modification::Ambiguous {
                    group,
                    modification,
                    localisation_score,
                    ..
                } => {
                    write!(
                        modifications_of_unknown_position,
                        "{}{modification}\x23{group}{}",
                        if !modifications_of_unknown_position.is_empty() {
                            ", "
                        } else {
                            ""
                        },
                        localisation_score.map_or(String::new(), |s| format!("({s:3})"))
                    )
                    .unwrap();
                }
                Modification::CrossLink { peptide, name, .. } => {
                    let xl_index = cross_link_lookup
                        .iter()
                        .position(|xl| *xl == *name)
                        .unwrap_or_else(|| {
                            cross_link_lookup.push(name.clone());
                            cross_link_lookup.len() - 1
                        });
                    write!(classes, " c{xl_index}").unwrap();
                    xl_indices.push(xl_index + 1);
                    xl_names.push(name);
                    if *peptide != peptidoform_index {
                        xl_peptides.push(*peptide + 1);
                    }
                }
                Modification::Simple(m) => match &**m {
                    SimpleModificationInner::Gno {
                        composition: GnoComposition::Topology(structure),
                        id,
                        ..
                    } => {
                        let svg = render_full_glycan(
                            structure,
                            false,
                            true,
                            theme,
                            glycan_footnotes,
                            false,
                            false,
                            0,
                            0,
                        );
                        glycans.push_str(&svg);
                        write!(
                            modifications,
                            "{}{id}",
                            if !modifications.is_empty() { ", " } else { "" },
                        )
                        .unwrap()
                    }
                    SimpleModificationInner::GlycanStructure(structure) => {
                        let svg = render_full_glycan(
                            structure,
                            false,
                            true,
                            theme,
                            glycan_footnotes,
                            false,
                            false,
                            0,
                            0,
                        );
                        glycans.push_str(&svg);
                    }
                    other => {
                        modification = true;
                        write!(
                            modifications,
                            "{}{other}",
                            if !modifications.is_empty() { ", " } else { "" },
                        )
                        .unwrap()
                    }
                },
            }
        }

        let cross_links = format!(
            "xl.{}{}",
            xl_indices.iter().join(","),
            if xl_peptides.is_empty() {
                String::new()
            } else {
                format!("p{}", xl_peptides.iter().join(","))
            }
        );

        let cross_links_compact = format!("x{}", xl_indices.iter().join(","),);

        if modification {
            write!(classes, " modification").unwrap();
        }
        if !xl_indices.is_empty() {
            write!(classes, " cross-link").unwrap();
        }
        if pos.modifications.iter().any(|m| m.is_ambiguous()) {
            write!(classes, " possible-modification").unwrap();
        }
        if !classes.is_empty() {
            classes = format!(" class='{}'", classes.trim());
        }
        write!(
            output,
            "<span data-pos='{peptidoform_ion_index}-{peptidoform_index}-{index}' data-cross-links='{cross_links}' data-cross-links-compact='{cross_links_compact}'{classes} tabindex='0' title='N terminal position: {}, C terminal position: {}{}{}{}' style='--confidence:{}'>{}",
            index + 1,
            peptidoform.len() - index,
            if modifications.is_empty() {
                String::new()
            } else {
                format!(", Modifications: {modifications}")
            },
            if modifications_of_unknown_position.is_empty() {
                String::new()
            } else {
                format!(", Modifications of unknown position: {modifications_of_unknown_position}")
            },
            if xl_names.is_empty() {
                String::new()
            } else {
                format!(", Cross-link{}: {}", if xl_names.len() == 1{""} else {"s"}, xl_names.iter().join(", "))
            },
            confidence.copied().unwrap_or_default(),
            pos.aminoacid.pro_forma_definition(),
        )
        .unwrap();
        if !glycans.is_empty() {
            write!(output, "<div class='glycans'>{glycans}</div>").unwrap();
        }
        if let Some(ions) = ions {
            for (ion, intensity) in ions {
                if !matches!(
                    ion,
                    FragmentType::Immonium(_, _)
                        | FragmentType::PrecursorSideChainLoss(_, _)
                        | FragmentType::Diagnostic(_)
                ) {
                    write!(
                        output,
                        "<span class='corner {}' style='--intensity:{intensity}'></span>",
                        ion.kind()
                    )
                    .unwrap();
                }
            }
        }
        write!(output, "</span>").unwrap();
    }
    if !peptidoform.get_c_term().is_empty() {
        let mut possible_modifications = Vec::new();
        let mut cross_link = Vec::new();
        let mut modifications = Vec::new();
        let mut glycans = String::new();
        for m in peptidoform.get_c_term() {
            match m {
                Modification::Ambiguous {
                    group,
                    modification,
                    localisation_score,
                    ..
                } => possible_modifications.push(format!(
                    "{modification}\x23{group}{}",
                    localisation_score.map_or(String::new(), |s| format!("({s:3})"))
                )),
                Modification::CrossLink { peptide, name, .. } => {
                    let xl_index = cross_link_lookup
                        .iter()
                        .position(|xl| *xl == *name)
                        .unwrap_or_else(|| {
                            cross_link_lookup.push(name.clone());
                            cross_link_lookup.len() - 1
                        });
                    cross_link.push((
                        format!("c{xl_index}"),
                        format!("xl.{}p{peptide}", xl_index + 1),
                        format!("x{}", xl_index + 1),
                        format!("{name}"),
                    ));
                }
                Modification::Simple(m) => match &**m {
                    SimpleModificationInner::Gno {
                        composition: GnoComposition::Topology(structure),
                        id,
                        ..
                    } => {
                        let svg = render_full_glycan(
                            structure,
                            false,
                            true,
                            theme,
                            glycan_footnotes,
                            false,
                            false,
                            0,
                            0,
                        );
                        glycans.push_str(&svg);
                        modifications.push(id.to_string());
                    }
                    SimpleModificationInner::GlycanStructure(structure) => {
                        let svg = render_full_glycan(
                            structure,
                            false,
                            true,
                            theme,
                            glycan_footnotes,
                            false,
                            false,
                            0,
                            0,
                        );
                        glycans.push_str(&svg);
                    }
                    other => modifications.push(other.to_string()),
                },
            }
        }
        write!(
            output,
            "<span class='{} term' data-cross-links='{}' data-cross-links-compact='{}' title='C-term {}{}{}{}{}{}{}{}'>{}</span>", 
            (!possible_modifications.is_empty()).then_some("possible-modification").into_iter()
                .chain((!cross_link.is_empty()).then_some("cross-link").into_iter())
                .chain((!modifications.is_empty()).then_some("modification").into_iter())
                .chain(cross_link.iter().map(|c| c.0.as_str())).join(" "),
                cross_link.iter().map(|c| c.1.as_str()).join(" "),
                cross_link.iter().map(|c| c.2.as_str()).join(" "),
                if modifications.is_empty() {""} else {"Modification: "}, 
                modifications.join(", "),
                if !modifications.is_empty() && !possible_modifications.is_empty() {", "} else {""},
                if possible_modifications.is_empty() {""} else{"Modification of unknown position: "}, 
                possible_modifications.join(", "),
                if (!modifications.is_empty() || !possible_modifications.is_empty()) && !cross_link.is_empty() {", "} else {""},
                if cross_link.is_empty() {""} else {"Cross-link: "}, 
                cross_link.iter().map(|c| c.3.as_str()).join(", "),
                if !glycans.is_empty() {
                    format!("<div class='glycans'>{glycans}</div>")
                } else {String::new()},
        )
        .unwrap();
    }
    if let Some(charge_carriers) = &peptidoform.get_charge_carriers() {
        write!(
            output,
            "<span class='charge-carriers'>/{charge_carriers}</span>",
        )
        .unwrap();
    }
    if let Some((_, c_flanking_)) = flanking_sequences {
        write!(
            output,
            "<span class='flanking'>{}</span>",
            match c_flanking_ {
                FlankingSequence::Unknown => "?".to_string(),
                FlankingSequence::Terminal => "Terminus".to_string(),
                FlankingSequence::AminoAcid(aa) => aa.to_string(),
                FlankingSequence::Sequence(seq) => seq.to_string(),
            }
        )
        .unwrap();
    }
    write!(output, "</div>").unwrap();
}
