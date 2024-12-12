use std::{collections::HashSet, fmt::Write};

use itertools::Itertools;
use rustyms::{
    fragment::*, modification::CrossLinkName, CompoundPeptidoform, LinearPeptide, Linked,
    Modification,
};

use crate::html_builder::HtmlTag;

/// Render a peptidoform with the given ions present. It returns a lookup with the unique ids for all peptides.
pub fn render_peptide(
    output: &mut String,
    compound_peptidoform: &CompoundPeptidoform,
    overview: Option<super::PositionCoverage>,
    local_confidence: Option<Vec<Vec<Vec<f64>>>>,
) -> Vec<(usize, usize)> {
    let mut unique_peptide_lookup = Vec::new();
    let multiple_peptidoforms = compound_peptidoform.peptidoforms().len() > 1;
    let multiple_peptides = compound_peptidoform
        .peptidoforms()
        .iter()
        .any(|p| p.peptides().len() > 1);
    write!(
        output,
        "<div class='complex-peptide{}{}'>",
        overview.as_ref().map_or("", |_| " with-fragments"),
        local_confidence
            .as_ref()
            .map_or("", |_| " with-local-confidence")
    )
    .unwrap();
    let mut cross_link_lookup = Vec::new();
    for (peptidoform_index, peptidoform) in compound_peptidoform.peptidoforms().iter().enumerate() {
        for (peptide_index, peptide) in peptidoform.peptides().iter().enumerate() {
            render_linear_peptide(
                output,
                peptide,
                overview
                    .as_ref()
                    .map(|o| o[peptidoform_index][peptide_index].as_slice()),
                local_confidence
                    .as_ref()
                    .map(|o| o[peptidoform_index][peptide_index].as_slice()),
                peptidoform_index,
                peptide_index,
                unique_peptide_lookup.len(),
                multiple_peptidoforms,
                multiple_peptides,
                &mut cross_link_lookup,
            );
            unique_peptide_lookup.push((peptidoform_index, peptide_index));
        }
    }
    write!(output, "</div>").unwrap();
    unique_peptide_lookup
}

fn render_linear_peptide(
    output: &mut String,
    peptide: &LinearPeptide<Linked>,
    overview: Option<&[HashSet<FragmentType>]>,
    local_confidence: Option<&[f64]>,
    peptidoform_index: usize,
    peptide_index: usize,
    unique_peptide_index: usize,
    multiple_peptidoforms: bool,
    multiple_peptides: bool,
    cross_link_lookup: &mut Vec<CrossLinkName>,
) {
    write!(
        output,
        "<div class='peptide pu{unique_peptide_index} p{peptidoform_index}'>"
    )
    .unwrap();
    if multiple_peptides || multiple_peptidoforms {
        write!(
            output,
            "{}",
            HtmlTag::span
                .new()
                .class("name")
                .data([("pos", format!("{peptidoform_index}-{peptide_index}"))])
                .content(if multiple_peptidoforms && multiple_peptides {
                    format!("{}.{}", peptidoform_index + 1, peptide_index + 1)
                } else if multiple_peptidoforms {
                    (peptidoform_index + 1).to_string()
                } else if multiple_peptides {
                    (peptide_index + 1).to_string()
                } else {
                    String::new()
                })
        )
        .unwrap();
    }
    if peptide.get_n_term().is_some() {
        let (class, xl_long, xl_compact) = match peptide.get_n_term() {
            Some(Modification::Ambiguous { .. }) => (
                "possible-modification".to_string(),
                String::new(),
                String::new(),
            ),
            Some(Modification::CrossLink { peptide, name, .. }) => {
                let xl_index = cross_link_lookup
                    .iter()
                    .position(|xl| *xl == *name)
                    .unwrap_or_else(|| {
                        cross_link_lookup.push(name.clone());
                        cross_link_lookup.len() - 1
                    });
                (
                    format!("cross-link c{xl_index}"),
                    format!("xl.{}p{peptide}", xl_index + 1),
                    format!("x{}", xl_index + 1),
                )
            }
            Some(Modification::Simple(_)) => {
                ("modification".to_string(), String::new(), String::new())
            }
            None => (String::new(), String::new(), String::new()),
        };
        write!(
            output,
            "<span class='{class} term' data-cross-links='{xl_long}' data-cross-links-compact='{xl_compact}'></span>",
        )
        .unwrap();
    }
    for (index, ((pos, ions), confidence)) in peptide
        .sequence()
        .iter()
        .zip(
            overview
                .map(|o| o.iter().map(Some).collect_vec())
                .unwrap_or(vec![None; peptide.len()]),
        )
        .zip(
            local_confidence
                .map(|o| o.iter().map(Some).collect_vec())
                .unwrap_or(vec![None; peptide.len()]),
        )
        .enumerate()
    {
        let mut classes = String::new();
        let mut xl_indices = Vec::new();
        let mut xl_names = Vec::new();
        let mut xl_peptides = Vec::new();
        for m in &pos.modifications {
            if let Modification::CrossLink { peptide, name, .. } = m {
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
                if *peptide != peptide_index {
                    xl_peptides.push(*peptide + 1);
                }
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

        if pos
            .modifications
            .iter()
            .any(|m| matches!(m, Modification::Simple(_)))
        {
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
            "<span data-pos='{peptidoform_index}-{peptide_index}-{index}' data-cross-links='{cross_links}' data-cross-links-compact='{cross_links_compact}'{classes} tabindex='0' title='N terminal position: {}, C terminal position: {}{}' style='--confidence:{}'>{}",
            index + 1,
            peptide.len() - index,
            if xl_names.is_empty() {
                String::new()
            } else {
                format!(", Cross-link{}: {}", if xl_names.len() == 1{""} else {"s"}, xl_names.iter().join(", "))
            },
            confidence.copied().unwrap_or_default(),
            pos.aminoacid.char(),
        )
        .unwrap();
        if let Some(ions) = ions {
            for ion in ions {
                if !matches!(
                    ion,
                    FragmentType::Immonium(_, _)
                        | FragmentType::PrecursorSideChainLoss(_, _)
                        | FragmentType::Diagnostic(_)
                ) {
                    write!(
                        output,
                        "<span class='corner {}'></span>",
                        ion.label().trim_end_matches('·')
                    )
                    .unwrap();
                }
            }
        }
        write!(output, "</span>").unwrap();
    }
    if peptide.get_c_term().is_some() {
        let (class, xl_long, xl_compact) = match peptide.get_c_term() {
            Some(Modification::Ambiguous { .. }) => (
                "possible-modification".to_string(),
                String::new(),
                String::new(),
            ),
            Some(Modification::CrossLink { peptide, name, .. }) => {
                let xl_index = cross_link_lookup
                    .iter()
                    .position(|xl| *xl == *name)
                    .unwrap_or_else(|| {
                        cross_link_lookup.push(name.clone());
                        cross_link_lookup.len() - 1
                    });
                (
                    format!("cross-link c{xl_index}"),
                    format!("xl.{}p{peptide}", xl_index + 1),
                    format!("x{}", xl_index + 1),
                )
            }
            Some(Modification::Simple(_)) => {
                ("modification".to_string(), String::new(), String::new())
            }
            None => (String::new(), String::new(), String::new()),
        };
        write!(
            output,
            "<span class='{class} term' data-cross-links='{xl_long}' data-cross-links-compact='{xl_compact}'></span>",
        )
        .unwrap();
    }
    if let Some(charge_carriers) = &peptide.get_charge_carriers() {
        write!(
            output,
            "<span class='charge-carriers'>/{charge_carriers}</span>",
        )
        .unwrap();
    }
    write!(output, "</div>").unwrap();
}
