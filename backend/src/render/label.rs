use itertools::Itertools;
use rustyms::{
    fragment::{FragmentType, GlycanBreakPos},
    AmbiguousLabel, CompoundPeptidoform, Fragment, SequencePosition,
};
use std::fmt::Write;

pub fn get_label(
    compound_peptidoform: &CompoundPeptidoform,
    annotations: &[Fragment],
    multiple_peptidoforms: bool,
    multiple_peptides: bool,
    multiple_glycans: bool,
) -> String {
    if annotations.is_empty() {
        String::new()
    } else {
        let mut shared_charge = Some(annotations[0].charge);
        let mut shared_ion = Some(annotations[0].ion.label());
        let mut shared_pos = Some(annotations[0].ion.position_label());
        let mut shared_peptidoform = Some(annotations[0].peptidoform_index);
        let mut shared_peptide = Some(annotations[0].peptide_index);
        let mut shared_glycan = Some(annotations[0].ion.glycan_position().map(|g| g.attachment()));
        let mut shared_loss = Some(annotations[0].neutral_loss.clone());
        let mut shared_xl = Some(get_xl(&annotations[0]));
        let mut shared_ambiguous_amino_acids = Some(get_ambiguous_amino_acids(
            &annotations[0],
            multiple_peptides,
        ));
        let mut shared_modifications = Some(get_modifications(
            &annotations[0],
            multiple_peptides,
            compound_peptidoform,
        ));
        let mut shared_charge_carriers = Some(get_charge_carriers(&annotations[0]));

        for a in annotations {
            if let Some(charge) = shared_charge {
                if charge != a.charge {
                    shared_charge = None;
                }
            }
            if let Some(ion) = &shared_ion {
                if *ion != a.ion.label() {
                    shared_ion = None;
                }
            }
            if let Some(pos) = &shared_pos {
                if *pos != a.ion.position_label() {
                    shared_pos = None;
                }
            }
            if let Some(peptidoform) = shared_peptidoform {
                if peptidoform != a.peptidoform_index {
                    shared_peptidoform = None;
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
            if let Some(xl) = &shared_xl {
                if *xl != get_xl(a) {
                    shared_xl = None;
                }
            }
            if let Some(aaa) = &shared_ambiguous_amino_acids {
                if *aaa != get_ambiguous_amino_acids(a, multiple_peptides) {
                    shared_ambiguous_amino_acids = None;
                }
            }
            if let Some(sm) = &shared_modifications {
                if *sm != get_modifications(a, multiple_peptides, compound_peptidoform) {
                    shared_modifications = None;
                }
            }
            if let Some(cc) = &shared_charge_carriers {
                if *cc != get_charge_carriers(a) {
                    shared_charge_carriers = None;
                }
            }
        }

        if shared_charge.is_none()
            && shared_ion.is_none()
            && shared_pos.is_none()
            && shared_peptidoform.is_none()
            && shared_peptide.is_none()
            && shared_glycan.is_none()
            && shared_loss.is_none()
            && shared_xl.is_none()
            && shared_ambiguous_amino_acids.is_none()
            && shared_modifications.is_none()
            && shared_charge_carriers.is_none()
        {
            "*".to_string()
        } else {
            let charge_str = shared_charge
                .map(|charge| format!("{:+}", charge.value))
                .unwrap_or("*".to_string());
            let ion_str = shared_ion
                .map(|c| c.into_owned())
                .unwrap_or("*".to_string());
            let pos_str = shared_pos
                .map(|pos| pos.unwrap_or_default())
                .unwrap_or("*".to_string());
            let peptidoform_str = shared_peptidoform
                .map(|pep| (pep + 1).to_string())
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
            let xl_str = shared_xl.unwrap_or("*".to_string());
            let aaa_str = shared_ambiguous_amino_acids.unwrap_or("*".to_string());
            let sm_str = shared_modifications.unwrap_or("*".to_string());
            let cc_str = shared_charge_carriers.unwrap_or("*".to_string());

            let multi = if annotations.len() > 1 {
                let mut multi = String::new();
                for annotation in annotations {
                    write!(
                        multi,
                        "{}",
                        get_single_label(
                            annotation,
                            multiple_peptidoforms,
                            multiple_peptides,
                            multiple_glycans,
                            compound_peptidoform,
                        )
                    )
                    .unwrap();
                }
                format!("<span class='multi'>{multi}</span>")
            } else {
                String::new()
            };
            let single_internal_glycan =
                matches!(annotations[0].ion, FragmentType::Oxonium(_)) && annotations.len() == 1;

            if single_internal_glycan {
                get_single_label(
                    &annotations[0],
                    multiple_peptidoforms,
                    multiple_peptides,
                    multiple_glycans,
                    compound_peptidoform,
                )
            } else {
                format!(
                    "<span>{}<sup class='charge'>{}</sup><sub style='--charge-width:{};'><span class='series'>{}</span><span class='glycan-id'>{}</span><span class='peptide-id'>{}</span></sub><span class='neutral-losses'>{}</span><span class='cross-links'>{}</span><span class='ambiguous-amino-acids'>{}</span><span class='modifications'>{}</span><span class='charge-carriers'>{}</span></span>{}",
                    ion_str,
                    charge_str,
                    charge_str.len(),
                    pos_str,
                    if multiple_glycans {
                        glycan_str
                    } else {
                        String::new()
                    },
                    if multiple_peptidoforms && multiple_peptides {
                        format!("p{}.{}", peptidoform_str, peptide_str)
                    }else if multiple_peptidoforms {
                        format!("p{}", peptidoform_str)
                    } else if multiple_peptides {
                        format!("p{}", peptide_str)
                    } else {
                        String::new()
                    },
                    loss_str,
                    xl_str,
                    aaa_str,
                    sm_str,
                    cc_str,
                    multi,
                )
            }
        }
    }
}

fn get_single_label(
    annotation: &Fragment,
    multiple_peptidoforms: bool,
    multiple_peptides: bool,
    multiple_glycans: bool,
    compound_peptidoform: &CompoundPeptidoform,
) -> String {
    let ch = format!("{:+}", annotation.charge.value);
    format!(
        "<span>{}<sup class='charge'>{}</sup><sub style='--charge-width:{};'><span class='series'>{}</span><span class='glycan-id'>{}</span><span class='peptide-id'>{}</span></sub><span class='neutral-losses'>{}</span><span class='cross-links'>{}</span><span class='ambiguous-amino-acids'>{}</span><span class='modifications'>{}</span><span class='charge-carriers'>{}</span></span>",
        if let FragmentType::Oxonium(breakages) = &annotation.ion {
            breakages
            .iter()
            .filter(|b| !matches!(b, GlycanBreakPos::End(_)))
            .map(|b| format!(
                "{}<sub>{}</sub>",
                b.label(),
                b.position().label()
            ))
            .join("")
        } else {
            annotation.ion.label().to_string()
        },
        ch,
        ch.len(),
        annotation.ion.position_label().unwrap_or_default(),
        if multiple_glycans {
            if let FragmentType::Oxonium(breakages) = &annotation.ion {
                breakages[0].position().attachment()
            } else {
                annotation.ion.glycan_position().map(|g| g.attachment()).unwrap_or_default()
            }
        } else {
            String::new()
        },
        if multiple_peptidoforms && multiple_peptides {
            format!("p{}.{}", annotation.peptidoform_index+1, annotation.peptide_index + 1)
        }else if multiple_peptidoforms {
            format!("p{}", annotation.peptidoform_index+1)
        } else if multiple_peptides {
            format!("p{}", annotation.peptide_index + 1)
        } else {
            String::new()
        },
        annotation.neutral_loss.as_ref().map(|n| n.hill_notation_html()).unwrap_or_default(),
        get_xl(annotation),
        get_ambiguous_amino_acids(annotation, multiple_peptides),
        get_modifications(annotation, multiple_peptides, compound_peptidoform),
        get_charge_carriers(annotation),
    )
}

fn get_xl(annotation: &Fragment) -> String {
    let bound = annotation
        .formula
        .labels()
        .iter()
        .filter_map(|l| match l {
            AmbiguousLabel::CrossLinkBound(name) => Some(name),
            _ => None,
        })
        .unique()
        .sorted()
        .collect_vec();
    let broken = annotation
        .formula
        .labels()
        .iter()
        .filter_map(|l| match l {
            AmbiguousLabel::CrossLinkBroken(name, _) => Some(name),
            _ => None,
        })
        .unique()
        .sorted()
        .collect_vec();
    let mut output = String::new();
    if !bound.is_empty() {
        write!(
            &mut output,
            "<span class='cs-xl-intact'></span><sub>{}</sub>",
            bound.iter().join(",")
        )
        .unwrap();
    }
    if !broken.is_empty() {
        write!(
            &mut output,
            "<span class='cs-xl-broken'></span><sub>{}</sub>",
            broken.iter().join(",")
        )
        .unwrap();
    }
    output
}

fn get_ambiguous_amino_acids(annotation: &Fragment, multiple_peptides: bool) -> String {
    annotation
        .formula
        .labels()
        .iter()
        .flat_map(|label| {
            if let AmbiguousLabel::AminoAcid {
                option,
                sequence_index,
                peptide_index,
            } = label
            {
                Some(format!(
                    "{option}<sub>{}</sub>{}",
                    sequence_index + 1,
                    if multiple_peptides {
                        format!("<sub class='peptide-id'>p{}</sub>", peptide_index + 1)
                    } else {
                        String::new()
                    }
                ))
            } else {
                None
            }
        })
        .join("")
}

fn get_modifications(
    annotation: &Fragment,
    multiple_peptides: bool,
    compound_peptidoform: &CompoundPeptidoform,
) -> String {
    annotation
        .formula
        .labels()
        .iter()
        .flat_map(|label| {
            if let AmbiguousLabel::Modification {
                id,
                sequence_index,
                peptide_index,
            } = label
            {
                Some(format!(
                    "{}<sub>{}</sub>{}",
                    compound_peptidoform.peptidoforms()[annotation.peptidoform_index].peptides()
                        [annotation.peptide_index][*sequence_index]
                        .possible_modifications
                        .iter()
                        .find(|m| m.id == *id)
                        .unwrap()
                        .group,
                    display_sequence_index(*sequence_index),
                    if multiple_peptides {
                        format!("<sub class='peptide-id'>p{}</sub>", peptide_index + 1)
                    } else {
                        String::new()
                    }
                ))
            } else {
                None
            }
        })
        .join("")
}

fn get_charge_carriers(annotation: &Fragment) -> String {
    let str = annotation
        .formula
        .labels()
        .iter()
        .flat_map(|label| {
            if let AmbiguousLabel::ChargeCarrier(charge) = label {
                Some(charge.hill_notation_html())
            } else {
                None
            }
        })
        .join(",");
    if str.is_empty() {
        str
    } else {
        format!("[{str}]")
    }
}

pub fn display_sequence_index(sequence_index: SequencePosition) -> String {
    match sequence_index {
        SequencePosition::NTerm => "N-terminal".to_string(),
        SequencePosition::Index(i) => (i + 1).to_string(),
        SequencePosition::CTerm => "C-terminal".to_string(),
    }
}
