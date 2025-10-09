use itertools::Itertools;
use mzannotate::{fragment::FragmentType, prelude::*};
use mzcore::{
    chemistry::AmbiguousLabel,
    prelude::*,
    sequence::{GnoComposition, Modification, SimpleModificationInner},
};
use std::fmt::Write;

use crate::Theme;

use super::render_glycan_fragment;

pub fn get_label(
    compound_peptidoform: &CompoundPeptidoformIon,
    annotations: &[Fragment],
    multiple_peptidoform_ions: bool,
    multiple_peptidoforms: bool,
    multiple_glycans: bool,
    glycan_footnotes: &mut Vec<String>,
    theme: Theme,
) -> String {
    if annotations.is_empty() {
        String::new()
    } else {
        let mut shared_charge = Some(annotations[0].charge);
        let mut shared_ion = Some(annotations[0].ion.label());
        let mut shared_pos = Some(annotations[0].ion.position_label());
        let mut shared_peptidoform_ion = Some(annotations[0].peptidoform_ion_index);
        let mut shared_peptidoform = Some(annotations[0].peptidoform_index);
        let mut shared_glycan = Some(annotations[0].ion.glycan_position().map(|g| g.attachment()));
        let first_glycan_image = get_glycan_figure(
            compound_peptidoform,
            &annotations[0],
            theme,
            glycan_footnotes,
        );
        let mut shared_glycan_figures = Some(first_glycan_image.clone());
        let mut shared_loss = Some(annotations[0].neutral_loss.clone());
        let mut shared_xl = Some(get_xl(&annotations[0]));
        let mut shared_ambiguous_amino_acids = Some(get_ambiguous_amino_acids(
            &annotations[0],
            multiple_peptidoforms,
            multiple_peptidoform_ions,
        ));
        let mut shared_modifications = Some(get_modifications(
            &annotations[0],
            multiple_peptidoforms,
            multiple_peptidoform_ions,
            compound_peptidoform,
        ));
        let mut shared_charge_carriers = Some(get_charge_carriers(&annotations[0]));
        let mut shared_glycan_peptide_fragments =
            Some(get_glycan_peptide_fragments(&annotations[0]));

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
            if let Some(peptidoform) = shared_peptidoform_ion {
                if peptidoform != a.peptidoform_ion_index {
                    shared_peptidoform_ion = None;
                }
            }
            if let Some(peptide) = shared_peptidoform {
                if peptide != a.peptidoform_index {
                    shared_peptidoform = None;
                }
            }
            if let Some(glycan) = &shared_glycan {
                if *glycan != a.ion.glycan_position().map(|g| g.attachment()) {
                    shared_glycan = None;
                }
            }
            if let Some(glycan) = &shared_glycan_figures {
                if *glycan != get_glycan_figure(compound_peptidoform, a, theme, glycan_footnotes) {
                    shared_glycan_figures = None;
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
                if *aaa
                    != get_ambiguous_amino_acids(
                        a,
                        multiple_peptidoforms,
                        multiple_peptidoform_ions,
                    )
                {
                    shared_ambiguous_amino_acids = None;
                }
            }
            if let Some(sm) = &shared_modifications {
                if *sm
                    != get_modifications(
                        a,
                        multiple_peptidoforms,
                        multiple_peptidoform_ions,
                        compound_peptidoform,
                    )
                {
                    shared_modifications = None;
                }
            }
            if let Some(cc) = &shared_charge_carriers {
                if *cc != get_charge_carriers(a) {
                    shared_charge_carriers = None;
                }
            }
            if let Some(cc) = &shared_glycan_peptide_fragments {
                if *cc != get_glycan_peptide_fragments(a) {
                    shared_glycan_peptide_fragments = None;
                }
            }
        }

        if shared_charge.is_none()
            && shared_ion.is_none()
            && shared_pos.is_none()
            && shared_peptidoform_ion.is_none()
            && shared_peptidoform.is_none()
            && shared_glycan.is_none()
            && shared_loss.is_none()
            && shared_xl.is_none()
            && shared_ambiguous_amino_acids.is_none()
            && shared_modifications.is_none()
            && shared_charge_carriers.is_none()
            && shared_glycan_peptide_fragments.is_none()
        {
            "*".to_string()
        } else {
            let charge_str = shared_charge
                .map(|charge| format!("{:+}", charge.value))
                .unwrap_or("*".to_string());
            let ion_str = shared_ion
                .map(|(sup, c)| {
                    if let Some(sup) = sup {
                        format!("<sup>{sup}</sup>{c}")
                    } else {
                        c.to_string()
                    }
                })
                .unwrap_or("*".to_string());
            let pos_str = shared_pos
                .map(|pos| pos.unwrap_or_default())
                .unwrap_or("*".to_string());
            let peptidoform_str = shared_peptidoform_ion
                .flatten()
                .map(|pep| (pep + 1).to_string())
                .unwrap_or("*".to_string());
            let peptide_str = shared_peptidoform
                .flatten()
                .map(|pep| (pep + 1).to_string())
                .unwrap_or("*".to_string());
            let glycan_str = shared_glycan
                .unwrap_or(Some("*".to_string()))
                .unwrap_or_default();
            let glycan_figure_str = shared_glycan_figures
                .unwrap_or(first_glycan_image)
                .unwrap_or_default();
            let loss_str = shared_loss
                .map(|o| o.iter().map(|n| n.hill_notation_html()).join(""))
                .unwrap_or("*".to_string());
            let xl_str = shared_xl.unwrap_or("*".to_string());
            let aaa_str = shared_ambiguous_amino_acids.unwrap_or("*".to_string());
            let sm_str = shared_modifications.unwrap_or("*".to_string());
            let cc_str = shared_charge_carriers.unwrap_or("*".to_string());
            let glycan_peptide_fragments_str =
                shared_glycan_peptide_fragments.unwrap_or("*".to_string());

            let multi = if annotations.len() > 1 {
                let mut multi = String::new();
                for annotation in annotations {
                    write!(
                        multi,
                        "{}",
                        get_single_label(
                            annotation,
                            multiple_peptidoform_ions,
                            multiple_peptidoforms,
                            multiple_glycans,
                            compound_peptidoform,
                            glycan_footnotes,
                            theme,
                        )
                    )
                    .unwrap();
                }
                format!("<span class='multi'>{multi}</span>")
            } else {
                String::new()
            };
            let single_internal_glycan =
                matches!(annotations[0].ion, FragmentType::B { .. }) && annotations.len() == 1;

            if single_internal_glycan {
                get_single_label(
                    &annotations[0],
                    multiple_peptidoform_ions,
                    multiple_peptidoforms,
                    multiple_glycans,
                    compound_peptidoform,
                    glycan_footnotes,
                    theme,
                )
            } else {
                format!(
                    "{}<span>{}<sup class='charge'>{}</sup><sub style='--charge-width:{};'><span class='series'>{}</span><span class='glycan-id'>{}</span><span class='peptide-id'>{}</span></sub><span class='neutral-losses'>{}</span><span class='cross-links'>{}</span><span class='ambiguous-amino-acids'>{}</span><span class='modifications'>{}</span><span class='charge-carriers'>{}</span><span class='glycan-peptide-fragment'>{}</span></span>{}",
                    glycan_figure_str,
                    ion_str,
                    charge_str,
                    charge_str.len(),
                    pos_str,
                    if multiple_glycans {
                        glycan_str
                    } else {
                        String::new()
                    },
                    if multiple_peptidoform_ions && multiple_peptidoforms {
                        format!("p{peptidoform_str}.{peptide_str}")
                    } else if multiple_peptidoform_ions {
                        format!("p{peptidoform_str}")
                    } else if multiple_peptidoforms {
                        format!("p{peptide_str}")
                    } else {
                        String::new()
                    },
                    loss_str,
                    xl_str,
                    aaa_str,
                    sm_str,
                    cc_str,
                    glycan_peptide_fragments_str,
                    multi,
                )
            }
        }
    }
}

fn get_glycan_figure(
    compound_peptidoform: &CompoundPeptidoformIon,
    annotation: &Fragment,
    theme: Theme,
    glycan_footnotes: &mut Vec<String>,
) -> Option<String> {
    let multiple_peptidoform_ions = compound_peptidoform.peptidoform_ions().len() > 1;
    let multiple_peptidoforms = compound_peptidoform
        .peptidoform_ions()
        .iter()
        .any(|p| p.peptidoforms().len() > 1);
    annotation
        .ion
        .glycan_break_positions()
        .and_then(|g| {
            annotation.peptidoform_ion_index.and_then(|pii| {
                annotation.peptidoform_index.and_then(|pi| {
                    g.0.and_then(|index| {
                        match index {
                            SequencePosition::NTerm => compound_peptidoform.peptidoform_ions()[pii]
                                .peptidoforms()[pi]
                                .get_n_term(),
                            SequencePosition::Index(i) => {
                                compound_peptidoform.peptidoform_ions()[pii].peptidoforms()[pi]
                                    .sequence()[i]
                                    .modifications
                                    .as_slice()
                            }
                            SequencePosition::CTerm => compound_peptidoform.peptidoform_ions()[pii]
                                .peptidoforms()[pi]
                                .get_c_term(),
                        }
                        .iter()
                        .find_map(|m| match m.simple().map(|m| m.as_ref()) {
                            Some(SimpleModificationInner::Gno {
                                composition: GnoComposition::Topology(structure),
                                ..
                            }) => Some(structure),
                            Some(SimpleModificationInner::GlycanStructure(structure)) => {
                                Some(structure)
                            }
                            _ => None,
                        })
                        .map(|s| (s, g.1, pii, pi))
                    })
                })
            })
        })
        .map(|(structure, selection, pii, pi)| {
            render_glycan_fragment(
                structure,
                selection,
                theme,
                glycan_footnotes,
                multiple_peptidoform_ions,
                multiple_peptidoforms,
                pii,
                pi,
            )
        })
}

fn get_single_label(
    annotation: &Fragment,
    multiple_peptidoform_ions: bool,
    multiple_peptidoforms: bool,
    multiple_glycans: bool,
    compound_peptidoform: &CompoundPeptidoformIon,
    glycan_footnotes: &mut Vec<String>,
    theme: Theme,
) -> String {
    let ch = format!("{:+}", annotation.charge.value);
    format!(
        "{}<span>{}<sup class='charge'>{}</sup><sub style='--charge-width:{};'><span class='series'>{}</span><span class='glycan-id'>{}</span><span class='peptide-id'>{}</span></sub><span class='neutral-losses'>{}</span><span class='cross-links'>{}</span><span class='ambiguous-amino-acids'>{}</span><span class='modifications'>{}</span><span class='charge-carriers'>{}</span><span class='glycan-peptide-fragment'>{}</span></span>",
        get_glycan_figure(compound_peptidoform, annotation, theme, glycan_footnotes)
            .unwrap_or_default(),
        if let FragmentType::B { b, y, .. } = &annotation.ion {
            format!("B<sub>{}</sub>", b.label())
                + &y.iter()
                    .map(|b| format!("Y<sub>{}</sub>", b.label()))
                    .join("")
        } else {
            let (sup, label) = annotation.ion.label();
            if let Some(sup) = sup {
                format!("<sup>{sup}</sup>{label}")
            } else {
                label.to_string()
            }
        },
        ch,
        ch.len(),
        if let FragmentType::B { .. } = &annotation.ion {
            String::new()
        } else {
            annotation.ion.position_label().unwrap_or_default()
        },
        if multiple_glycans {
            if let FragmentType::B { b, .. } = &annotation.ion {
                b.attachment()
            } else {
                annotation
                    .ion
                    .glycan_position()
                    .map(|g| g.attachment())
                    .unwrap_or_default()
            }
        } else {
            String::new()
        },
        if multiple_peptidoform_ions && multiple_peptidoforms {
            format!(
                "p{}.{}",
                annotation
                    .peptidoform_ion_index
                    .map_or("?".to_string(), |i| (i + 1).to_string()),
                annotation
                    .peptidoform_index
                    .map_or("?".to_string(), |i| (i + 1).to_string())
            )
        } else if multiple_peptidoform_ions {
            format!(
                "p{}",
                annotation
                    .peptidoform_ion_index
                    .map_or("?".to_string(), |i| (i + 1).to_string())
            )
        } else if multiple_peptidoforms {
            format!(
                "p{}",
                annotation
                    .peptidoform_index
                    .map_or("?".to_string(), |i| (i + 1).to_string())
            )
        } else {
            String::new()
        },
        annotation
            .neutral_loss
            .iter()
            .map(|n| n.hill_notation_html())
            .join(""),
        get_xl(annotation),
        get_ambiguous_amino_acids(annotation, multiple_peptidoforms, multiple_peptidoform_ions),
        get_modifications(
            annotation,
            multiple_peptidoforms,
            multiple_peptidoform_ions,
            compound_peptidoform,
        ),
        get_charge_carriers(annotation),
        get_glycan_peptide_fragments(annotation),
    )
}

fn get_xl(annotation: &Fragment) -> String {
    let bound = annotation
        .formula
        .iter()
        .flat_map(|f| f.labels())
        .filter_map(|l| match l {
            AmbiguousLabel::CrossLinkBound(name) => Some(name),
            _ => None,
        })
        .unique()
        .sorted()
        .collect_vec();
    let broken = annotation
        .formula
        .iter()
        .flat_map(|f| f.labels())
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

fn get_ambiguous_amino_acids(
    annotation: &Fragment,
    multiple_peptidoforms: bool,
    multiple_peptidoform_ions: bool,
) -> String {
    annotation
        .formula
        .iter()
        .flat_map(|f| f.labels())
        .filter_map(|label| {
            if let AmbiguousLabel::AminoAcid {
                option,
                sequence_index,
                peptidoform_index,
                peptidoform_ion_index,
            } = label
            {
                Some(format!(
                    "{option}<sub>{}</sub>{}",
                    sequence_index + 1,
                    if multiple_peptidoforms && multiple_peptidoform_ions {
                        format!(
                            "<sub class='peptide-id'>p{}.{}</sub>",
                            peptidoform_ion_index + 1,
                            peptidoform_index + 1
                        )
                    } else if multiple_peptidoforms {
                        format!("<sub class='peptide-id'>p{}</sub>", peptidoform_index + 1)
                    } else if multiple_peptidoform_ions {
                        format!(
                            "<sub class='peptide-id'>p{}</sub>",
                            peptidoform_ion_index + 1
                        )
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
    multiple_peptidoforms: bool,
    multiple_peptidoform_ions: bool,
    compound_peptidoform: &CompoundPeptidoformIon,
) -> String {
    annotation
        .formula
        .iter()
        .flat_map(|f| f.labels())
        .filter_map(|label| {
            if let AmbiguousLabel::Modification {
                id,
                sequence_index,
                peptidoform_index,
                peptidoform_ion_index,
            } = label
            {
                Some(format!(
                    "{}<sub>{}</sub>{}",
                    match sequence_index {
                        SequencePosition::NTerm => compound_peptidoform.peptidoform_ions()
                            [annotation.peptidoform_ion_index?]
                            .peptidoforms()[annotation.peptidoform_index?]
                            .get_n_term(),
                        SequencePosition::Index(i) => compound_peptidoform.peptidoform_ions()
                            [annotation.peptidoform_ion_index?]
                            .peptidoforms()[annotation.peptidoform_index?]
                            .sequence()[*i]
                            .modifications
                            .as_slice(),
                        SequencePosition::CTerm => compound_peptidoform.peptidoform_ions()
                            [annotation.peptidoform_index?]
                            .peptidoforms()[annotation.peptidoform_index?]
                            .get_c_term(),
                    }
                    .iter()
                    .find_map(
                        |m| if let Modification::Ambiguous { id: mid, group, .. } = m {
                            (*mid == *id).then_some(group)
                        } else {
                            None
                        }
                    )
                    .unwrap(), // Can crash
                    display_sequence_index(*sequence_index),
                    if multiple_peptidoforms && multiple_peptidoform_ions {
                        format!(
                            "<sub class='peptide-id'>p{}.{}</sub>",
                            peptidoform_ion_index + 1,
                            peptidoform_index + 1
                        )
                    } else if multiple_peptidoforms {
                        format!("<sub class='peptide-id'>p{}</sub>", peptidoform_index + 1)
                    } else if multiple_peptidoform_ions {
                        format!(
                            "<sub class='peptide-id'>p{}</sub>",
                            peptidoform_ion_index + 1
                        )
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
        .iter()
        .flat_map(|f| f.labels())
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

fn get_glycan_peptide_fragments(annotation: &Fragment) -> String {
    annotation
        .formula
        .iter()
        .flat_map(|f| f.labels())
        .flat_map(|label| {
            if matches!(
                label,
                AmbiguousLabel::GlycanFragment(_) | AmbiguousLabel::GlycanFragmentComposition(_)
            ) {
                Some(format!("+{label}"))
            } else {
                None
            }
        })
        .join("")
}

pub fn display_sequence_index(sequence_index: SequencePosition) -> String {
    match sequence_index {
        SequencePosition::NTerm => "N-terminal".to_string(),
        SequencePosition::Index(i) => (i + 1).to_string(),
        SequencePosition::CTerm => "C-terminal".to_string(),
    }
}
