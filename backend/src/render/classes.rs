use itertools::Itertools;
use rustyms::{fragment::FragmentKind, Fragment};

/// Get all applicable classes for a set of annotations.
/// These are:
///   * The ion type(s) (eg 'precursor', 'y') (and 'multi' if there are multiple)
///   * The peptidoform(s) (eg 'p0', 'p1') (and 'mp' if there are multiple peptides)
///   * The peptide(s) (eg 'p0-2', 'p2-1') (and 'mpp' if there are multiple peptides)
///   * The position(s) (eg 'p0-2-3', 'p2-1-123')
///   * The unique peptide index (eg 'pu1', 'pu12')
///   * Other auxiliary classes ('oxonium', 'neutral-loss', & 'diagnostic')
pub fn get_classes(annotations: &[Fragment], unique_peptide_lookup: &[(usize, usize)]) -> String {
    let mut output = Vec::new();
    let mut shared_ion = annotations.first().map(|a| a.ion.kind());
    let mut first_peptidoform_index = None;
    let mut first_peptide_index = None;
    for annotation in annotations {
        output.push(annotation.ion.label().to_string());
        if annotation.ion.label() == "z·" {
            output.push("z".to_string());
        }
        output.push(format!(
            "p{}",
            annotation
                .peptidoform_index
                .map_or("?".to_string(), |i| i.to_string())
        ));
        output.push(format!(
            "p{}-{}",
            annotation
                .peptidoform_index
                .map_or("?".to_string(), |i| i.to_string()),
            annotation
                .peptide_index
                .map_or("?".to_string(), |i| i.to_string())
        ));
        if let Some(num) = first_peptidoform_index {
            if num != annotation.peptidoform_index && !output.contains(&"mp".to_string()) {
                output.push("mp".to_string());
            }
        } else {
            first_peptidoform_index = Some(annotation.peptidoform_index);
        }
        if let (Some(first_peptidoform_index), Some(fist_peptide_index)) =
            (first_peptidoform_index, first_peptide_index)
        {
            if first_peptidoform_index != annotation.peptidoform_index
                || fist_peptide_index != annotation.peptide_index
            {
                if !output.contains(&"mpp".to_string()) {
                    output.push("mpp".to_string())
                };
                let pu = format!(
                    "pu{}",
                    unique_peptide_lookup
                        .iter()
                        .position(|id| (Some(id.0), Some(id.1))
                            == (annotation.peptidoform_index, annotation.peptide_index))
                        .unwrap()
                );
                if !output.contains(&pu) {
                    output.push(pu)
                };
            }
        } else {
            first_peptide_index = Some(annotation.peptide_index);
            output.push(format!(
                "pu{}",
                unique_peptide_lookup
                    .iter()
                    .position(|id| (Some(id.0), Some(id.1))
                        == (annotation.peptidoform_index, annotation.peptide_index))
                    .unwrap()
            ))
        }
        if let Some(pos) = annotation.ion.position() {
            output.push(format!(
                "p{}-{}-{}",
                annotation
                    .peptidoform_index
                    .map_or("?".to_string(), |i| i.to_string()),
                annotation
                    .peptide_index
                    .map_or("?".to_string(), |i| i.to_string()),
                pos.sequence_index
            ));
        }
        if annotation.ion.kind() == FragmentKind::Oxonium {
            output.push("oxonium".to_string());
        }
        if annotation.ion.kind() == FragmentKind::diagnostic {
            output.push("diagnostic".to_string());
        }
        if !annotation.neutral_loss.is_empty() {
            output.push("neutral-loss".to_string());
        }
        if let Some(ion) = &shared_ion {
            if *ion != annotation.ion.kind() {
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
