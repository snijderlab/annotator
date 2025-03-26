use std::num::{IntErrorKind, ParseIntError};

use itertools::Itertools;
use mzsignal::text;
use rustyms::{
    error::{Context, CustomError},
    placement_rule::PlacementRule,
    AminoAcid, MolecularFormula, NeutralLoss,
};

use crate::render::{display_formula, display_neutral_loss, display_placement_rule, display_stubs};

#[tauri::command]
pub fn validate_molecular_formula(text: String) -> Result<String, CustomError> {
    text.parse::<f64>()
        .map(MolecularFormula::with_additional_mass)
        .or_else(|_| MolecularFormula::from_pro_forma(&text, .., true, true, true))
        .map(|f| display_formula(&f, true))
}

pub fn parse_amino_acid(text: &str) -> Result<AminoAcid, CustomError> {
    text.parse::<AminoAcid>().map_err(|()| {
        CustomError::error(
            "Invalid amino acid",
            "This is not a valid amino acid code, use `A` not `Ala` or `Alanine`",
            Context::Show {
                line: text.to_string(),
            },
        )
    })
}

#[tauri::command]
pub fn validate_amino_acid(text: String) -> Result<String, CustomError> {
    parse_amino_acid(&text).map(|f| f.to_string())
}

#[tauri::command]
pub fn validate_neutral_loss(text: String) -> Result<String, CustomError> {
    text.parse::<NeutralLoss>()
        .map(|f| display_neutral_loss(&f))
}

pub fn parse_aa_neutral_loss(
    text: &str,
) -> Result<(Vec<AminoAcid>, Vec<NeutralLoss>), CustomError> {
    if let Some(index) = text.find(':') {
        let aas = parse_aa_list(&text[..index])?;

        let mut pos = index + 1;
        let mut losses = Vec::new();
        for loss in text[index + 1..].split(',') {
            losses.push(loss.parse::<NeutralLoss>().map_err(|err| {
                err.with_context(Context::Line {
                    line_index: None,
                    line: text.to_string(),
                    offset: pos,
                    length: loss.len(),
                })
            })?);
            pos += loss.len() + 1;
        }

        Ok((aas, losses))
    } else {
        Err(CustomError::error(
            "Invalid amino acid neutral loss",
            "An amino acid neutral loss should be specified with 'A,C:-O1H2,+H2'",
            Context::full_line(0, text),
        ))
    }
}

pub fn display_aa_neutral_loss(aa: &[AminoAcid], loss: &[NeutralLoss]) -> String {
    aa.iter().join(",") + ":" + &loss.iter().map(display_neutral_loss).join(",")
}

#[tauri::command]
pub fn validate_aa_neutral_loss(text: String) -> Result<String, CustomError> {
    parse_aa_neutral_loss(&text).map(|s| display_aa_neutral_loss(&s.0, &s.1))
}

pub fn parse_aa_list(text: &str) -> Result<Vec<AminoAcid>, CustomError> {
    let mut pos = 0;
    let mut aas = Vec::new();
    for aa in text.split(',') {
        aas.push(parse_amino_acid(aa).map_err(|err| {
            err.with_context(Context::Line {
                line_index: None,
                line: text.to_string(),
                offset: pos,
                length: aa.len(),
            })
        })?);
        pos += aa.len() + 1;
    }
    Ok(aas)
}

pub fn parse_satellite_ion(text: &str) -> Result<(Vec<AminoAcid>, u8), CustomError> {
    if let Some(index) = text.find(':') {
        let aas = parse_aa_list(&text[..index])?;

        let maximal_distance = text[index + 1..].parse::<u8>().map_err(|err| {
            CustomError::error(
                "Invalid maximal distance",
                format!("The maximal distance {}", explain_number_error(&err)),
                Context::Line {
                    line_index: None,
                    line: text.to_string(),
                    offset: index + 1,
                    length: text[index + 1..].len(),
                },
            )
        })?;

        Ok((aas, maximal_distance))
    } else {
        Err(CustomError::error(
            "Invalid satellite ion location",
            "A satellite ion location should be specified with 'A,C:2'",
            Context::full_line(0, text),
        ))
    }
}

pub fn display_satellite_ion(aa: &[AminoAcid], distance: u8) -> String {
    aa.iter().join(",") + ":" + &distance.to_string()
}

#[tauri::command]
pub fn validate_satellite_ion(text: String) -> Result<String, CustomError> {
    parse_satellite_ion(&text).map(|s| display_satellite_ion(&s.0, s.1))
}

#[tauri::command]
pub fn validate_placement_rule(text: String) -> Result<String, CustomError> {
    text.parse::<PlacementRule>()
        .map(|r| display_placement_rule(&r, false))
}

pub fn parse_stub(text: &str) -> Result<(MolecularFormula, MolecularFormula), CustomError> {
    if let Some(index) = text.find(':') {
        let f1 = text[..index]
            .parse::<f64>()
            .map(MolecularFormula::with_additional_mass)
            .or_else(|_| MolecularFormula::from_pro_forma(text, ..index, true, true, true))?;
        let f2 = text[index + 1..]
            .parse::<f64>()
            .map(MolecularFormula::with_additional_mass)
            .or_else(|_| MolecularFormula::from_pro_forma(text, index + 1.., true, true, true))?;
        Ok((f1, f2))
    } else {
        Err(CustomError::error(
            "Invalid breakage",
            "A breakage should be specified with 'formula1:formula2'",
            Context::full_line(0, text),
        ))
    }
}

#[tauri::command]
pub fn validate_stub(text: String) -> Result<String, CustomError> {
    parse_stub(&text).map(|s| display_stubs(&s, true))
}

#[tauri::command]
pub fn validate_custom_single_specificity(
    placement_rules: Vec<String>,
    neutral_losses: Vec<String>,
    diagnostic_ions: Vec<String>,
) -> Result<String, CustomError> {
    let rules = placement_rules
        .into_iter()
        .map(|text| text.parse::<PlacementRule>())
        .collect::<Result<Vec<_>, _>>()?;
    let rules = if rules.is_empty() {
        vec![PlacementRule::Anywhere]
    } else {
        rules
    };

    let neutral_losses = neutral_losses
        .into_iter()
        .map(|text| text.parse::<NeutralLoss>())
        .collect::<Result<Vec<_>, _>>()?;

    let diagnostic_ions = diagnostic_ions
        .into_iter()
        .map(|text| {
            text.parse::<f64>()
                .map(MolecularFormula::with_additional_mass)
                .or_else(|_| MolecularFormula::from_pro_forma(&text, .., true, true, true))
        })
        .collect::<Result<Vec<_>, _>>()?;
    Ok(format!(
        "<span data-value='{{\"placement_rules\":[{}],\"neutral_losses\":[{}],\"diagnostic_ions\":[{}]}}'>Placement rules: {}{}{}</span>",
        rules.iter().map(|r| format!("\"{}\"", display_placement_rule(r, false))).join(","),
        neutral_losses.iter().map(|n| format!("\"{}\"", n.hill_notation())).join(","),
        diagnostic_ions.iter().map(|n| format!("\"{}\"", n.hill_notation())).join(","),
        rules.iter().map(|p|display_placement_rule(p,true)).join(", "),
        if neutral_losses.is_empty() {
            String::new()
        } else {
            ", Neutral losses: ".to_string() + &neutral_losses.iter().map(display_neutral_loss).join(", ")
        },
        if diagnostic_ions.is_empty() {
            String::new()
        } else {
            ", Diagnostic ions: ".to_string() + &diagnostic_ions.iter().map(|f| display_formula(f, true)).join(", ")
        },
    ))
}

#[tauri::command]
pub fn validate_custom_linker_specificity(
    asymmetric: bool,
    placement_rules: Vec<String>,
    secondary_placement_rules: Vec<String>,
    stubs: Vec<String>,
    diagnostic_ions: Vec<String>,
) -> Result<String, CustomError> {
    let rules1 = placement_rules
        .into_iter()
        .map(|text| text.parse::<PlacementRule>())
        .collect::<Result<Vec<_>, _>>()?;
    let rules1 = if rules1.is_empty() {
        vec![PlacementRule::Anywhere]
    } else {
        rules1
    };
    let rules2 = secondary_placement_rules
        .into_iter()
        .map(|text| text.parse::<PlacementRule>())
        .collect::<Result<Vec<_>, _>>()?;
    let rules2 = if rules2.is_empty() {
        vec![PlacementRule::Anywhere]
    } else {
        rules2
    };

    let stubs = stubs
        .into_iter()
        .map(|text| parse_stub(&text))
        .collect::<Result<Vec<_>, _>>()?;

    let diagnostic_ions = diagnostic_ions
        .into_iter()
        .map(|text| {
            text.parse::<f64>()
                .map(MolecularFormula::with_additional_mass)
                .or_else(|_| MolecularFormula::from_pro_forma(&text, .., true, true, true))
        })
        .collect::<Result<Vec<_>, _>>()?;
    Ok(format!(
        "<span data-value='{{\"asymmetric\":{asymmetric},\"placement_rules\":[{}],\"secondary_placement_rules\":[{}],\"stubs\":[{}],\"diagnostic_ions\":[{}]}}'>Placement rules: {}{}{}{}</span>",
        rules1.iter().map(|r| format!("\"{}\"", display_placement_rule(r, false))).join(","),
        rules2.iter().map(|r| format!("\"{}\"", display_placement_rule(r, false))).join(","),
        stubs.iter().map(|s| format!("\"{}\"", display_stubs(s, false))).join(","),
        diagnostic_ions.iter().map(|n| format!("\"{}\"", n.hill_notation())).join(","),
        rules1.iter().map(|p|display_placement_rule(p,true)).join(", "),
        if asymmetric {
            ", Secondary placement rules: ".to_string() + &rules2.iter().map(|p|display_placement_rule(p,true)).join(", ")
        } else {
            String::new()
        },
        if stubs.is_empty() {
            String::new()
        } else {
            ", Breakage: ".to_string() + &stubs.iter().map(|s| display_stubs(s, true)).join(", ")
        },
        if diagnostic_ions.is_empty() {
            String::new()
        } else {
            ", Diagnostic ions: ".to_string() + &diagnostic_ions.iter().map(|f| display_formula(f, true)).join(", ")
        },
    ))
}

/// To be used as `The xx number ` + the explanation from here (does not have a dot).
pub const fn explain_number_error(error: &ParseIntError) -> &'static str {
    match error.kind() {
        IntErrorKind::Empty => "is empty",
        IntErrorKind::InvalidDigit => "contains an invalid character",
        IntErrorKind::NegOverflow => "is too small to fit in the internal representation",
        IntErrorKind::PosOverflow => "is too big to fit in the internal representation",
        IntErrorKind::Zero => "is zero, which is not allowed here",
        _ => "is not a valid number",
    }
}
