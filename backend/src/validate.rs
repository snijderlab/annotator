use std::num::{IntErrorKind, ParseIntError};

    use custom_error::{BoxedError, Context, CustomErrorTrait};
use itertools::Itertools;
use rustyms::{ fragment::{FragmentKind, NeutralLoss}, glycan::MonoSaccharide, prelude::*, sequence::PlacementRule};

use crate::render::{display_formula, display_neutral_loss, display_placement_rule, display_stubs};

#[tauri::command]
pub fn validate_molecular_formula(text: String) -> Result<String, BoxedError<'static>> {
    text.parse::<f64>()
        .map(MolecularFormula::with_additional_mass)
        .or_else(|_| MolecularFormula::from_pro_forma(&text, .., true, true, true, false))
        .map(|f| display_formula(&f, true))
        .map_err(BoxedError::to_owned)
}

pub fn parse_amino_acid(text: &str) -> Result<AminoAcid, BoxedError<'_>> {
    text.parse::<AminoAcid>().map_err(|()| {
        BoxedError::error(
            "Invalid amino acid",
            "This is not a valid amino acid code, use `A` not `Ala` or `Alanine`",
            Context::show(text.to_string()),
        )
    })
}

#[tauri::command]
pub fn validate_amino_acid(text: String) -> Result<String, BoxedError<'static>> {
    parse_amino_acid(&text).map(|f| f.to_string()).map_err(BoxedError::to_owned)
}

pub fn parse_fragment_kind(text: &str) -> Result<FragmentKind, BoxedError> {
    match text {
        "a" => Ok(FragmentKind::a),
        "b" => Ok(FragmentKind::b),
        "c" => Ok(FragmentKind::c),
        "d" => Ok(FragmentKind::d),
        "v" => Ok(FragmentKind::v),
        "w" => Ok(FragmentKind::w),
        "x" => Ok(FragmentKind::x),
        "y" => Ok(FragmentKind::y),
        "z" => Ok(FragmentKind::z),
        "immonium" => Ok(FragmentKind::immonium),
        _ => Err(BoxedError::error("Invalid fragment kind", "Use any of a/b/c/d/v/w/x/y/z/immonium, notice that this is case sensitive.", Context::none()))
    }
}

#[tauri::command]
pub fn validate_fragment_kind(text: String) -> Result<String, BoxedError<'static>> {
    parse_fragment_kind(&text).map(|f| f.to_string()).map_err(BoxedError::to_owned)
}

#[tauri::command]
pub fn validate_neutral_loss(text: String) -> Result<String, BoxedError<'static>> {
    text.parse::<NeutralLoss>()
        .map(|f| display_neutral_loss(&f, false))
}

pub fn parse_aa_neutral_loss(
    text: &str,
) -> Result<(Vec<AminoAcid>, Vec<NeutralLoss>), BoxedError<'_>> {
    if let Some(index) = text.find(':') {
        let aas = parse_aa_list(&text[..index])?;

        let mut pos = index + 1;
        let mut losses = Vec::new();
        for loss in text[index + 1..].split(',') {
            losses.push(loss.parse::<NeutralLoss>().map_err(|err| {
                err.add_context(Context::line(
                    None,
                    text.to_string(),
                    pos,
                    loss.len(),
                ))
            })?);
            pos += loss.len() + 1;
        }

        Ok((aas, losses))
    } else {
        Err(BoxedError::error(
            "Invalid amino acid neutral loss",
            "An amino acid neutral loss should be specified with 'A,C:-O1H2,+H2'",
            Context::full_line(0, text),
        ))
    }
}

pub fn display_aa_neutral_loss(aa: &[AminoAcid], loss: &[NeutralLoss]) -> String {
    aa.iter().join(",") + ":" + &loss.iter().map(|l| l.hill_notation()).join(",")
}

#[tauri::command]
pub fn validate_aa_neutral_loss(text: String) -> Result<String, BoxedError<'static>> {
    parse_aa_neutral_loss(&text).map(|s| display_aa_neutral_loss(&s.0, &s.1)).map_err(BoxedError::to_owned)
}

pub fn parse_monosaccharide_neutral_loss(text: &str) -> Result<(MonoSaccharide, bool, Vec<NeutralLoss>), BoxedError<'_>> {
    let (monosaccharide, losses, specific) = text.split_once(':').map(|(a, b)| (a, b, false)).or(text.split_once('=').map(|(a, b)| (a, b, true))).ok_or_else(|| BoxedError::error("Invalid monosaccharide neutral loss", "The monosaccharide neutral loss has the be in the form of 'sugar:loss,loss' or 'sugar=loss,loss' bot neither ':' not '=' are found in the text.", Context::none()))?;

    let mut pos = monosaccharide.len() + 1;
    let mut found_losses = Vec::new();
    for loss in losses.split(',') {
        found_losses.push(loss.parse::<NeutralLoss>().map_err(|err| {
            err.add_context(Context::line(
                None,
                text.to_string(),
                pos,
                loss.len(),
            ))
        })?);
        pos += loss.len() + 1;
    }

    Ok((MonoSaccharide::from_short_iupac(monosaccharide, 0, 0).and_then(|(sugar, amount_parsed)| 
        if amount_parsed == monosaccharide.len() {
            Ok(sugar.with_name(monosaccharide))
        } else {
            Err(BoxedError::error("Could not parse monosaccharide", format!("The monosaccharide was interpreted to mean '{sugar}' but this left text unaccounted for"), Context::line_range(None, text, amount_parsed..monosaccharide.len())))
        }
    )?, specific, found_losses))
}

pub fn display_monosaccharide_neutral_loss(rule: (MonoSaccharide, bool, Vec<NeutralLoss>)) -> String {
    format!("{}{}{}", rule.0, if rule.1 {"="} else {":"}, rule.2.iter().join(","))
}

#[tauri::command]
pub fn validate_monosaccharide_neutral_loss(text: String) -> Result<String, BoxedError<'static>> {
    parse_monosaccharide_neutral_loss(&text).map(display_monosaccharide_neutral_loss).map_err(BoxedError::to_owned)
}

pub fn parse_aa_list(text: &str) -> Result<Vec<AminoAcid>, BoxedError<'_>> {
    let mut pos = 0;
    let mut aas = Vec::new();
    for aa in text.split(',') {
        aas.push(parse_amino_acid(aa).map_err(|err| {
            err.add_context(Context::line(
                None,
                text.to_string(),
                pos,
                aa.len(),
            ))
        })?);
        pos += aa.len() + 1;
    }
    Ok(aas)
}

pub fn parse_satellite_ion(text: &str) -> Result<(Vec<AminoAcid>, u8), BoxedError<'_>> {
    if let Some(index) = text.find(':') {
        let aas = parse_aa_list(&text[..index])?;

        let maximal_distance = text[index + 1..].parse::<u8>().map_err(|err| {
            BoxedError::error(
                "Invalid maximal distance",
                format!("The maximal distance {}", explain_number_error(&err)),
                Context::line(
                    None,
                    text,
                    index + 1,
                    text[index + 1..].len(),
                ),
            )
        })?;

        Ok((aas, maximal_distance))
    } else {
        Err(BoxedError::error(
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
pub fn validate_satellite_ion(text: String) -> Result<String, BoxedError<'static>> {
    parse_satellite_ion(&text).map(|s| display_satellite_ion(&s.0, s.1)).map_err(BoxedError::to_owned)
}

#[tauri::command]
pub fn validate_placement_rule(text: String) -> Result<String, BoxedError<'static>> {
    text.parse::<PlacementRule>()
        .map(|r| display_placement_rule(&r, false))
}

pub fn parse_stub(text: &str) -> Result<(MolecularFormula, MolecularFormula), BoxedError<'static>> {
    if let Some(index) = text.find(':') {
        let f1 = text[..index]
            .parse::<f64>()
            .map(MolecularFormula::with_additional_mass)
            .or_else(|_| MolecularFormula::from_pro_forma(text, ..index, true, true, true, false))
            .map_err(BoxedError::to_owned)?;
        let f2 = text[index + 1..]
            .parse::<f64>()
            .map(MolecularFormula::with_additional_mass)
            .or_else(|_| MolecularFormula::from_pro_forma(text, index + 1.., true, true, true, false))
            .map_err(BoxedError::to_owned)?;
        Ok((f1, f2))
    } else {
        Err(BoxedError::error(
            "Invalid breakage",
            "A breakage should be specified with 'formula1:formula2'",
            Context::full_line(0, text.to_string()),
        ))
    }
}

#[tauri::command]
pub fn validate_stub(text: String) -> Result<String, BoxedError<'static>> {
    parse_stub(&text).map(|s| display_stubs(&s, true)).map_err(BoxedError::to_owned)
}

#[tauri::command]
pub fn validate_custom_single_specificity(
    placement_rules: Vec<String>,
    neutral_losses: Vec<String>,
    diagnostic_ions: Vec<String>,
) -> Result<String, BoxedError<'static>> {
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
                .or_else(|_| MolecularFormula::from_pro_forma(&text, .., true, true, true, false))
                .map_err(BoxedError::to_owned)
        })
        .collect::<Result<Vec<_>, _>>()?;
    Ok(format!(
        "<span data-value='{{\"placement_rules\":[{}],\"neutral_losses\":[{}],\"diagnostic_ions\":[{}]}}'>Placement rules: {}{}{}</span>",
        rules
            .iter()
            .map(|r| format!("\"{}\"", display_placement_rule(r, false)))
            .join(","),
        neutral_losses
            .iter()
            .map(|n| format!("\"{}\"", n.hill_notation()))
            .join(","),
        diagnostic_ions
            .iter()
            .map(|n| format!("\"{}\"", n.hill_notation()))
            .join(","),
        rules
            .iter()
            .map(|p| display_placement_rule(p, true))
            .join(", "),
        if neutral_losses.is_empty() {
            String::new()
        } else {
            ", Neutral losses: ".to_string()
                + &neutral_losses.iter().map(|n| display_neutral_loss(n, true)).join(", ")
        },
        if diagnostic_ions.is_empty() {
            String::new()
        } else {
            ", Diagnostic ions: ".to_string()
                + &diagnostic_ions
                    .iter()
                    .map(|f| display_formula(f, true))
                    .join(", ")
        },
    ))
}

#[tauri::command]
pub fn validate_custom_linker_specificity(
    asymmetric: bool,
    placement_rules: Vec<String>,
    secondary_placement_rules: Vec<String>,
    stubs: Vec<String>,
    neutral_losses: Vec<String>,
    diagnostic_ions: Vec<String>,
) -> Result<String, BoxedError<'static>> {
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
        .collect::<Result<Vec<_>, _>>()
        .map_err(BoxedError::to_owned)?;

    let neutral_losses = neutral_losses
        .into_iter()
        .map(|text| text.parse::<NeutralLoss>())
        .collect::<Result<Vec<_>, _>>()?;

    let diagnostic_ions = diagnostic_ions
        .into_iter()
        .map(|text| {
            text.parse::<f64>()
                .map(MolecularFormula::with_additional_mass)
                .or_else(|_| MolecularFormula::from_pro_forma(&text, .., true, true, true, false)).map_err(BoxedError::to_owned)
        })
        .collect::<Result<Vec<_>, _>>()?;
    Ok(format!(
        "<span data-value='{{\"asymmetric\":{asymmetric},\"placement_rules\":[{}],\"secondary_placement_rules\":[{}],\"stubs\":[{}],\"neutral_losses\":[{}],\"diagnostic_ions\":[{}]}}'>Placement rules: {}{}{}{}{}</span>",
        rules1
            .iter()
            .map(|r| format!("\"{}\"", display_placement_rule(r, false)))
            .join(","),
        rules2
            .iter()
            .map(|r| format!("\"{}\"", display_placement_rule(r, false)))
            .join(","),
        stubs
            .iter()
            .map(|s| format!("\"{}\"", display_stubs(s, false)))
            .join(","),
        neutral_losses
            .iter()
            .map(|s| format!("\"{}\"", display_neutral_loss(s, false)))
            .join(","),
        diagnostic_ions
            .iter()
            .map(|n| format!("\"{}\"", n.hill_notation()))
            .join(","),
        rules1
            .iter()
            .map(|p| display_placement_rule(p, true))
            .join(", "),
        if asymmetric {
            ", Secondary placement rules: ".to_string()
                + &rules2
                    .iter()
                    .map(|p| display_placement_rule(p, true))
                    .join(", ")
        } else {
            String::new()
        },
        if stubs.is_empty() {
            String::new()
        } else {
            ", Breakage: ".to_string() + &stubs.iter().map(|s| display_stubs(s, true)).join(", ")
        },
        if neutral_losses.is_empty() {
            String::new()
        } else {
            ", Neutral losses: ".to_string() + &neutral_losses.iter().map(|n| display_neutral_loss(n, true)).join(", ")
        },
        if diagnostic_ions.is_empty() {
            String::new()
        } else {
            ", Diagnostic ions: ".to_string()
                + &diagnostic_ions
                    .iter()
                    .map(|f| display_formula(f, true))
                    .join(", ")
        },
    ))
}

#[tauri::command]
pub fn validate_glycan_fragments(
    fallback: bool,
    aa: Vec<String>,
    kind: Vec<String>,
    full: bool,
    core: Option<(u8, u8)>,
) -> Result<String, BoxedError<'static>> {
    let aas = aa
        .iter()
        .map(|aa| parse_amino_acid(aa))
        .collect::<Result<Vec<_>, _>>()
        .map_err(BoxedError::to_owned)?;
    let kind = kind
        .iter()
        .map(|v| parse_fragment_kind(v))
        .collect::<Result<Vec<_>, _>>()
        .map_err(BoxedError::to_owned)?;

    if !fallback && aas.is_empty() {
        Err(BoxedError::error(
            "Invalid glycan fragments",
            "At least one amino acid has to be provided for the selection",
            Context::none(),
        ))
    } else if !full && core.is_none() {
        Err(BoxedError::error(
            "Invalid glycan fragments",
            "At least one fragment type has to be specified (one of full or core)",
            Context::none(),
        ))
    } else {
        Ok(format!(
            "<span data-value='[[{}], [{}], {{\"full\": {full}, \"core\": {} }}, {fallback}]'>{}, Full glycan: {full}, Core: {}",
            aas.iter().map(|aa| format!("\"{aa}\"")).join(","),
            kind.iter().map(|value| format!("\"{value}\"")).join(","),
            core.map_or("null".to_string(), |(min, max)| format!("[{min}, {max}]")),
            if fallback {
                "All undefined".to_string()
            } else {
                format!("Attachment: {}", aas.iter().join(", "),) + &if kind.is_empty() {String::new()} else {
                    format!(", Fragments: {}", kind.iter().join(", "),)
                }
            },
            core.map_or("false".to_string(), |(min, max)| format!("{min}—{max}")),
        ))
    }
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
