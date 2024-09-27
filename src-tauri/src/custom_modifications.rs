use std::{io::BufWriter, sync::Mutex};

use itertools::Itertools;
use ordered_float::OrderedFloat;
use rustyms::{
    error::{Context, CustomError},
    modification::{LinkerSpecificity, ModificationId, Ontology, SimpleModification},
    placement_rule::PlacementRule,
    DiagnosticIon, MolecularFormula, NeutralLoss,
};
use serde::{Deserialize, Serialize};
use tauri::Manager;

use crate::{
    render::{display_formula, display_neutral_loss, display_placement_rule, display_stubs},
    state::State,
    ModifiableState,
};

#[tauri::command]
pub fn validate_molecular_formula(text: String) -> Result<String, CustomError> {
    text.parse::<f64>()
        .map(MolecularFormula::with_additional_mass)
        .or_else(|_| MolecularFormula::from_pro_forma(&text, .., false, true))
        .map(|f| display_formula(&f, true))
}

#[tauri::command]
pub fn validate_neutral_loss(text: String) -> Result<String, CustomError> {
    text.parse::<NeutralLoss>()
        .map(|f| display_neutral_loss(&f))
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
            .or_else(|_| MolecularFormula::from_pro_forma(text, ..index, false, true))?;
        let f2 = text[index + 1..]
            .parse::<f64>()
            .map(MolecularFormula::with_additional_mass)
            .or_else(|_| MolecularFormula::from_pro_forma(text, index + 1.., false, true))?;
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
                .or_else(|_| MolecularFormula::from_pro_forma(&text, .., false, true))
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
                .or_else(|_| MolecularFormula::from_pro_forma(&text, .., false, true))
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
            ", Neutral losses: ".to_string() + &stubs.iter().map(|s| display_stubs(s, true)).join(", ")
        },
        if diagnostic_ions.is_empty() {
            String::new()
        } else {
            ", Diagnostic ions: ".to_string() + &diagnostic_ions.iter().map(|f| display_formula(f, true)).join(", ")
        },
    ))
}

#[tauri::command]
pub fn get_custom_modifications(
    state: ModifiableState,
) -> Result<Vec<(usize, String)>, &'static str> {
    let state = state.lock().map_err(|_| "Could not lock mutex")?;
    Ok(state
        .database
        .iter()
        .map(|(index, _, modification)| {
            (
                index.unwrap_or_default(),
                crate::search_modification::render_modification(modification).to_string(),
            )
        })
        .collect())
}

#[tauri::command]
pub fn delete_custom_modification(id: usize, state: ModifiableState) -> Result<(), &'static str> {
    let mut state = state.lock().map_err(|_| "Could not lock mutex")?;
    if let Some(index) = state
        .database
        .iter()
        .position(|p| p.0.is_some_and(|i| i == id))
    {
        state.database.remove(index);
        Ok(())
    } else {
        Err("Could not find specified id")
    }
}

#[tauri::command]
pub fn get_custom_modification(
    id: usize,
    state: ModifiableState,
) -> Result<CustomModification, &'static str> {
    let state = state.lock().map_err(|_| "Could not lock mutex")?;
    if let Some(index) = state
        .database
        .iter()
        .position(|p| p.0.is_some_and(|i| i == id))
    {
        match &state.database[index].2 {
            SimpleModification::Database {
                specificities,
                formula,
                id,
            } => Ok(CustomModification {
                id: id.id.unwrap_or_default(),
                name: id.name.clone(),
                formula: formula.to_string(),
                description: id.description.clone(),
                synonyms: id.synonyms.clone(),
                cross_ids: id
                    .cross_ids
                    .iter()
                    .map(|(a, b)| format!("{a}:{b}"))
                    .collect(),
                linker: false,
                single_specificities: specificities
                    .iter()
                    .map(|(rules, neutral_losses, diagnostic_ions)| {
                        (
                            rules
                                .iter()
                                .map(|p| display_placement_rule(p, false))
                                .collect(),
                            neutral_losses.iter().map(|n| n.hill_notation()).collect(),
                            diagnostic_ions
                                .iter()
                                .map(|n| n.0.hill_notation())
                                .collect(),
                        )
                    })
                    .collect(),
                linker_specificities: Vec::new(),
                linker_length: None,
            }),
            SimpleModification::Linker {
                specificities,
                formula,
                id,
                length,
            } => Ok(CustomModification {
                id: id.id.unwrap_or_default(),
                name: id.name.clone(),
                formula: formula.to_string(),
                description: id.description.clone(),
                synonyms: id.synonyms.clone(),
                cross_ids: id
                    .cross_ids
                    .iter()
                    .map(|(a, b)| format!("{a}:{b}"))
                    .collect(),
                linker: true,
                single_specificities: Vec::new(),
                linker_specificities: specificities
                    .iter()
                    .map(|spec| match spec {
                        LinkerSpecificity::Asymmetric(
                            (rules, secondary_rules),
                            stubs,
                            diagnostic_ions,
                        ) => (
                            true,
                            rules
                                .iter()
                                .map(|p| display_placement_rule(p, false))
                                .collect(),
                            secondary_rules
                                .iter()
                                .map(|p| display_placement_rule(p, false))
                                .collect(),
                            stubs.iter().map(|s| display_stubs(s, false)).collect(),
                            diagnostic_ions
                                .iter()
                                .map(|n| n.0.hill_notation())
                                .collect(),
                        ),
                        LinkerSpecificity::Symmetric(rules, stubs, diagnostic_ions) => (
                            false,
                            rules
                                .iter()
                                .map(|p| display_placement_rule(p, false))
                                .collect(),
                            Vec::new(),
                            stubs.iter().map(|s| display_stubs(s, false)).collect(),
                            diagnostic_ions
                                .iter()
                                .map(|n| n.0.hill_notation())
                                .collect(),
                        ),
                    })
                    .collect(),
                linker_length: length.map(|n| n.0),
            }),
            _ => Err("Invalid custom modification type"),
        }
    } else {
        Err("Given index does not exist")
    }
}

/// To be sent back and forth with the JS side.
/// This structure is way easier to handle over there.
#[derive(Deserialize, Serialize)]
pub struct CustomModification {
    id: usize,
    name: String,
    formula: String,
    description: String,
    synonyms: Vec<String>,
    cross_ids: Vec<String>,
    linker: bool,
    single_specificities: Vec<(Vec<String>, Vec<String>, Vec<String>)>,
    linker_specificities: Vec<(bool, Vec<String>, Vec<String>, Vec<String>, Vec<String>)>,
    linker_length: Option<f64>,
}

#[tauri::command]
pub async fn update_modification(
    custom_modification: CustomModification,
    app: tauri::AppHandle,
) -> Result<(), CustomError> {
    let formula = custom_modification
        .formula
        .parse::<f64>()
        .map(MolecularFormula::with_additional_mass)
        .or_else(|_| {
            MolecularFormula::from_pro_forma(&custom_modification.formula, .., false, true)
        })?;
    let id = ModificationId {
        ontology: Ontology::Custom,
        name: custom_modification.name.clone(),
        id: Some(custom_modification.id),
        description: custom_modification.description,
        synonyms: custom_modification.synonyms,
        cross_ids: custom_modification
            .cross_ids
            .iter()
            .filter_map(|id| id.split_once(':'))
            .map(|(a, b)| (a.to_string(), b.to_string()))
            .collect(),
    };
    let modification = (
        Some(custom_modification.id),
        custom_modification.name.to_lowercase(),
        if custom_modification.linker {
            SimpleModification::Linker {
                specificities: custom_modification
                    .linker_specificities
                    .iter()
                    .map(
                        |(
                            asymmetric,
                            placement_rules,
                            secondary_placement_rules,
                            stubs,
                            diagnostic_ions,
                        )| {
                            let stubs = stubs
                                .iter()
                                .map(|text| parse_stub(text))
                                .collect::<Result<Vec<_>, _>>()?;
                            let diagnostic_ions = diagnostic_ions
                                .iter()
                                .map(|d| {
                                    d.parse::<f64>()
                                        .map(MolecularFormula::with_additional_mass)
                                        .or_else(|_| {
                                            MolecularFormula::from_pro_forma(d, .., false, true)
                                        })
                                        .map(DiagnosticIon)
                                })
                                .collect::<Result<Vec<_>, _>>()?;
                            if *asymmetric {
                                Ok(LinkerSpecificity::Asymmetric(
                                    (
                                        placement_rules
                                            .iter()
                                            .map(|r| r.parse::<PlacementRule>())
                                            .collect::<Result<Vec<_>, _>>()?,
                                        secondary_placement_rules
                                            .iter()
                                            .map(|r| r.parse::<PlacementRule>())
                                            .collect::<Result<Vec<_>, _>>()?,
                                    ),
                                    stubs,
                                    diagnostic_ions,
                                ))
                            } else {
                                Ok(LinkerSpecificity::Symmetric(
                                    placement_rules
                                        .iter()
                                        .map(|r| r.parse::<PlacementRule>())
                                        .collect::<Result<Vec<_>, _>>()?,
                                    stubs,
                                    diagnostic_ions,
                                ))
                            }
                        },
                    )
                    .collect::<Result<Vec<_>, _>>()?,
                formula,
                id,
                length: custom_modification.linker_length.map(OrderedFloat::from),
            }
        } else {
            SimpleModification::Database {
                specificities: custom_modification
                    .single_specificities
                    .iter()
                    .map(|(placement_rules, neutral_losses, diagnostic_ions)| {
                        Ok((
                            placement_rules
                                .iter()
                                .map(|r| r.parse::<PlacementRule>())
                                .collect::<Result<Vec<_>, _>>()?,
                            neutral_losses
                                .iter()
                                .map(|n| n.parse::<NeutralLoss>())
                                .collect::<Result<Vec<_>, _>>()?,
                            diagnostic_ions
                                .iter()
                                .map(|d| {
                                    d.parse::<f64>()
                                        .map(MolecularFormula::with_additional_mass)
                                        .or_else(|_| {
                                            MolecularFormula::from_pro_forma(d, .., false, true)
                                        })
                                        .map(DiagnosticIon)
                                })
                                .collect::<Result<Vec<_>, _>>()?,
                        ))
                    })
                    .collect::<Result<Vec<_>, _>>()?,
                formula,
                id,
            }
        },
    );

    if let Ok(mut state) = app.state::<Mutex<State>>().lock() {
        // Update state
        if let Some(index) = state.database.iter().position(|p| p.0 == modification.0) {
            state.database[index] = modification;
        } else {
            state.database.push(modification);
        }

        // Store mods config file
        let path = app
            .path_resolver()
            .app_config_dir()
            .map(|dir| dir.join(crate::CUSTOM_MODIFICATIONS_FILE))
            .ok_or(CustomError::error(
                "Cannot find app data directory",
                "",
                Context::None,
            ))?;
        let parent = path.parent().ok_or_else(|| {
            CustomError::error(
                "Custom modifications configuration does not have a valid directory",
                "Please report",
                Context::show(path.to_string_lossy()),
            )
        })?;
        std::fs::create_dir_all(parent).map_err(|err| {
            CustomError::error(
                "Could not create parent directories for custom modifications configuration file",
                err,
                Context::show(parent.to_string_lossy()),
            )
        })?;
        let file = BufWriter::new(std::fs::File::create(&path).map_err(|err| {
            CustomError::error(
                "Could not open custom modifications configuration file",
                err,
                Context::show(path.to_string_lossy()),
            )
        })?);
        serde_json::to_writer_pretty(file, &state.database).map_err(|err| {
            CustomError::error(
                "Could not write custom modifications to configuration file",
                err,
                Context::None,
            )
        })?;

        Ok(())
    } else {
        Err(CustomError::error(
            "State locked",
            "Cannot unlock the mutable state, are you doing many things in parallel?",
            Context::None,
        ))
    }
}
