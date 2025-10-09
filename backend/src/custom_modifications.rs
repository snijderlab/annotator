use std::{
    io::BufWriter,
    sync::{Arc, Mutex},
};

use context_error::{BasicKind, BoxedError, Context, CreateError, FullErrorContent};
use mzcore::{
    chemistry::{DiagnosticIon, NeutralLoss},
    ontology::Ontology,
    prelude::*,
    sequence::{LinkerSpecificity, ModificationId, PlacementRule, SimpleModificationInner},
};
use ordered_float::OrderedFloat;
use serde::{Deserialize, Serialize};
use tauri::Manager;

use crate::{
    ModifiableState, Theme,
    render::{display_neutral_loss, display_placement_rule, display_stubs},
    state::State,
    validate::parse_stub,
};

#[tauri::command]
pub fn get_custom_modifications(
    state: ModifiableState,
    theme: Theme,
) -> Result<(Vec<(usize, String)>, Option<(String, Vec<String>)>), &'static str> {
    let state = state.lock().map_err(|_| "Could not lock mutex")?;
    Ok((
        state
            .custom_modifications
            .iter()
            .map(|(index, _, modification)| {
                (
                    index.unwrap_or_default(),
                    crate::search_modification::render_modification(modification, theme)
                        .to_string(),
                )
            })
            .collect(),
        state.custom_modifications_error.clone(),
    ))
}

#[tauri::command]
pub fn duplicate_custom_modification(
    id: usize,
    new_id: usize,
    state: ModifiableState,
) -> Result<CustomModification, &'static str> {
    let mut locked_state = state.lock().map_err(|_| "Could not lock mutex")?;
    if let Some(index) = locked_state
        .custom_modifications
        .iter()
        .position(|p| p.0.is_some_and(|i| i == id))
    {
        let mut modification = locked_state.custom_modifications[index].clone();
        modification.0 = Some(new_id);
        modification.2 = match modification.2.as_ref().clone() {
            SimpleModificationInner::Database {
                mut id,
                specificities,
                formula,
            } => {
                id.id = Some(new_id);
                SimpleModificationInner::Database {
                    specificities,
                    formula,
                    id,
                }
            }
            SimpleModificationInner::Gno {
                mut id,
                composition,
                structure_score,
                subsumption_level,
                motif,
                taxonomy,
                glycomeatlas,
            } => {
                id.id = Some(new_id);
                SimpleModificationInner::Gno {
                    id,
                    composition,
                    structure_score,
                    subsumption_level,
                    motif,
                    taxonomy,
                    glycomeatlas,
                }
            }
            SimpleModificationInner::Linker {
                mut id,
                specificities,
                formula,
                length,
            } => {
                id.id = Some(new_id);
                SimpleModificationInner::Linker {
                    id,
                    specificities,
                    formula,
                    length,
                }
            }
            full => full,
        }
        .into();
        locked_state.custom_modifications.push(modification);
        drop(locked_state);
        get_custom_modification(new_id, state)
    } else {
        Err("Could not find specified id")
    }
}

#[tauri::command]
pub fn delete_custom_modification(id: usize, state: ModifiableState) -> Result<(), &'static str> {
    let mut state = state.lock().map_err(|_| "Could not lock mutex")?;
    if let Some(index) = state
        .custom_modifications
        .iter()
        .position(|p| p.0.is_some_and(|i| i == id))
    {
        state.custom_modifications.remove(index);
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
        .custom_modifications
        .iter()
        .position(|p| p.0.is_some_and(|i| i == id))
    {
        match &*state.custom_modifications[index].2 {
            SimpleModificationInner::Database {
                specificities,
                formula,
                id,
            } => Ok(CustomModification {
                id: id.id.unwrap_or_default(),
                name: id.name.clone(),
                formula: formula.to_string(),
                description: id.description.clone(),
                synonyms: id.synonyms.to_vec(),
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
            SimpleModificationInner::Linker {
                specificities,
                formula,
                id,
                length,
            } => Ok(CustomModification {
                id: id.id.unwrap_or_default(),
                name: id.name.clone(),
                formula: formula.to_string(),
                description: id.description.clone(),
                synonyms: id.synonyms.to_vec(),
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
                        LinkerSpecificity::Asymmetric {
                            rules: (rules, secondary_rules),
                            stubs,
                            neutral_losses,
                            diagnostic,
                        } => (
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
                            neutral_losses
                                .iter()
                                .map(|n| display_neutral_loss(n, false))
                                .collect(),
                            diagnostic.iter().map(|n| n.0.hill_notation()).collect(),
                        ),
                        LinkerSpecificity::Symmetric {
                            rules,
                            stubs,
                            neutral_losses,
                            diagnostic,
                        } => (
                            false,
                            rules
                                .iter()
                                .map(|p| display_placement_rule(p, false))
                                .collect(),
                            Vec::new(),
                            stubs.iter().map(|s| display_stubs(s, false)).collect(),
                            neutral_losses
                                .iter()
                                .map(|n| display_neutral_loss(n, false))
                                .collect(),
                            diagnostic.iter().map(|n| n.0.hill_notation()).collect(),
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
    linker_specificities: Vec<(
        bool,
        Vec<String>,
        Vec<String>,
        Vec<String>,
        Vec<String>,
        Vec<String>,
    )>,
    linker_length: Option<f64>,
}

#[tauri::command]
pub async fn update_modification(
    custom_modification: CustomModification,
    app: tauri::AppHandle,
) -> Result<(), String> {
    let formula = custom_modification
        .formula
        .parse::<f64>()
        .map(MolecularFormula::with_additional_mass)
        .or_else(|_| {
            MolecularFormula::from_pro_forma(
                &custom_modification.formula,
                ..,
                true,
                true,
                true,
                false,
            )
        })
        .map_err(|e| e.to_html())?;
    let id = ModificationId {
        ontology: Ontology::Custom,
        name: custom_modification.name.clone(),
        id: Some(custom_modification.id),
        description: custom_modification.description,
        synonyms: custom_modification.synonyms.into(),
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
        Arc::new(if custom_modification.linker {
            SimpleModificationInner::Linker {
                specificities: custom_modification
                    .linker_specificities
                    .iter()
                    .map(
                        |(
                            asymmetric,
                            placement_rules,
                            secondary_placement_rules,
                            stubs,
                            neutral_losses,
                            diagnostic,
                        )| {
                            let stubs = stubs
                                .iter()
                                .map(|text| parse_stub(text))
                                .collect::<Result<Vec<_>, _>>()?;
                            let neutral_losses = neutral_losses
                                .iter()
                                .map(|text| text.parse())
                                .collect::<Result<Vec<_>, _>>()?;
                            let diagnostic = diagnostic
                                .iter()
                                .map(|d| {
                                    d.parse::<f64>()
                                        .map(MolecularFormula::with_additional_mass)
                                        .or_else(|_| {
                                            MolecularFormula::from_pro_forma(
                                                d,
                                                ..,
                                                true,
                                                true,
                                                true,
                                                false,
                                            )
                                        })
                                        .map(DiagnosticIon)
                                })
                                .collect::<Result<Vec<_>, _>>()?;
                            if *asymmetric {
                                Ok(LinkerSpecificity::Asymmetric {
                                    rules: (
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
                                    neutral_losses,
                                    diagnostic,
                                })
                            } else {
                                Ok(LinkerSpecificity::Symmetric {
                                    rules: placement_rules
                                        .iter()
                                        .map(|r| r.parse::<PlacementRule>())
                                        .collect::<Result<Vec<_>, _>>()?,
                                    stubs,
                                    neutral_losses,
                                    diagnostic,
                                })
                            }
                        },
                    )
                    .collect::<Result<Vec<_>, BoxedError<'_, BasicKind>>>()
                    .map_err(|e| e.to_html())?,
                formula,
                id,
                length: custom_modification.linker_length.map(OrderedFloat::from),
            }
        } else {
            SimpleModificationInner::Database {
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
                                            MolecularFormula::from_pro_forma(
                                                d,
                                                ..,
                                                true,
                                                true,
                                                true,
                                                false,
                                            )
                                        })
                                        .map(DiagnosticIon)
                                })
                                .collect::<Result<Vec<_>, _>>()?,
                        ))
                    })
                    .collect::<Result<Vec<_>, BoxedError<'_, BasicKind>>>()
                    .map_err(|e| e.to_html())?,
                formula,
                id,
            }
        }),
    );

    if let Ok(mut state) = app.state::<Mutex<State>>().lock() {
        // Update state
        if let Some(index) = state
            .custom_modifications
            .iter()
            .position(|p| p.0 == modification.0)
        {
            state.custom_modifications[index] = modification;
        } else {
            state.custom_modifications.push(modification);
        }

        // Store mods config file
        let path = app
            .path()
            .app_config_dir()
            .map(|dir| dir.join(crate::CUSTOM_MODIFICATIONS_FILE))
            .map_err(|e| {
                BoxedError::new(
                    BasicKind::Error,
                    "Cannot find app data directory",
                    e.to_string(),
                    Context::none(),
                )
                .to_html()
            })?;
        let parent = path.parent().ok_or_else(|| {
            BoxedError::new(
                BasicKind::Error,
                "Custom modifications configuration does not have a valid directory",
                "Please report",
                Context::show(path.to_string_lossy()).to_owned(),
            )
            .to_html()
        })?;
        std::fs::create_dir_all(parent).map_err(|err| {
            BoxedError::new(
                BasicKind::Error,
                "Could not create parent directories for custom modifications configuration file",
                err.to_string(),
                Context::show(parent.to_string_lossy()).to_owned(),
            )
            .to_html()
        })?;
        let file = BufWriter::new(std::fs::File::create(&path).map_err(|err| {
            BoxedError::new(
                BasicKind::Error,
                "Could not open custom modifications configuration file",
                err.to_string(),
                Context::show(path.to_string_lossy()).to_owned(),
            )
            .to_html()
        })?);
        serde_json::to_writer_pretty(file, &state.custom_modifications).map_err(|err| {
            BoxedError::new(
                BasicKind::Error,
                "Could not write custom modifications to configuration file",
                err.to_string(),
                Context::none(),
            )
            .to_html()
        })?;

        Ok(())
    } else {
        Err(BoxedError::new(
            BasicKind::Error,
            "State locked",
            "Cannot unlock the mutable state, are you doing many things in parallel?",
            Context::none(),
        )
        .to_html())
    }
}
