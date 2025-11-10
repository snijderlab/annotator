use std::sync::{Arc, Mutex};

use context_error::{BasicKind, BoxedError, Context, CreateError, FullErrorContent};
use mzcore::{
    chemistry::{DiagnosticIon, NeutralLoss},
    ontology::Ontology,
    prelude::*,
    sequence::{LinkerSpecificity, ModificationId, PlacementRule, SimpleModificationInner},
};
use mzcv::SynonymScope;
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
) -> Result<(Vec<(u32, String)>, Option<(String, Vec<String>)>), &'static str> {
    let state = state.lock().map_err(|_| "Could not lock mutex")?;

    Ok((
        state
            .ontologies
            .custom()
            .data()
            .map(|modification| {
                (
                    modification
                        .description()
                        .and_then(|d| d.id())
                        .unwrap_or_default(),
                    crate::search_modification::render_modification(
                        modification,
                        theme,
                        &state.ontologies,
                    )
                    .to_string(),
                )
            })
            .collect(),
        state.custom_modifications_error.clone(),
    ))
}

#[tauri::command]
pub fn duplicate_custom_modification(
    id: u32,
    new_id: u32,
    state: ModifiableState,
) -> Result<CustomModification, String> {
    let mut locked_state = state.lock().map_err(|_| "Could not lock mutex")?;
    if let Some(modification) = locked_state.ontologies.custom().data().find(|p| {
        p.description()
            .and_then(|d| d.id())
            .is_some_and(|i| i == id)
    }) {
        let new_modification = match modification.as_ref().clone() {
            SimpleModificationInner::Database {
                mut id,
                specificities,
                formula,
            } => {
                id.set_id(Some(new_id));
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
                id.set_id(Some(new_id));
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
                id.set_id(Some(new_id));
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
        locked_state.ontologies.custom_mut().add(new_modification);
        locked_state
            .ontologies
            .custom()
            .save_to_cache()
            .map_err(|e| e.to_html(true))?;
        drop(locked_state);
        get_custom_modification(new_id, state).map_err(|e| e.to_string())
    } else {
        Err("Could not find specified id".to_string())
    }
}

#[tauri::command]
pub fn delete_custom_modification(id: u32, state: ModifiableState) -> Result<(), &'static str> {
    let mut state = state.lock().map_err(|_| "Could not lock mutex")?;
    if state.ontologies.custom_mut().remove(&id) {
        Ok(())
    } else {
        Err("Could not find specified id")
    }
}

#[tauri::command]
pub fn get_custom_modification(
    id: u32,
    state: ModifiableState,
) -> Result<CustomModification, &'static str> {
    let state = state.lock().map_err(|_| "Could not lock mutex")?;
    match state.ontologies.custom().get_by_index(&id).as_deref() {
        Some(SimpleModificationInner::Database {
            specificities,
            formula,
            id,
        }) => Ok(CustomModification {
            id: id.id().unwrap_or_default(),
            name: id.name.to_string(),
            formula: formula.to_string(),
            description: id.description.to_string(),
            synonyms: id.synonyms.iter().map(|s| s.1.to_string()).collect(),
            cross_ids: id
                .cross_ids
                .iter()
                .map(|(a, b)| a.as_ref().map_or(b.to_string(), |a| format!("{a}:{b}")))
                .collect(),
            linker: false,
            single_specificities: specificities
                .iter()
                .map(|(rules, neutral_losses, diagnostic_ions)| {
                    (
                        rules
                            .iter()
                            .map(|p| display_placement_rule(p, false, &state.ontologies))
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
        Some(SimpleModificationInner::Linker {
            specificities,
            formula,
            id,
            length,
        }) => Ok(CustomModification {
            id: id.id().unwrap_or_default(),
            name: id.name.to_string(),
            formula: formula.to_string(),
            description: id.description.to_string(),
            synonyms: id.synonyms.iter().map(|s| s.1.to_string()).collect(),
            cross_ids: id
                .cross_ids
                .iter()
                .map(|(a, b)| a.as_ref().map_or(b.to_string(), |a| format!("{a}:{b}")))
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
                            .map(|p| display_placement_rule(p, false, &state.ontologies))
                            .collect(),
                        secondary_rules
                            .iter()
                            .map(|p| display_placement_rule(p, false, &state.ontologies))
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
                            .map(|p| display_placement_rule(p, false, &state.ontologies))
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
        Some(_) => Err("Invalid custom modification type"),
        None => Err("Given index does not exist"),
    }
}

/// To be sent back and forth with the JS side.
/// This structure is way easier to handle over there.
#[derive(Deserialize, Serialize)]
pub struct CustomModification {
    id: u32,
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
        .or_else(|_| MolecularFormula::pro_forma::<true, true>(&custom_modification.formula))
        .map_err(|e| e.to_html(false))?;
    let id = ModificationId::new(
        Ontology::Custom,
        custom_modification.name.clone().into_boxed_str(),
        Some(custom_modification.id),
        custom_modification.description.into_boxed_str(),
        custom_modification
            .synonyms
            .iter()
            .map(|s| (SynonymScope::Exact, s.clone().into_boxed_str()))
            .collect(),
        custom_modification
            .cross_ids
            .iter()
            .map(|id| {
                id.split_once(':')
                    .map_or((None, id.clone().into_boxed_str()), |(a, b)| {
                        (Some(a.into()), b.into())
                    })
            })
            .collect::<Vec<_>>()
            .into(),
    );
    let modification = Arc::new(if custom_modification.linker {
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
                                    .or_else(|_| MolecularFormula::pro_forma::<true, true>(d))
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
                .map_err(|e| e.to_html(false))?,
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
                                    .or_else(|_| MolecularFormula::pro_forma::<true, true>(d))
                                    .map(DiagnosticIon)
                            })
                            .collect::<Result<Vec<_>, _>>()?,
                    ))
                })
                .collect::<Result<Vec<_>, BoxedError<'_, BasicKind>>>()
                .map_err(|e| e.to_html(false))?,
            formula,
            id,
        }
    });

    if let Ok(mut state) = app.state::<Mutex<State>>().lock() {
        // Update state (if this modification was not stored before the remove will just do nothing)
        let custom = state.ontologies.custom_mut();
        custom.remove(&custom_modification.id);
        custom.add(modification);
        custom.save_to_cache().map_err(|e| e.to_html(true))?;

        Ok(())
    } else {
        Err(BoxedError::new(
            BasicKind::Error,
            "State locked",
            "Cannot unlock the mutable state, are you doing many things in parallel?",
            Context::none(),
        )
        .to_html(false))
    }
}
