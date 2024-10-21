use crate::{
    html_builder,
    metadata_render::OptionalString,
    render::{display_formula, display_masses, display_placement_rule},
    ModifiableState,
};
use itertools::Itertools;
use modification::ModificationId;
use rustyms::{
    error::*,
    modification::{GnoComposition, LinkerSpecificity, SimpleModification},
    placement_rule::PlacementRule,
    system::{dalton, Mass},
    ReturnModification, *,
};

use crate::{
    html_builder::{HtmlElement, HtmlTag},
    render::{display_neutral_loss, display_stubs, link_modification},
};

#[tauri::command]
pub async fn search_modification(
    text: &str,
    tolerance: f64,
    state: ModifiableState<'_>,
) -> Result<String, CustomError> {
    let state = state.lock().unwrap();
    let modification = if text.is_empty() {
        Err(CustomError::error(
            "Invalid modification",
            "The modification is empty",
            Context::None,
        ))
    } else {
        SimpleModification::try_from(
            text,
            0..text.len(),
            &mut Vec::new(),
            &mut Vec::new(),
            Some(&state.database),
        )
        .map(|m| match m {
            ReturnModification::Defined(d) => Ok(d),
            _ => Err(CustomError::error(
                "Invalid modification",
                "Can not define ambiguous modifications for the modifications parameter",
                Context::None,
            )),
        })
    }??;
    let tolerance = Tolerance::new_absolute(Mass::new::<dalton>(tolerance));

    match modification {
        SimpleModification::Mass(m)
        | SimpleModification::Gno {
            composition: GnoComposition::Weight(m),
            ..
        } => Ok(html_builder::HtmlElement::table(
            Some(&[
                "Name".to_string(),
                "Id".to_string(),
                MassMode::Monoisotopic.to_string(),
                "Formula".to_string(),
            ]),
            modification_search_mass(
                m.into_inner(),
                tolerance,
                None,
                MassMode::Monoisotopic,
                Some(&state.database),
            )
            .map(|(ontology, id, name, modification)| {
                [
                    modification.to_string(),
                    link_modification(ontology, id, &name),
                    display_masses(&modification.formula()).to_string(),
                    display_formula(&modification.formula(), true),
                ]
            }),
        )
        .to_string()),
        SimpleModification::Formula(f) => Ok(html_builder::HtmlElement::table(
            Some(&["Name", "Id"]),
            modification_search_formula(&f, Some(&state.database)).map(
                |(ontology, id, name, modification)| {
                    [
                        modification.to_string(),
                        link_modification(ontology, id, &name),
                    ]
                },
            ),
        )
        .to_string()),
        SimpleModification::Glycan(ref g)
        | SimpleModification::Gno {
            composition: GnoComposition::Composition(ref g),
            ..
        } => Ok(HtmlTag::div
            .new()
            .content(render_modification(&modification))
            .content(html_builder::HtmlElement::table(
                Some(&["Name", "Structure"]),
                modification_search_glycan(g, true).map(|(ontology, id, name, modification)| {
                    [
                        link_modification(ontology, id, &name),
                        if let SimpleModification::Gno {
                            composition: GnoComposition::Topology(structure),
                            ..
                        } = modification
                        {
                            structure.to_string()
                        } else {
                            "-".to_string()
                        },
                    ]
                }),
            ))
            .to_string()),
        modification => Ok(render_modification(&modification).to_string()),
    }
}

pub fn render_modification(modification: &SimpleModification) -> HtmlElement {
    let mut output = HtmlTag::div.new();
    output.class("modification");

    output.content(HtmlTag::p.new().content(format!(
        "Formula {} Mass {}",
        display_formula(&modification.formula(), true),
        display_masses(&modification.formula()),
    )));

    if let modification::SimpleModification::Database {
        specificities, id, ..
    } = &modification
    {
        output.content(render_modification_id(id));
        output.content(HtmlTag::p.new().content("Placement rules"));
        let mut ul = HtmlTag::ul.new();

        for rule in specificities {
            ul.content(HtmlTag::li.new().content(format!(
                "Positions: {}{}{}{}{}",
                render_places(&rule.0),
                if rule.1.is_empty() {
                    ""
                } else {
                    ", Neutral losses: "
                },
                rule.1.iter().map(display_neutral_loss).join(", "),
                if rule.2.is_empty() {
                    ""
                } else {
                    ", Diagnostic ions: "
                },
                rule.2.iter().map(|d| &d.0).map(|f| display_formula(f, true)).join(", ")
            )));
        }
        output.content(ul);
    } else if let SimpleModification::Gno {
        composition,
        id,
        structure_score,
        subsumption_level,
        motif,
        taxonomy,
        glycomeatlas,
    } = &modification
    {
        output.content(render_modification_id(id));
        output.content(
            HtmlTag::p
                .new()
                .content(format!("Subsumption level: {subsumption_level}")),
        );
        output.maybe_content(structure_score.map(|score| {
            HtmlTag::p
                .new()
                .content(format!("Structure score: {score}"))
                .clone()
        }));
        output.maybe_content(motif.as_ref().map(|(name, id)| {
            HtmlTag::p
                .new()
                .content(format!(
                    "Motif: {name} {}",
                    link_modification(modification::Ontology::Gnome, None, id)
                ))
                .clone()
        }));
        output.content(HtmlTag::p.new().content(format!(
                "Taxonomy: {}",
                taxonomy
                    .iter()
                    .map(|(name, id)| format!("<span title='ID: {id}'>{name}</span>"))
                    .join(", ")
            )));

        output.content(HtmlTag::p.new().content("Glycomeatlas"));
        output.content(
            HtmlTag::ul
                .new()
                .children(glycomeatlas.iter().map(|(species, places)| {
                    HtmlTag::li
                        .new()
                        .content(format!("{species}: "))
                        .children(places.iter().map(|(name, id)| {
                            HtmlTag::span
                                .new()
                                .title(format!("ID: {id}"))
                                .content(name)
                                .clone()
                        }))
                        .clone()
                }))
                .clone(),
        );

        match composition {
            GnoComposition::Weight(_) => {
                output.content(HtmlTag::p.new().content("Only weight known"));
            }
            GnoComposition::Composition(composition) => {
                output.content(HtmlTag::p.new().content(format!(
                        "Composition: {}",
                        composition
                            .iter()
                            .fold(String::new(), |acc, m| acc + &format!("{}{}", m.0, m.1))
                    )));
            }
            GnoComposition::Topology(structure) => {
                output.extend([
                    HtmlTag::p.new().content(format!(
                        "Composition: {}",
                        structure
                            .composition()
                            .iter()
                            .fold(String::new(), |acc, m| acc + &format!("{}{}", m.0, m.1))
                    )),
                    HtmlTag::p.new().content(format!("Structure: {structure}")),
                ]);
            }
        }
    } else if let SimpleModification::Linker {
        specificities,
        id,
        length,
        ..
    } = &modification
    {
        output.content(render_modification_id(id));
        output.content(HtmlTag::p.new().content(format!(
            "Length: {}",
            length.map_or("-".to_string(), |l| l.to_string()),
        )));
        output.content(HtmlTag::p.new().content("Placement rules"));
        let mut ul = HtmlTag::ul.new();
        ul.class("placement-rules");

        for specificity in specificities {
            match specificity {
                LinkerSpecificity::Symmetric(places, stubs, diagnostic_ions) => {
                    ul.content(
                        HtmlTag::li
                            .new()
                            .content(format!("Positions: {}", render_places(places)))
                            .maybe_content((!stubs.is_empty()).then(|| {
                                format!(
                                    ", Breakages: {}",
                                    stubs.iter().map(|s| display_stubs(s, true)).join(", ")
                                )
                            }))
                            .maybe_content((!diagnostic_ions.is_empty()).then(|| {
                                format!(
                                    ", Diagnostic ions: {}",
                                    diagnostic_ions
                                        .iter()
                                        .map(|d| display_formula(&d.0, true))
                                        .join(", ")
                                )
                            })),
                    );
                }
                LinkerSpecificity::Asymmetric(
                    (left_places, right_places),
                    stubs,
                    diagnostic_ions,
                ) => {
                    ul.content(
                        HtmlTag::li
                            .new()
                            .content(format!(
                                "First positions: {}, Second positions: {}",
                                render_places(left_places),
                                render_places(right_places)
                            ))
                            .maybe_content((!stubs.is_empty()).then(|| {
                                format!(
                                    ", Breakages: {}",
                                    stubs.iter().map(|s| display_stubs(s, true)).join(", ")
                                )
                            }))
                            .maybe_content((!diagnostic_ions.is_empty()).then(|| {
                                format!(
                                    ", Diagnostic ions: {}",
                                    diagnostic_ions
                                        .iter()
                                        .map(|d| display_formula(&d.0, true))
                                        .join(", ")
                                )
                            })),
                    );
                }
            }
        }
        output.content(ul);
    }

    output
}

fn render_modification_id(id: &ModificationId) -> HtmlElement {
    let mut text = HtmlTag::div.new();
    text
        .content(HtmlTag::p.new().content(format!(
            "Ontology: <span class='ontology'>{}</span>, name: <span class='name'>{}</span>, index: <span class='index'>{}</span>{}",
            id.ontology, id.name, id.id.to_optional_string(), if let Some(url) = id.url() {
                format!(", {}", HtmlTag::a.new()
                .content("view online")
                .header("href", url)
                .header("target", "_blank"))
            } else {
                String::new()
            }
        )))
        .maybe_content((!id.description.is_empty()).then_some(HtmlTag::p.new().class("description").content(&id.description)));
    if !id.synonyms.is_empty() {
        text.content(
            HtmlTag::ul
                .new()
                .class("synonyms")
                .children([HtmlTag::li.new().class("title").content("Synonyms")])
                .children(
                    id.synonyms
                        .iter()
                        .map(|syn| HtmlTag::li.new().content(syn).clone()),
                ),
        );
    }
    if !id.cross_ids.is_empty() {
        text.content(
            HtmlTag::ul
                .new()
                .class("cross-ids")
                .children([HtmlTag::li.new().class("title").content("Cross IDs")])
                .children(id.cross_ids.iter().map(|(meta, id)| {
                    HtmlTag::li
                        .new()
                        .content(format!("<span class='cross-id-system'>{meta}</span>: {id}"))
                        .clone()
                })),
        );
    }
    text
}

fn render_places(places: &[PlacementRule]) -> String {
    if places.is_empty() {
        "Anywhere".to_string()
    } else {
        places
            .iter()
            .map(|p| display_placement_rule(p, true))
            .join(", ")
    }
}
