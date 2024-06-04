use crate::{
    html_builder,
    render::{display_formula, display_mass, display_placement_rule},
};
use itertools::Itertools;
use modification::ModificationId;
use rustyms::{
    error::*,
    modification::{
        GnoComposition, LinkerSpecificity, ModificationSearchResult, Ontology, SimpleModification,
    },
    placement_rule::PlacementRule,
    system::{dalton, Mass},
    ReturnModification, *,
};

use crate::{
    html_builder::{HtmlElement, HtmlTag},
    render::{display_neutral_loss, display_stubs, link_modification},
};

#[tauri::command]
pub async fn search_modification(text: &str, tolerance: f64) -> Result<String, CustomError> {
    fn render_places(places: &[PlacementRule]) -> String {
        if places.is_empty() {
            "Anywhere".to_string()
        } else {
            places.iter().map(display_placement_rule).join(",")
        }
    }

    let modification = if text.is_empty() {
        Err(CustomError::error(
            "Invalid modification",
            "The modification is empty",
            Context::None,
        ))
    } else {
        SimpleModification::try_from(text, 0..text.len(), &mut Vec::new(), &mut Vec::new(), None)
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
    let result = SimpleModification::search(&modification, tolerance, None);

    match result {
        ModificationSearchResult::Single(modification) => {
            let mut output = HtmlElement::new(HtmlTag::div);

            output = output.content(HtmlElement::new(HtmlTag::p).content(format!(
                "Formula {} monoisotopic mass {} average mass {} most abundant mass {}",
                display_formula(&modification.formula()),
                display_mass(modification.formula().monoisotopic_mass()),
                display_mass(modification.formula().average_weight()),
                display_mass(modification.formula().most_abundant_mass()),
            )));

            if let modification::SimpleModification::Database {
                specificities, id, ..
            } = &modification
            {
                output = output.content(render_modification_id(id));
                output = output.content(HtmlElement::new(HtmlTag::p).content("Placement rules"));
                let mut ul = HtmlElement::new(HtmlTag::ul);

                for rule in specificities {
                    ul = ul.content(HtmlElement::new(HtmlTag::li).content(format!(
                        "Positions: {}{}{}{}{}",
                        render_places(&rule.0),
                        if rule.1.is_empty() {
                            ""
                        } else {
                            ", Neutral losses: "
                        },
                        rule.1.iter().map(display_neutral_loss).join(","),
                        if rule.2.is_empty() {
                            ""
                        } else {
                            ", Diagnostic ions: "
                        },
                        rule.2.iter().map(|d| &d.0).map(display_formula).join(",")
                    )));
                }
                output = output.content(ul);
            } else if let SimpleModification::Gno(composition, name) = &modification {
                output = output.content(
                    HtmlElement::new(HtmlTag::p)
                        .content(format!("Ontology: Gnome, name: {name}, "))
                        .maybe_content(modification.ontology_url().map(|url| {
                            HtmlElement::new(HtmlTag::a)
                                .content("ontology link")
                                .header("href", url)
                                .header("target", "_blank")
                        })),
                );
                match composition {
                    GnoComposition::Mass(_) => {
                        output =
                            output.content(HtmlElement::new(HtmlTag::p).content("Only mass known"));
                    }
                    GnoComposition::Structure(structure) => {
                        output.extend([
                            HtmlElement::new(HtmlTag::p).content(format!(
                                "Composition: {}",
                                structure
                                    .composition()
                                    .iter()
                                    .fold(String::new(), |acc, m| acc + &format!("{}{}", m.0, m.1))
                            )),
                            HtmlElement::new(HtmlTag::p).content(format!("Structure: {structure}")),
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
                output = output.content(render_modification_id(id));
                output = output.content(HtmlElement::new(HtmlTag::p).content(format!(
                    "Length: {}",
                    length.map_or("-".to_string(), |l| l.to_string()),
                )));
                output = output.content(HtmlElement::new(HtmlTag::p).content("Placement rules"));
                let mut ul = HtmlElement::new(HtmlTag::ul);

                for specificity in specificities {
                    match specificity {
                        LinkerSpecificity::Symmetric(places, stubs, diagnostic_ions) => {
                            ul = ul.content(HtmlElement::new(HtmlTag::li).content(format!(
                                    "Positions: {}{}{}{}{}",
                                    render_places(places),
                                    if stubs.is_empty() { "" } else { ", stubs: " },
                                    stubs.iter().map(display_stubs).join(","),
                                    if diagnostic_ions.is_empty() {
                                        ""
                                    } else {
                                        ", diagnostic ions: "
                                    },
                                    diagnostic_ions
                                        .iter()
                                        .map(|d| &d.0)
                                        .map(display_formula)
                                        .join(",")
                                )));
                        }
                        LinkerSpecificity::Asymmetric(
                            (left_places, right_places),
                            stubs,
                            diagnostic_ions,
                        ) => {
                            ul = ul
                                .content(HtmlElement::new(HtmlTag::li).content(format!(
                                    "Positions: {}{}{}<br>",
                                    render_places(left_places),
                                    if stubs.is_empty() { "" } else { ", stubs: " },
                                    stubs.iter().map(|(s, _)| display_formula(s)).join(",")
                                )))
                                .content(format!(
                                    "Positions: {}{}{}",
                                    render_places(right_places),
                                    if stubs.is_empty() { "" } else { ", stubs: " },
                                    stubs.iter().map(|(_, s)| display_formula(s)).join(",")
                                ))
                                .content(format!(
                                    "{}{}",
                                    if diagnostic_ions.is_empty() {
                                        ""
                                    } else {
                                        "<br>Diagnostic ions: "
                                    },
                                    diagnostic_ions
                                        .iter()
                                        .map(|d| &d.0)
                                        .map(display_formula)
                                        .join(",")
                                ));
                        }
                    }
                }
                output = output.content(ul);
            }

            Ok(output.to_string())
        }
        ModificationSearchResult::Mass(_, _, modifications) => {
            Ok(html_builder::HtmlElement::table(
                Some(&["Name", "Id", "Monoisotopic mass", "Formula"]),
                modifications
                    .iter()
                    .map(|(ontology, id, name, modification)| {
                        [
                            modification.to_string(),
                            link_modification(*ontology, *id, name),
                            display_mass(modification.formula().monoisotopic_mass()).to_string(),
                            display_formula(&modification.formula()),
                        ]
                    }),
            )
            .to_string())
        }
        ModificationSearchResult::Formula(_, modifications) => {
            Ok(html_builder::HtmlElement::table(
                Some(&["Name", "Id"]),
                modifications
                    .iter()
                    .map(|(ontology, id, name, modification)| {
                        [
                            modification.to_string(),
                            link_modification(*ontology, *id, name),
                        ]
                    }),
            )
            .to_string())
        }
        ModificationSearchResult::Glycan(_, modifications) => Ok(html_builder::HtmlElement::table(
            Some(&["Name", "Structure"]),
            modifications
                .iter()
                .map(|(ontology, id, name, modification)| {
                    [
                        link_modification(*ontology, *id, name),
                        if let SimpleModification::Gno(GnoComposition::Structure(structure), _) =
                            modification
                        {
                            structure.to_string()
                        } else {
                            "-".to_string()
                        },
                    ]
                }),
        )
        .to_string()),
    }
}

fn render_modification_id(id: &ModificationId) -> HtmlElement {
    let mut text = HtmlTag::div
        .empty()
        .content(HtmlTag::p.empty().content(format!(
            "Ontology: {}, name: {}, index: {}, ",
            id.ontology, id.name, id.id
        )))
        .maybe_content((!id.description.is_empty()).then_some(&id.description));
    if !id.synonyms.is_empty() {
        text = text.content(HtmlTag::p.empty().content("Synonyms"));
        text = text.content(
            HtmlTag::ul.empty().children(
                id.synonyms
                    .iter()
                    .map(|syn| HtmlTag::li.empty().content(syn)),
            ),
        )
    }
    if !id.cross_ids.is_empty() {
        text = text.content(HtmlTag::p.empty().content("Cross IDs"));
        text = text.content(
            HtmlTag::ul.empty().children(
                id.cross_ids
                    .iter()
                    .map(|(meta, id)| HtmlTag::li.empty().content(format!("{meta}: {id}"))),
            ),
        )
    }
    text.maybe_content(id.url().map(|url| {
        HtmlElement::new(HtmlTag::a)
            .content("ontology link")
            .header("href", url)
            .header("target", "_blank")
    }))
}
