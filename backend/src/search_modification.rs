use crate::{
    ModifiableState, Theme, html_builder,
    metadata_render::OptionalString,
    render::{display_formula, display_masses, display_placement_rule, render_full_glycan},
};
use context_error::{BasicKind, BoxedError, CreateError, FullErrorContent};
use itertools::Itertools;
use mzcore::{
    ontology::{Ontologies, Ontology},
    prelude::*,
    quantities::Tolerance,
    sequence::{
        GnoComposition, LinkerLength, LinkerSpecificity, ModificationId, PlacementRule,
        ReturnModification, SimpleModification, SimpleModificationInner,
        modification_search_formula, modification_search_glycan, modification_search_mass,
    },
    system::{Mass, dalton},
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
    theme: Theme,
) -> Result<String, String> {
    let state = state.lock().await;
    let mut glycan_footnotes = Vec::new();
    let modification = if text.is_empty() {
        Err(BoxedError::small(
            BasicKind::Error,
            "Invalid modification",
            "The modification is empty",
        )
        .to_html(false))
    } else {
        SimpleModificationInner::pro_forma(
            text,
            &mut Vec::new(),
            &mut Vec::new(),
            &state.ontologies,
        )
        .map_err(|err| {
            BoxedError::small(
                BasicKind::Error,
                "Invalid modification",
                "Could not parse the modification",
            )
            .add_underlying_errors(err)
            .to_html(false)
        })
        .map(|((m, _), _)| match m {
            ReturnModification::Defined(d) => Ok(d),
            _ => Err(BoxedError::small(
                BasicKind::Error,
                "Invalid modification",
                "Can not define ambiguous modifications for the modifications parameter",
            )
            .to_html(false)),
        })
    }??;
    let tolerance = Tolerance::new_absolute(Mass::new::<dalton>(tolerance));

    match &*modification {
        SimpleModificationInner::Mass(_, m, _)
        | SimpleModificationInner::Gno {
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
                &state.ontologies,
            )
            .map(|modification| {
                [
                    modification.to_string(),
                    link_modification(modification.clone()),
                    display_masses(&modification.formula()).to_string(),
                    display_formula(&modification.formula(), true),
                ]
            }),
        )
        .to_string()),
        SimpleModificationInner::Formula(f) => Ok(html_builder::HtmlElement::table(
            Some(&["Name", "Id"]),
            modification_search_formula(f, &state.ontologies)
                .map(|modification| [modification.to_string(), link_modification(modification)]),
        )
        .to_string()),
        SimpleModificationInner::Glycan(g)
        | SimpleModificationInner::Gno {
            composition: GnoComposition::Composition(g),
            ..
        } => Ok(HtmlTag::div
            .new()
            .content(render_modification(
                modification.clone(),
                theme,
                &state.ontologies,
            ))
            .content(html_builder::HtmlElement::table(
                Some(&["Name", "Structure"]),
                modification_search_glycan(g, true, &state.ontologies).map(|modification| {
                    [
                        link_modification(modification.clone()),
                        if let SimpleModificationInner::Gno {
                            composition: GnoComposition::Topology(structure),
                            ..
                        } = &*modification
                        {
                            render_full_glycan(
                                structure,
                                false,
                                false,
                                theme,
                                &mut glycan_footnotes,
                                false,
                                false,
                                0,
                                0,
                            )
                        } else {
                            "-".to_string()
                        },
                    ]
                }),
            ))
            .children(glycan_footnotes.iter().enumerate().map(|(index, note)| {
                HtmlTag::span
                    .new()
                    .class("glycan-footnote")
                    .content(format!(
                        "{}{}: {note}",
                        if index != 0 { ", " } else { "" },
                        index + 1
                    ))
                    .clone()
            }))
            .to_string()),
        _ => Ok(render_modification(modification, theme, &state.ontologies).to_string()),
    }
}

pub fn render_modification(
    modification: SimpleModification,
    theme: Theme,
    ontologies: &Ontologies,
) -> HtmlElement {
    let mut output = HtmlTag::div.new();
    output.class("modification");

    output.content(HtmlTag::p.new().content(format!(
        "Formula {} Mass {}",
        display_formula(&modification.formula(), true),
        display_masses(&modification.formula()),
    )));

    if let SimpleModificationInner::Database {
        specificities, id, ..
    } = &*modification
    {
        output.content(render_modification_id(id, ontologies));
        output.content(HtmlTag::p.new().content("Placement rules"));
        let mut ul = HtmlTag::ul.new();

        for rule in specificities {
            ul.content(HtmlTag::li.new().content(format!(
                "Positions: {}{}{}{}{}",
                render_places(&rule.0, ontologies),
                if rule.1.is_empty() {
                    ""
                } else {
                    ", Neutral losses: "
                },
                rule.1.iter().map(|n| display_neutral_loss(n, true)).join(", "),
                if rule.2.is_empty() {
                    ""
                } else {
                    ", Diagnostic ions: "
                },
                rule.2.iter().map(|d| &d.0).map(|f| display_formula(f, true)).join(", ")
            )));
        }
        output.content(ul);
    } else if let SimpleModificationInner::Gno {
        composition,
        id,
        structure_score,
        subsumption_level,
        motif,
        taxonomy,
        glycomeatlas,
    } = &*modification.clone()
    {
        output.content(render_modification_id(id, ontologies));
        output.content(
            HtmlTag::p
                .new()
                .content(format!("Subsumption level: {subsumption_level}")),
        );
        if *structure_score != u16::MAX {
            output.content(
                HtmlTag::p
                    .new()
                    .content(format!("Structure score: {structure_score}"))
                    .clone(),
            );
        }
        output.maybe_content(motif.as_ref().map(|(name, _)| {
            HtmlTag::p
                .new()
                .content(format!("Motif: {name} {}", link_modification(modification)))
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
                let mut glycan_footnotes = Vec::new();
                output.extend(
                    [
                        HtmlTag::p
                            .new()
                            .content(format!(
                                "Composition: {}",
                                structure
                                    .composition()
                                    .iter()
                                    .fold(String::new(), |acc, m| acc + &format!("{}{}", m.0, m.1))
                            ))
                            .clone(),
                        HtmlTag::p
                            .new()
                            .content(render_full_glycan(
                                structure,
                                true,
                                false,
                                theme,
                                &mut glycan_footnotes,
                                false,
                                false,
                                0,
                                0,
                            ))
                            .clone(),
                    ]
                    .into_iter()
                    .chain(glycan_footnotes.iter().enumerate().map(|(index, note)| {
                        HtmlTag::span
                            .new()
                            .class("glycan-footnote")
                            .content(format!(
                                "{}{}: {note}",
                                if index != 0 { ", " } else { "" },
                                index + 1
                            ))
                            .clone()
                    })),
                );
            }
        }
    } else if let SimpleModificationInner::Linker {
        specificities,
        id,
        length,
        ..
    } = &*modification
    {
        output.content(render_modification_id(id, ontologies));
        match length {
            LinkerLength::Unknown => (),
            LinkerLength::Discreet(v) => {
                output.content(
                    HtmlTag::p
                        .new()
                        .content(format!("Length: {}", v.iter().map(|v| v.0).join(", "),)),
                );
            }
            LinkerLength::InclusiveRange(s, e) => {
                output.content(HtmlTag::p.new().content(format!("Length: {}—{}", s.0, e.0)));
            }
        }
        output.content(HtmlTag::p.new().content("Placement rules"));
        let mut ul = HtmlTag::ul.new();
        ul.class("placement-rules");

        for specificity in specificities {
            match specificity {
                LinkerSpecificity::Symmetric {
                    rules,
                    stubs,
                    neutral_losses,
                    diagnostic,
                } => {
                    ul.content(
                        HtmlTag::li
                            .new()
                            .content(format!("Positions: {}", render_places(rules, ontologies)))
                            .maybe_content((!stubs.is_empty()).then(|| {
                                format!(
                                    ", Breakages: {}",
                                    stubs.iter().map(|s| display_stubs(s, true)).join(", ")
                                )
                            }))
                            .maybe_content((!neutral_losses.is_empty()).then(|| {
                                format!(
                                    ", Neutral losses: {}",
                                    neutral_losses
                                        .iter()
                                        .map(|n| display_neutral_loss(n, true))
                                        .join(", ")
                                )
                            }))
                            .maybe_content((!diagnostic.is_empty()).then(|| {
                                format!(
                                    ", Diagnostic ions: {}",
                                    diagnostic
                                        .iter()
                                        .map(|d| display_formula(&d.0, true))
                                        .join(", ")
                                )
                            })),
                    );
                }
                LinkerSpecificity::Asymmetric {
                    rules: (left_places, right_places),
                    stubs,
                    neutral_losses,
                    diagnostic,
                } => {
                    ul.content(
                        HtmlTag::li
                            .new()
                            .content(format!(
                                "First positions: {}, Second positions: {}",
                                render_places(left_places, ontologies),
                                render_places(right_places, ontologies)
                            ))
                            .maybe_content((!stubs.is_empty()).then(|| {
                                format!(
                                    ", Breakages: {}",
                                    stubs.iter().map(|s| display_stubs(s, true)).join(", ")
                                )
                            }))
                            .maybe_content((!neutral_losses.is_empty()).then(|| {
                                format!(
                                    ", Neutral losses: {}",
                                    neutral_losses
                                        .iter()
                                        .map(|n| display_neutral_loss(n, true))
                                        .join(", ")
                                )
                            }))
                            .maybe_content((!diagnostic.is_empty()).then(|| {
                                format!(
                                    ", Diagnostic ions: {}",
                                    diagnostic
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

fn render_modification_id(id: &ModificationId, ontologies: &Ontologies) -> HtmlElement {
    let mut text = HtmlTag::div.new();
    text
        .content(HtmlTag::p.new().content(format!(
            "Ontology: <span class='ontology'>{}</span>, name: <span class='name'>{}</span>, index: <span class='index'>{}</span>{}",
            id.ontology, if id.ontology == Ontology::Gnome {id.name.to_ascii_uppercase()} else {id.name.to_string()}, id.id().to_optional_string(), if let Some(url) = id.url() {
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
                        .map(|(_, syn)| HtmlTag::li.new().content(syn).clone()),
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
                        .content(match meta.as_deref() {
                            Some("RESID")
                                if id.starts_with("AA")
                                    && id.len() > 2
                                    && id[2..].parse::<usize>().is_ok() =>
                            {
                                if let Ok(index) = id[2..].parse::<u32>()
                                    && let Some(modification) =
                                        ontologies.resid().get_by_index(&index)
                                {
                                    link_modification(modification)
                                } else {
                                    String::new()
                                }
                            }
                            Some("PSI-MOD") if id.parse::<usize>().is_ok() => {
                                if let Ok(index) = id.parse::<u32>()
                                    && let Some(modification) =
                                        ontologies.psimod().get_by_index(&index)
                                {
                                    link_modification(modification)
                                } else {
                                    String::new()
                                }
                            }
                            Some("Unimod") if id.parse::<usize>().is_ok() => {
                                if let Ok(index) = id.parse::<u32>()
                                    && let Some(modification) =
                                        ontologies.unimod().get_by_index(&index)
                                {
                                    link_modification(modification)
                                } else {
                                    String::new()
                                }
                            }
                            Some("ChEBI") if id.parse::<usize>().is_ok() => HtmlTag::a
                                .new()
                                .header(
                                    "href",
                                    format!(
                                        "https://www.ebi.ac.uk/chebi/chebiOntology.do?chebiId={id}"
                                    ),
                                )
                                .header("target", "_blank")
                                .content(format!(
                                    "<span class='cross-id-system'>ChEBI</span>: {id}"
                                ))
                                .to_string(),
                            Some("PubMed" | "PMID") if id.parse::<usize>().is_ok() => HtmlTag::a
                                .new()
                                .header("href", format!("https://pubmed.ncbi.nlm.nih.gov/{id}"))
                                .header("target", "_blank")
                                .content(format!(
                                    "<span class='cross-id-system'>PubMed</span>: {id}"
                                ))
                                .to_string(),
                            Some("DOI") => HtmlTag::a
                                .new()
                                .header("href", format!("https://doi.org/{id}"))
                                .header("target", "_blank")
                                .content(format!("<span class='cross-id-system'>DOI</span>: {id}"))
                                .to_string(),
                            Some("FindMod") => HtmlTag::a
                                .new()
                                .header("href", format!("https://web.expasy.org/findmod/{id}.html"))
                                .header("target", "_blank")
                                .content(format!(
                                    "<span class='cross-id-system'>FindMod</span>: {id}"
                                ))
                                .to_string(),
                            Some("PDBHET") => HtmlTag::a
                                .new()
                                .header("href", format!("https://www.rcsb.org/ligand/{id}"))
                                .header("target", "_blank")
                                .content(format!(
                                    "<span class='cross-id-system'>PDBHET</span>: {id}"
                                ))
                                .to_string(),
                            _ => {
                                if id.starts_with("http://") || id.starts_with("https://") {
                                    HtmlTag::a
                                        .new()
                                        .header("href", id.to_string())
                                        .header("target", "_blank")
                                        .content(meta.as_deref().unwrap_or("URL"))
                                        .to_string()
                                } else if let Some(meta) = meta {
                                    format!("<span class='cross-id-system'>{meta}</span>: {id}")
                                } else {
                                    id.to_string()
                                }
                            }
                        })
                        .clone()
                })),
        );
    }
    text
}

fn render_places(places: &[PlacementRule], ontologies: &Ontologies) -> String {
    if places.is_empty() {
        "Anywhere".to_string()
    } else {
        places
            .iter()
            .map(|p| display_placement_rule(p, true, ontologies))
            .join(", ")
    }
}
