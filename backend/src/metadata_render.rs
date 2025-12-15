use itertools::Itertools;
use mzannotate::annotation::model::BuiltInFragmentationModel;
use mzdata::prelude::SpectrumLike;
use mzident::{
    CVTerm, MSFraggerOpenModification, PSM, PSMData, PSMMetaData, ProteinMetaData, SpectrumIds
};

use crate::{
    Theme,
    html_builder::{HtmlContent, HtmlElement, HtmlTag},
    render::{display_mass, display_masses, render_peptide},
    spectra::spectrum_description,
};

pub trait RenderToHtml {
    fn to_html(&self, theme: Theme) -> HtmlElement;
}

pub trait RenderToTable {
    fn to_table(&self, theme: Theme) -> Vec<(&'static str, String)>;
}

impl<C, A> RenderToHtml for PSM<C, A> {
    fn to_html(&self, theme: Theme) -> HtmlElement {
        // Render the peptide with its local confidence
        let peptide = if let Some(peptide) = self.compound_peptidoform_ion() {
            let mut glycan_footnotes = Vec::new();
            let mut buffer = String::new();
            render_peptide(
                &mut buffer,
                &peptide,
                None,
                self.local_confidence
                    .as_ref()
                    .map(|lc| vec![vec![lc.clone()]]),
                Some(self.flanking_sequences()),
                &mut glycan_footnotes,
                theme,
            );
            for (index, footnote) in glycan_footnotes.into_iter().enumerate() {
                buffer.push_str(&format!(
                    "{}<span class='glycan-footnote'>{}: {footnote}</span>",
                    if index != 0 { ", " } else { "" },
                    index + 1
                ));
            }

            HtmlContent::Text(buffer)
        } else {
            HtmlContent::Html(HtmlTag::p.new().content("No peptide").clone())
        };

        let formula = self
            .compound_peptidoform_ion()
            .map(|p| p.formulas()[0].clone());
        HtmlTag::div.new()
            .children([
                HtmlTag::p.new()
                    .content(format!(
                        "<span class='colour-dot {reliability}' title='Reliability: {reliability}'></span> Score:&nbsp;{}, ALC:&nbsp;{}, Original:&nbsp;{}, Length:&nbsp;{}, Mass:&nbsp;{}, Charge:&nbsp;{}, m/z:&nbsp;{}",
                        self.score.map_or(String::from("-"), |s| format!("{s:.3}")),
                        self.local_confidence.as_ref().map_or(String::from("-"), |lc| format!("{:.3}", lc.iter().sum::<f64>() / lc.len() as f64)),
                        self.original_confidence().map_or(String::from("-"), |(s, t)| format!("<span title='{}|{}'>{s:.3}</span>", t.accession, t.name)),
                        self.compound_peptidoform_ion().and_then(|p| p.singular_peptidoform_ref().map(|p| p.len()))
                            .to_optional_string(),
                        formula
                            .as_ref()
                            .map(display_masses)
                            .to_optional_string(),
                        self.charge()
                            .map_or("-".to_string(), |c| c.value.to_string()),
                        self.charge()
                            .and_then(|c| formula.map(|f| (f, c)))
                            .map(|(f, c)| format!(
                                "{:.3}&nbsp;Th",
                                f.monoisotopic_mass().value / c.value as f64
                            ))
                            .to_optional_string(),
                        reliability = self.reliability().map_or("unknown".to_string(), |r| r.to_string()),
                    ))
                    .clone(),
                HtmlTag::p.new()
                    .content(format!(
                        "Experimental&nbsp;mass:&nbsp;{}, Experimental&nbsp;m/z:&nbsp;{}, Mass&nbsp;error:&nbsp;{}&nbsp;ppm&nbsp;/&nbsp;{}, Mode:&nbsp;{}, RT:&nbsp;{}",
                        self.experimental_mass()
                            .map(|m| display_mass(m, None))
                            .to_optional_string(),
                        self.experimental_mz().map(|mz| format!(
                            "{:.3}&nbsp;Th",
                            mz.value
                        )).to_optional_string(),
                        self.ppm_error().map(|ppm| format!(
                            "{:.2}",
                            ppm.value * 10e6
                        )).to_optional_string(),
                        self.mass_error()
                            .map(|m| display_mass(m, None))
                            .to_optional_string(),
                        self.fragmentation_model().map_or("-".to_string(), |f| match (f, self.mode()) {
                            (BuiltInFragmentationModel::All | BuiltInFragmentationModel::None, Some(m)) if !m.is_empty() => format!("{f}&nbsp;({m})"),
                            _ => f.to_string(),
                        }),
                        self.retention_time()
                            .map_or("-".to_string(), |c| format!("{:.3}&nbsp;min", c.get::<mzcore::system::time::min>())),
                    ))
                    .clone(),
                HtmlTag::p.new()
                    .content(format!(
                        "Protein&nbsp;id:&nbsp;{}, Name:&nbsp;{}, Database:&nbsp;{}, Location:&nbsp;{}, Unique:&nbsp{}",
                        self.proteins().iter().filter_map(|n| n.numerical_id()).join(";"),
                        self.proteins().iter().map(|n| n.id().name().to_string()).join(";"),
                        self.database().map(|(db, version)| format!("{db}{}", version.map(|v| format!(" ({v})")).unwrap_or_default())).to_optional_string(),
                        self.protein_location().map(|s| format!("{} — {}", s.start, s.end)).to_optional_string(),
                        self.unique().to_optional_string(),
                    ))
                    .clone(),
                HtmlTag::ul.new().children(match self.scans() {
                    SpectrumIds::None => vec![HtmlTag::p.new().content("No spectrum reference").clone()],
                    SpectrumIds::FileNotKnown(scans) => vec![HtmlTag::li.new().content("Scans: ").content(scans.iter().join(";")).clone()],
                    SpectrumIds::FileKnown(scans) => scans.iter().map(|(file, scans)| HtmlTag::li.new().content("File: ").content(HtmlTag::span.new().title(file.to_string_lossy()).content(file.file_name().map_or(String::new(), |s| s.to_string_lossy().to_string())).content(" Scans: ").content(scans.iter().join(";"))).clone()).collect(),
                }).clone(),
                HtmlTag::p.new().maybe_content(self.annotated_spectrum().map(|a| spectrum_description(
                                &a.description,
                                &a.peaks().fetch_summaries(),
                            ))).clone(),
                ]
            ).content(peptide).children([
                HtmlTag::p.new().content(format!("Additional MetaData <span title='{}'>{}</span> ID: {}", 
                    self.search_engine().map(|s| format!("{}|{}", s.accession, s.name)).unwrap_or_default(), 
                    self.format(), 
                    self.id())).clone(),
                {
                    HtmlElement::table::<HtmlContent, String>(
                        None,
                        self.data.to_table(theme).into_iter().map(|(h,v)| [h.to_string(), v]),
                    )
                }
            ])
            .clone()
    }
}

impl RenderToTable for PSMData {
    fn to_table(&self, theme: Theme) -> Vec<(&'static str, String)> {
        match self {
            PSMData::Peaks(data) => vec![
                ("Fraction", data.fraction.to_optional_string()),
                (
                    "Feature",
                    data.feature_tryp_cid
                        .and_then(|f_cid| data.feature_tryp_ead.map(|f_ead| (f_cid, f_ead)))
                        .map_or(
                            data.feature.as_ref().to_optional_string(),
                            |(f_cid, f_ead)| {
                                data.feature
                                    .as_ref()
                                    .map(|f| format!("{f} Tryp CID: {f_cid} Tryp EAD: {f_ead}"))
                                    .to_optional_string()
                            },
                        ),
                ),
                ("Scores", {
                    let mut output = String::new();
                    if let Some(s) = data.de_novo_score {
                        output += &format!("de novo: {s}");
                    }
                    if let Some(s) = data.alc {
                        if !output.is_empty() {
                            output.push(' ');
                        }
                        output += &format!("alc: {s}");
                    }
                    if let Some(s) = data.logp {
                        if !output.is_empty() {
                            output.push(' ');
                        }
                        output += &format!("-10logP: {s}");
                    }
                    if let Some(s) = data.quality {
                        if !output.is_empty() {
                            output.push(' ');
                        }
                        output += &format!("quality: {s}");
                    }
                    if output.is_empty() {
                        output.push('-')
                    }
                    output
                }),
                (
                    "Predicted RT",
                    data.predicted_rt.map(|v| v.value).to_optional_string(),
                ),
                (
                    "RT range",
                    data.rt_begin
                        .and_then(|b| data.rt_end.map(|e| (b, e)))
                        .map(|(b, e)| format!("begin: {:.3} min end: {:.3} min", b.value, e.value))
                        .to_optional_string(),
                ),
                (
                    "Area",
                    data.area_tryp_ead.map_or(
                        data.area.map(|a| format!("{a:e}")).to_optional_string(),
                        |a_ead| {
                            data.area
                                .map(|a| format!("Tryp CID: {a:e} Tryp EAD: {a_ead:e}"))
                                .to_optional_string()
                        },
                    ),
                ),
                (
                    "Modification scores",
                    data.ascore.as_ref().to_optional_string(),
                ),
                ("Found by", data.found_by.as_ref().to_optional_string()),
                ("From Chimera", data.from_chimera.to_optional_string()),
                ("ID", data.id.to_optional_string()),
                ("Precursor ID", data.precursor_id.to_optional_string()),
                (
                    "Ion mobility",
                    data.k0_range.as_ref().map_or("-".to_string(), |range| {
                        format!("{} — {} 1/K<sub>0</sub>", range.start(), range.end())
                    }),
                ),
                (
                    "Protein",
                    data.protein_group
                        .and_then(|g| data.protein_id.map(|i| (g, i)))
                        .map(|(g, i)| format!("group: {g}, ID: {i}"))
                        .to_optional_string(),
                ),
            ],
            PSMData::Novor(data) => vec![
                ("Spectra ID", data.spectra_id.to_optional_string()),
                ("Fraction", data.fraction.to_optional_string()),
                ("Protein", data.protein.to_optional_string()),
                (
                    "Protein Origin",
                    data.protein_origin.as_ref().to_optional_string(),
                ),
                (
                    "Protein All",
                    data.protein_all.as_ref().to_optional_string(),
                ),
                (
                    "Database Sequence",
                    data.database_sequence.as_ref().to_optional_string(),
                ),
            ],
            PSMData::Opair(data) => vec![
                (
                    "Precursor scan number",
                    data.precursor_scan_number.to_string(),
                ),
                ("Glycan Mass", data.glycan_mass.value.to_string()),
                ("Accession", data.proteins()[0].accession.to_string()),
                ("Organism", data.proteins()[0].organism.to_string()),
                ("Rank", data.rank.to_string()),
                ("Matched ion counts", data.matched_ion_counts.to_string()),
                ("Kind", data.proteins()[0].kind.to_string()),
                ("Q value", data.q_value.to_string()),
                ("PEP", data.pep.to_string()),
                ("PEP Q value", data.pep_q_value.to_string()),
                ("Localisation score", data.localisation_score.to_string()),
                ("Yion score", data.yion_score.to_string()),
                (
                    "Diagnostic ion score",
                    data.diagnostic_ion_score.to_string(),
                ),
                (
                    "Plausible glycan number",
                    data.plausible_glycan_number.to_string(),
                ),
                (
                    "Total glycosylation sites",
                    data.total_glycosylation_sites.to_string(),
                ),
                (
                    "Plausible glycan composition",
                    data.plausible_glycan_composition.to_string(),
                ),
                (
                    "Plausible glycan structure",
                    data.plausible_glycan_structure.to_string(),
                ),
                (
                    "N glycan motif check passed",
                    data.n_glycan_motif.to_string(),
                ),
                ("R138/144", data.r138_144.to_string()),
                (
                    "Glycan localisation level",
                    data.glycan_localisation_level.to_string(),
                ),
                (
                    "Glycan peptide site specificity",
                    data.glycan_peptide_site_specificity.to_string(),
                ),
                (
                    "Glycan protein site specificity",
                    data.glycan_protein_site_specificity.to_string(),
                ),
                (
                    "All potential glycan localisations",
                    data.all_potential_glycan_localisations.to_string(),
                ),
                (
                    "All site specific localisation probabilities",
                    data.all_site_specific_localisation_probabilities
                        .to_string(),
                ),
            ],
            PSMData::MetaMorpheus(data) => vec![
                (
                    "Precursor scan number",
                    data.precursor_scan_number.to_string(),
                ),
                ("Accession", data.proteins()[0].protein_accession.join("|")),
                ("Organism", data.proteins()[0].organism.join("|")),
                ("Matched ion counts", data.matched_ion_counts.to_string()),
                ("Kind", data.proteins()[0].kind.iter().join("|")),
                ("Q value", data.q_value.to_string()),
                ("PEP", data.pep.to_string()),
                ("PEP Q value", data.pep_q_value.to_string()),

            ],
            PSMData::MaxQuant(data) => vec![
                ("Scan index", data.scan_index.to_optional_string()),
                (
                    "DN sequence",
                    data.dn_sequence.as_ref().to_optional_string(),
                ),
                (
                    "DN combined score",
                    data.dn_combined_score.to_optional_string(),
                ),
                (
                    "Mass analyser",
                    data.mass_analyser.as_ref().to_optional_string(),
                ),
                ("Type", data.ty.to_string()),
                (
                    "Scan event number",
                    data.scan_event_number.to_optional_string(),
                ),
                ("Pep", data.pep.to_string()),
                (
                    "Precursor full scan number",
                    data.precursor.to_optional_string(),
                ),
                (
                    "Precursor intensity",
                    data.precursor_intensity.to_optional_string(),
                ),
                (
                    "Precursor apex",
                    data.precursor_apex_function
                        .map(|function| {
                            format!(
                                "function: {function} offset: {} offset time: {}",
                                data.precursor_apex_offset.to_optional_string(),
                                data.precursor_apex_offset_time.to_optional_string(),
                            )
                        })
                        .to_optional_string(),
                ),
                (
                    "Missed cleavages",
                    data.missed_cleavages.to_optional_string(),
                ),
                ("Isotope index", data.isotope_index.to_optional_string()),
                (
                    "Simple mass error [ppm]",
                    data.simple_mass_error_ppm.to_optional_string(),
                ),
                (
                    "Number of matches",
                    data.number_of_matches.to_optional_string(),
                ),
                (
                    "Intensity coverage",
                    data.intensity_coverage.to_optional_string(),
                ),
                ("Peak coverage", data.peak_coverage.to_optional_string()),
                ("Delta score", data.delta_score.to_optional_string()),
                ("Score diff", data.score_diff.to_optional_string()),
                (
                    "Localization probability",
                    data.localisation_probability.to_optional_string(),
                ),
                (
                    "All modified sequences",
                    data.all_modified_sequences
                        .as_ref()
                        .map(|v| v.iter().map(ToString::to_string).join(","))
                        .to_optional_string(),
                ),
                (
                    "Protein group ids",
                    data.protein_group_ids
                        .as_ref()
                        .map(|v| v.iter().map(ToString::to_string).join(","))
                        .to_optional_string(),
                ),
                (
                    "IDs",
                    format!(
                        "id: {} peptide: {} mod. peptide: {} evidence: {}",
                        data.id.to_optional_string(),
                        data.peptide_id.to_optional_string(),
                        data.modified_peptide_id.to_optional_string(),
                        data.evidence_id.to_optional_string()
                    ),
                ),
                (
                    "Base peak intensity",
                    data.base_peak_intensity.to_optional_string(),
                ),
                (
                    "Total ion current",
                    data.total_ion_current.to_optional_string(),
                ),
                (
                    "Collision energy",
                    data.collision_energy.to_optional_string(),
                ),
                (
                    "DN mass left",
                    data.dn_n_mass
                        .map(|n| {
                            format!(
                                "N: {} C: {}",
                                display_mass(n, None),
                                data.dn_c_mass
                                    .map(|v| display_mass(v, None),)
                                    .to_optional_string(),
                            )
                        })
                        .to_optional_string(),
                ),
                (
                    "Ratio H/L",
                    data.ration_h_l
                        .map(|v| {
                            format!(
                                "{v} normalised: {}",
                                data.ration_h_l_normalised.to_optional_string()
                            )
                        })
                        .to_optional_string(),
                ),
                (
                    "Intensity",
                    data.intensity
                        .map(|i| {
                            format!(
                                "{i} H: {} L: {}",
                                data.intensity_h.to_optional_string(),
                                data.intensity_l.to_optional_string()
                            )
                        })
                        .to_optional_string(),
                ),
                ("Experiment", data.experiment.as_ref().to_optional_string()),
                ("Labeling state", data.labeling_state.to_optional_string()),
            ],
            PSMData::Sage(data) => vec![
                ("Rank", data.rank.to_string()),
                ("Decoy", data.decoy.to_string()),
                ("Missed cleavages", data.missed_cleavages.to_string()),
                ("Semi enzymatic", data.semi_enzymatic.to_string()),
                ("Isotope error", data.isotope_error.to_string()),
                (
                    "Fragment error (average ppm)",
                    data.fragment_ppm.get::<mzcore::system::ratio::ppm>().to_string(),
                ),
                ("Hyperscore", format!("{:.3} Δnext {:.3} Δbest {:.3}", data.hyperscore, data.delta_next, data.delta_best)),
                ("Retention time (fraction)", format!("{:.3} predicted {:.3} Δ {:.3}", 
                    data.aligned_rt.get::<mzcore::system::ratio::fraction>(), 
                    data.predicted_rt.get::<mzcore::system::ratio::fraction>(), 
                    data.delta_rt_model)),
                ("Ion mobility", format!("{:.3} predicted {:.3} Δ {:.3}", 
                    data.ion_mobility, 
                    data.predicted_mobility, 
                    data.delta_mobility)),
                ("Matched peaks", data.matched_peaks.to_string()),
                ("Longest b", data.longest_b.to_string()),
                ("Longest y", data.longest_y.to_string()),
                (
                    "Matched intensity",
                    data.matched_intensity_pct.value.to_string(),
                ),
                ("Scored candidates", data.scored_candidates.to_string()),
                ("Poisson", data.poisson.to_string()),
                (
                    "Sage discriminant score",
                    data.sage_discriminant_score.to_string(),
                ),
                ("Posterior error", data.posterior_error.to_string()),
                ("Spectrum q", data.spectrum_q.to_string()),
                ("Peptide q", data.peptide_q.to_string()),
                ("Protein q", data.protein_q.to_string()),
                ("MS2 intensity", format!("{:e}", data.ms2_intensity)),
            ],
            PSMData::MSFragger(data) => vec![
                (
                    "Open modification",
                    match &data.open_search_modification {
                        Some(MSFraggerOpenModification::Glycan(composition)) => composition
                            .iter()
                            .fold(String::new(), |acc, m| acc + &format!("{}{}", m.0, m.1)),
                        Some(MSFraggerOpenModification::Other(other)) => other.to_string(),
                        None => "-".to_string(),
                    },
                ),
                (
                    "Calibrated experimental mass",
                    data.calibrated_experimental_mass
                        .map(|m| display_mass(m, None))
                        .to_optional_string(),
                ),
                (
                    "Calibrated experimental m/z",
                    data.calibrated_experimental_mz
                        .map(|m| m.value)
                        .to_optional_string(),
                ),
                ("Expectation score", format!("{:e}", data.expectation_score)),
                ("Next score", data.next_score.to_string()),
                (
                    "Peptide Prophet probability",
                    data.peptide_prophet_probability.to_optional_string(),
                ),
                ("Missed cleavages", data.missed_cleavages.to_string()),
                (
                    "Number of enzymatic termini",
                    data.enzymatic_termini.to_optional_string(),
                ),
                (
                    "Intensity",
                    data.intensity
                        .map(|i| format!("{i:e}"))
                        .to_optional_string(),
                ),
                (
                    "Glycan",
                    format!(
                        "composition: {}, score: {}, q-value: {}",
                        data.total_glycan_composition
                            .as_ref()
                            .or(match data.open_search_modification.as_ref() {
                                Some(MSFraggerOpenModification::Glycan(composition)) =>
                                    Some(composition),
                                _ => None,
                            })
                            .map(|composition| composition
                                .iter()
                                .fold(String::new(), |acc, m| acc + &format!("{}{}", m.0, m.1)))
                            .to_optional_string(),
                        data.glycan_score.to_optional_string(),
                        data.glycan_q_value.to_optional_string(),
                    ),
                ),
                (
                    "Modification position scores",
                    data.open_search_position_scores
                        .as_ref()
                        .map(|v| v.iter().join(","))
                        .to_optional_string(),
                ),
                ("Purity", data.purity.to_optional_string()),
                ("Gene", data.proteins()[0].gene.as_ref().to_optional_string()),
                (
                    "Protein description",
                    data.proteins()[0].protein_description.as_ref().to_optional_string(),
                ),
                (
                    "Mapped genes",
                    data.proteins()[0].mapped_genes
                        .as_ref()
                        .map(|g| g.join(","))
                        .to_optional_string(),
                ),
                (
                    "Mapped proteins",
                    data.proteins()[0].mapped_proteins
                        .as_ref()
                        .map(|p| p.join(","))
                        .to_optional_string(),
                ),
                ("Condition", data.condition.as_ref().to_optional_string()),
                ("Group", data.group.as_ref().to_optional_string()),
            ],
            PSMData::MzTab(data) => vec![
                (
                    "Accession",
                    data.protein
                        .as_ref()
                        .map(|(acc, prot)| prot.as_ref().map(|p| &p.accession).unwrap_or(acc))
                        .to_optional_string(),
                ),
                (
                    "Search engine",
                    data.search_engine
                        .iter()
                        .map(|(engine, score, score_type)| {
                            format!(
                                "{} score: {} ({})",
                                engine.to_html(theme),
                                score.to_optional_string(),
                                score_type.to_html(theme),
                            )
                        })
                        .join("|"),
                ),
                ("Uri", data.uri.as_ref().to_optional_string()),
                (
                    "Other columns",
                    (!data.additional.is_empty())
                        .then(|| {
                            HtmlElement::table(
                                Some(&["Column", "Value"]),
                                data.additional.iter().map(|(k, v)| [k, v]),
                            )
                        })
                        .to_optional_string(),
                ),
            ],
            PSMData::PLink(data) => vec![
                ("Peptide type", data.peptide_type.to_string()),
                ("Score", format!("{:e}", data.score)),
                ("Refined score", data.refined_score.to_string()),
                ("SVM score", data.svm_score.to_string()),
                ("E-value", data.e_value.to_string()),
                ("Q-value", data.q_value.to_string()),
                ("Decoy", data.is_decoy.to_string()),
                (
                    "Proteins",
                    data.proteins
                        .iter()
                        .map(|(p1, pos1, p2, pos2)| {
                            format!(
                                "{p1}{}{}{}",
                                pos1.map_or(String::new(), |p| format!("({p})")),
                                p2.as_ref().map_or(String::new(), |p| format!("-{p}")),
                                pos2.map_or(String::new(), |p| format!("({p})"))
                            )
                        })
                        .join("/"),
                ),
                (
                    "Between different proteins",
                    data.is_different_protein.to_string(),
                ),
                ("Raw file ID", data.raw_file_id.to_string()),
                ("Complex satisfied", data.is_complex_satisfied.to_string()),
                ("In filter", data.is_filter_in.to_string()),
                ("Title", data.title.clone()),
            ],
            PSMData::PUniFind(data) => {
                vec![("Cosine similarity", data.cos_similarity.to_string())]
            }
            PSMData::Fasta(data) => vec![
                ("Description", data.description().to_string()),
                (
                    "Tags",
                    (data.tags().next().is_some())
                        .then(|| {
                            HtmlElement::table(
                                Some(&["Tag", "Value"]),
                                data.tags().map(|(k, v)| [k, v]),
                            )
                        })
                        .to_optional_string(),
                ),
                ("Full header", data.header().to_string()),
            ],
            PSMData::Proteoscape(data) => vec![
                (
                    "O/k0",
                    format!(
                        "Corrected: {:.3}, predicted: {:.3}",
                        data.corrected_o_over_k0, data.predicted_o_over_k0
                    ),
                ),
                (
                    "Scores",
                    format!(
                        "Xcorr: {:.3}, Δ CN: {:.3}, TIM: {:.3}",
                        data.xcorr_score, data.delta_cn_score, data.tim_score
                    ),
                ),
                ("Matched ions", data.matched_ions.to_string()),
            ],
            PSMData::NovoB(data) => vec![
                (
                    "Score",
                    format!(
                        "forward: {:.3} reverse: {:.3}",
                        data.score_forward, data.score_reverse
                    ),
                ),
                (
                    "PPM dif",
                    format!(
                        "forward: {:.3} reverse: {:.3}",
                        data.ppm_diff_forward.value * 1e6,
                        data.ppm_diff_reverse.value * 1e6
                    ),
                ),
                (
                    "Sequence forward",
                    data.peptide_forward.as_ref().to_optional_string(),
                ),
                (
                    "Sequence reverse",
                    data.peptide_reverse.as_ref().to_optional_string(),
                ),
            ],
            PSMData::PLGS(data) => vec![
                (
                    "Protein score",
                    format!(
                        "<span class='colour-dot {curate}' title='{curate}'></span> score: {} fpr: {} coverage: {}%",
                        data.proteins()[0].protein_score,
                        data.proteins()[0].protein_fpr,
                        data.proteins()[0].protein_sequence_coverage,
                        curate = data.proteins()[0].protein_auto_curate,
                    ),
                ),
                (
                    "Protein matched peptides",
                    format!(
                        "{} (∑I: {:.3e}) (top3: {:.3e})",
                        data.proteins()[0].protein_matched_peptides,
                        data.proteins()[0].protein_matched_peptide_intensity_sum,
                        data.proteins()[0].protein_matched_peptide_intensity_top3,
                    ),
                ),
                (
                    "Protein matched products",
                    format!(
                        "{} (∑I: {:.3e})",
                        data.proteins()[0].protein_matched_products, data.proteins()[0].protein_matched_product_intensity_sum,
                    ),
                ),
                (
                    "Peptide score",
                    format!(
                        "raw: {:.3} rank: {}",
                        data.peptide_raw_score,
                        data.peptide_rank,
                    ),
                ),
                (
                    "Peptide metrics",
                    format!(
                        "pI: {} volume: {} csa: {:.3} x-P bond: {} charge: {:+} ({:+.2})",
                        data.peptide_pi,
                        data.peptide_volume,
                        data.peptide_csa,
                        data.peptide_x_p_bond_identified.to_optional_string(),
                        data.precursor_z.value,
                        data.precursor_charge,
                    ),
                ),
                (
                    "Peptide matched products",
                    format!(
                        "{:.2}% ({}/{}) (∑I: {:.3e})",
                        data.peptide_matched_products as f64
                            / data.peptide_matched_product_theoretical as f64
                            * 100.0,
                        data.peptide_matched_products,
                        data.peptide_matched_product_theoretical,
                        data.peptide_matched_product_intensity,
                    ),
                ),
                (
                    "Peptide matched products",
                    data.peptide_matched_product_string.to_string(),
                ),
                ("Precursor id", data.precursor_le_id.to_string()),
                (
                    "Precursor peak fwhm",
                    format!(
                        "{:.3} (rms: {:.3})",
                        data.precursor_fwhm, data.precursor_rms_fwhm_delta,
                    ),
                ),
                (
                    "Precursor RT (min)",
                    format!(
                        "{:.3}<span style='color:var(--color-red)'>◞</span> {:.3} <span style='color:var(--color-red)'>◠</span> {:.3} <span style='color:var(--color-red)'>◟</span>{:.3}",
                        data.precursor_lift_off_rt.value,
                        data.precursor_inf_up_rt.value,
                        data.precursor_inf_down_rt.value,
                        data.precursor_touch_down_rt.value,
                    ),
                ),
                ("Product", {
                    let charge = data
                        .product_z
                        .map_or(String::new(), |c| format!("+{}", c.value));
                    format!(
                            "{}{}{}{} id: {} intensity: {} charge: {}",
                            data.fragment_type.as_ref().to_optional_string(),
                            data.product_charge
                                .map_or(String::new(), |_| format!("<sup>{charge}</sup>")),
                            data.fragment_index
                                .map_or(String::new(), |i| format!("<sub style='margin-left:-{w}ch;min-width:{w}ch;display:inline-block'>{i}</sub>", w=charge.len())),
                            data.fragment_neutral_loss
                                .as_ref()
                                .map_or(String::new(), |n| n.hill_notation_html()),
                            data.product_he_id.to_optional_string(),
                            data.product_intensity
                                .map(|v| format!("{v:.3e}"))
                                .to_optional_string(),
                            data.product_charge
                                .map(|v| format!("{v:+.2}"))
                                .to_optional_string(),
                        )
                }),
                ("Product peak fwhm", data.product_fwhm.to_optional_string()),
                (
                    "Product RT (min)",
                    format!(
                        "{:.3}<span style='color:var(--color-red)'>◞</span> {} <span style='color:var(--color-red)'>◠</span> {:.3} <span style='color:var(--color-red)'>◟</span>{:.3}",
                        data.product_lift_off_rt
                            .map(|v| format!("{:.3}", v.value))
                            .to_optional_string(),
                        data.product_inf_up_rt
                            .map(|v| format!("{:.3}", v.value))
                            .to_optional_string(),
                        data.product_inf_down_rt
                            .map(|v| format!("{:.3}", v.value))
                            .to_optional_string(),
                        data.product_touch_down_rt
                            .map(|v| format!("{:.3}", v.value))
                            .to_optional_string(),
                    ),
                ),
            ],
            PSMData::SpectrumSequenceList(data) => vec![
                (
                    "Retention range (min)",
                    data.start_time.map_or("-".to_string(), |st| {
                        format!(
                            "{:.3} — {}",
                            st.get::<mzcore::system::time::min>(),
                            data.end_time.map_or("-".to_string(), |c| format!(
                                "{:.3}",
                                c.get::<mzcore::system::time::min>()
                            ))
                        )
                    }),
                ),
                ("Score type", data.score_type.as_ref().to_optional_string()),
                ("Adduct", data.adduct.as_ref().to_optional_string()),
                (
                    "Molecule name",
                    data.moleculename.as_ref().to_optional_string(),
                ),
                ("Inchikey", data.inchikey.as_ref().to_optional_string()),
                ("Other keys", data.otherkeys.as_ref().to_optional_string()),
                (
                    "Ion Mobility",
                    data.ion_mobility.as_ref().map_or("-".to_string(), |im| {
                        format!(
                            "{im:.3} {}",
                            data.ion_mobility_units.as_ref().to_optional_string()
                        )
                    }),
                ),
                (
                    "CCS",
                    data.ccs
                        .map_or("-".to_string(), |ccs| format!("{ccs:.3} Å²")),
                ),
            ],
            PSMData::AnnotatedSpectrum(_) // TODO: make a nice display
            | PSMData::DeepNovoFamily(_)
            | PSMData::InstaNovo(_)
            | PSMData::PowerNovo(_)
            | PSMData::PepNet(_)
            | PSMData::BasicCSV(_)
            | PSMData::PiHelixNovo(_)
            | PSMData::PiPrimeNovo(_) => Vec::new(),
        }
    }
}

pub trait OptionalString {
    fn to_optional_string(self) -> String;
}

impl<T: ToString> OptionalString for Option<T> {
    fn to_optional_string(self) -> String {
        self.map_or("-".to_string(), |v| v.to_string())
    }
}

impl RenderToHtml for CVTerm {
    fn to_html(&self, _theme: Theme) -> HtmlElement {
        HtmlTag::span
            .new()
            .content(self.term.accession.to_string())
            .title(format!(
                "[{}, {}, {}, {}]",
                self.term.accession.cv, self.term.accession.accession, self.term.name, self.comment
            ))
            .clone()
    }
}
