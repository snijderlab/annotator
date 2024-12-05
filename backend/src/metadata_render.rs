use itertools::Itertools;
use rustyms::{
    identification::{CVTerm, IdentifiedPeptide, MetaData, ReturnedPeptide, SpectrumIds},
    MultiChemical,
};

use crate::{
    html_builder::{HtmlContent, HtmlElement, HtmlTag},
    render::{display_mass, display_masses, render_peptide},
};

pub trait RenderToHtml {
    fn to_html(&self) -> HtmlElement;
}

pub trait RenderToTable {
    fn to_table(&self) -> Vec<(&'static str, String)>;
}

impl RenderToHtml for IdentifiedPeptide {
    fn to_html(&self) -> HtmlElement {
        // Render the peptide with its local confidence
        let peptide = if let Some(peptide) = self.peptide() {
            let mut buffer = String::new();
            render_peptide(
                &mut buffer,
                &peptide.compound_peptidoform(),
                None,
                self.local_confidence
                    .as_ref()
                    .map(|lc| vec![vec![lc.clone()]]),
            );

            HtmlContent::Text(buffer)
        } else {
            HtmlContent::Html(HtmlTag::p.new().content("No peptide").clone())
        };

        let formula = self.peptide().map(|p| p.formulas()[0].clone());
        HtmlTag::div.new()
            .children([
                HtmlTag::p.new()
                    .content(format!(
                        "Score:&nbsp;{}, ALC:&nbsp;{}, Length:&nbsp;{}, Mass:&nbsp;{}, Charge:&nbsp;{}, m/z:&nbsp;{}, Mode:&nbsp;{}, RT:&nbsp;{}",
                        self.score.map_or(String::from("-"), |s| format!("{s:.3}")),
                        self.local_confidence.as_ref().map_or(String::from("-"), |lc| format!("{:.3}", lc.iter().sum::<f64>() / lc.len() as f64)),
                        self.peptide()
                            .map(|p| format!("{}&nbsp;AA", match p {
                                ReturnedPeptide::LinearSemiAmbiguous(p) => p.len().to_string(),
                                ReturnedPeptide::LinearSimpleLinear(p) => p.len().to_string(),
                                ReturnedPeptide::Peptidoform(p) => p.peptides().iter().map(|p| p.len()).join("&nbsp;+&nbsp;"),
                                ReturnedPeptide::CompoundPeptidoform(p) => p.peptides().map(|p| p.len()).join("&nbsp;+&nbsp;"),
                            }))
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
                        self.mode()
                            .map_or("-".to_string(), |c| c.to_string()),
                        self.retention_time()
                            .map_or("-".to_string(), |c| format!("{:.3}&nbsp;min", c.get::<rustyms::system::time::min>())),
                    ))
                    .clone(),
                HtmlTag::p.new()
                    .content(format!(
                        "Experimental&nbsp;mass:&nbsp;{}, Experimental&nbsp;m/z:&nbsp;{}, Mass&nbsp;error&nbsp;[ppm]:&nbsp;{}, Mass&nbsp;error&nbsp;[Abs]:&nbsp;{}",
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
                    ))
                    .clone(),
                HtmlTag::ul.new().children(match self.scans() {
                    SpectrumIds::None => vec![HtmlTag::p.new().content("No spectrum reference").clone()],
                    SpectrumIds::FileNotKnown(scans) => vec![HtmlTag::li.new().content("Scans: ").content(scans.iter().join(";")).clone()],
                    SpectrumIds::FileKnown(scans) => scans.iter().map(|(file, scans)| HtmlTag::li.new().content("File: ").content(HtmlTag::span.new().title(file.to_string_lossy()).content(file.file_name().map_or(String::new(), |s| s.to_string_lossy().to_string())).content(" Scans: ").content(scans.iter().join(";"))).clone()).collect(),
                }).clone(),]
            ).content(peptide).children([
                HtmlTag::p.new().content(format!("Additional MetaData {} {} ID: {}", self.format_name(), self.format_version(), self.id())).clone(),
                {
                    HtmlElement::table::<HtmlContent, String>(
                        None,
                        self.metadata.to_table().into_iter().map(|(h,v)| [h.to_string(), v]),
                    )
                }
            ])
            .clone()
    }
}

impl RenderToTable for MetaData {
    fn to_table(&self) -> Vec<(&'static str, String)> {
        match self {
            MetaData::Peaks(data) => vec![
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
                ("Accession", data.accession.as_ref().to_optional_string()),
                ("From Chimera", data.from_chimera.to_optional_string()),
                ("ID", data.id.to_optional_string()),
                (
                    "Protein",
                    data.protein_accession
                        .as_ref()
                        .and_then(|a| {
                            data.protein_group.and_then(|g| {
                                data.protein_id
                                    .and_then(|i| data.unique.map(|u| (a, g, i, u)))
                            })
                        })
                        .map(|(a, g, i, u)| {
                            format!("accession: {a}, group: {g}, ID: {i}, unique: {u}")
                        })
                        .to_optional_string(),
                ),
                (
                    "Protein location",
                    data.start
                        .and_then(|s| data.end.map(|e| (s, e)))
                        .map(|(s, e)| format!("{s} — {e}"))
                        .to_optional_string(),
                ),
                (
                    "Flanking residues",
                    data.peptide
                        .0
                        .or(data.peptide.2)
                        .map(|_| {
                            format!(
                                "{}_(seq)_{}",
                                data.peptide
                                    .0
                                    .as_ref()
                                    .map_or("Terminal".to_string(), |p| p.to_string()),
                                data.peptide
                                    .2
                                    .as_ref()
                                    .map_or("Terminal".to_string(), |p| p.to_string())
                            )
                        })
                        .to_optional_string(),
                ),
            ],
            MetaData::Novor(data) => vec![
                ("Score", data.score.to_string()),
                ("ID", data.id.to_optional_string()),
                ("Spectra ID", data.spectra_id.to_optional_string()),
                ("Fraction", data.fraction.to_optional_string()),
                ("Protein", data.protein.to_optional_string()),
                ("Protein location", data.protein_start.map(|s| format!("{s} — {}", s + data.peptide.len())).to_optional_string()),
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
            MetaData::Opair(data) => vec![
                (
                    "Precursor scan number",
                    data.precursor_scan_number.to_string(),
                ),
                ("Glycan Mass", data.glycan_mass.value.to_string()),
                ("Accession", data.accession.to_string()),
                ("Organism", data.organism.to_string()),
                ("Protein name", data.protein_name.to_string()),
                (
                    "Protein location",
                    format!("{} — {}", data.protein_location.0, data.protein_location.1),
                ),
                (
                    "Flanking residues",
                    format!(
                        "{}_(seq)_{}",
                        data.flanking_residues.0.char(),
                        data.flanking_residues.1.char()
                    ),
                ),
                ("Rank", data.rank.to_string()),
                ("Matched ion series", data.matched_ion_series.to_string()),
                (
                    "Matched ion mz ratios",
                    data.matched_ion_mz_ratios.to_string(),
                ),
                (
                    "Matched ion mass error",
                    data.matched_ion_mass_error.to_string(),
                ),
                (
                    "Matched ion mass error (ppm)",
                    data.matched_ion_ppm.to_string(),
                ),
                (
                    "Matched ion intensities",
                    data.matched_ion_intensities.to_string(),
                ),
                ("Matched ion counts", data.matched_ion_counts.to_string()),
                ("Kind", data.kind.to_string()),
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
            MetaData::MaxQuant(data) => vec![
                ("Scan index", data.scan_index.to_optional_string()),
                ("Modifications", data.modifications.to_string()),
                ("Proteins", data.proteins.to_string()),
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
                ("Score", data.score.to_string()),
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
                    "DN sequence",
                    data.dn_sequence.as_ref().to_optional_string(),
                ),
                (
                    "DN combined score",
                    data.dn_combined_score.to_optional_string(),
                ),
                (
                    "DN mass left",
                    data.dn_missing_mass
                        .map(|v| {
                            format!(
                                "N: {} other: {} C: {}",
                                data.dn_n_mass
                                    .map(|v| display_mass(v, None),)
                                    .to_optional_string(),
                                display_mass(v, None),
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
            MetaData::Sage(data) => vec![
                ("Proteins", data.proteins.join(";")),
                ("Rank", data.rank.to_string()),
                ("Decoy", data.decoy.to_string()),
                ("Missed cleavages", data.missed_cleavages.to_string()),
                ("Semi enzymatic", data.semi_enzymatic.to_string()),
                ("Isotope error", data.isotope_error.to_string()),
                (
                    "Fragment error (average ppm)",
                    data.fragment_ppm.value.to_string(),
                ),
                ("Hyperscore", data.hyperscore.to_string()),
                ("Hyperscore delta next", data.delta_next.to_string()),
                ("Hyperscore delta best", data.delta_best.to_string()),
                ("Retention time aligned", data.aligned_rt.value.to_string()),
                (
                    "Retention time predicted",
                    data.predicted_rt.value.to_string(),
                ),
                ("Retention time delta", data.delta_rt_model.to_string()),
                ("Ion mobility", data.ion_mobility.to_string()),
                (
                    "Ion mobility predicted",
                    data.predicted_mobility.to_string(),
                ),
                ("Ion mobility delta", data.delta_mobility.to_string()),
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
                ("MS2 intensity", data.ms2_intensity.to_string()),
            ],
            MetaData::MSFragger(data) => vec![
                (
                    "Extended Peptide",
                    format!(
                        "{}_(seq)_{}",
                        data.extended_peptide[0]
                            .as_ref()
                            .map_or("Terminal".to_string(), |p| p.to_string()),
                        data.extended_peptide[2]
                            .as_ref()
                            .map_or("Terminal".to_string(), |p| p.to_string())
                    ),
                ),
                (
                    "Assigned modifications",
                    data.assigned_modifications.to_string(),
                ),
                (
                    "Calibrated experimental mass",
                    display_mass(data.calibrated_experimental_mass, None).to_string(),
                ),
                (
                    "Calibrated experimental m/z",
                    data.calibrated_experimental_mz.value.to_string(),
                ),
                ("Expectation", data.expectation.to_string()),
                ("Hyperscore", data.hyperscore.to_string()),
                ("Next score", data.next_score.to_string()),
                (
                    "Peptide Prophet probability",
                    data.peptide_prophet_probability.to_string(),
                ),
                ("Missed cleavages", data.missed_cleavages.to_string()),
                (
                    "Number of enzymatic termini",
                    data.enzymatic_termini.to_string(),
                ),
                ("Protein location", format!("{} — {}", data.protein_start, data.protein_end)),
                ("Intensity", format!("{:e}", data.intensity)),
                ("Purity", data.purity.to_string()),
                ("Is unique", data.is_unique.to_string()),
                ("Protein", data.protein.to_string()),
                ("Protein ID", data.protein_id.to_string()),
                ("Entry name", data.entry_name.to_string()),
                ("Gene", data.gene.to_string()),
                ("Protein description", data.protein_description.to_string()),
                ("Mapped genes", data.mapped_genes.join(",")),
                ("Mapped proteins", data.mapped_proteins.join(",")),
                ("Condition", data.condition.as_ref().to_optional_string()),
                ("Group", data.group.as_ref().to_optional_string()),
            ],
            MetaData::MZTab(data) => vec![
                ("Accession", data.accession.as_ref().to_optional_string()),
                ("Unique", data.unique.to_optional_string()),
                (
                    "Database",
                    format!(
                        "db: {}, version: {}",
                        data.database.as_ref().to_optional_string(),
                        data.database_version.as_ref().to_optional_string()
                    ),
                ),
                (
                    "Search engine",
                    data.search_engine
                        .iter()
                        .map(|(engine, score, score_type)| {
                            format!(
                                "{} score: {} ({})",
                                engine.to_html(),
                                score.to_optional_string(),
                                score_type.to_html(),
                            )
                        })
                        .join("|"),
                ),
                (
                    "Reliability",
                    data.reliability.as_ref().to_optional_string(),
                ),
                ("Uri", data.uri.as_ref().to_optional_string()),
                (
                    "Protein location",
                    format!(
                        "{} — {}",
                        data.start.to_optional_string(),
                        data.end.to_optional_string()
                    ),
                ),
                (
                    "Flanking residues",
                    format!("{}_(seq)_{}", data.preceding_aa, data.following_aa),
                ),
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
            MetaData::PLink(data) => vec![
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
                                pos1.map_or(String::new(), |p| format!("({})", p)),
                                p2.as_ref().map_or(String::new(), |p| format!("-{}", p)),
                                pos2.map_or(String::new(), |p| format!("({})", p))
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
            MetaData::Fasta(data) => vec![
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
            MetaData::NovoB(data) => vec![
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
            MetaData::PLGS(data) => vec![
                (
                    "Protein",
                    format!(
                        "id: {} description: {}",
                        data.protein_id, data.protein_description
                    ),
                ),
                (
                    "Protein score",
                    format!(
                        "<span class='colour-dot {curate}' title='{curate}'></span> score: {} fpr: {} coverage: {}%",
                        data.protein_score,
                        data.protein_fpr,
                        data.protein_sequence_coverage,
                        curate=data.protein_auto_curate,
                    ),
                ),
                (
                    "Protein matched peptides",
                    format!(
                        "{} (∑I: {:.3e}) (top3: {:.3e})",
                        data.protein_matched_peptides,
                        data.protein_matched_peptide_intensity_sum,
                        data.protein_matched_peptide_intensity_top3,
                    ),
                ),
                (
                    "Protein matched products",
                    format!(
                        "{} (∑I: {:.3e})",
                        data.protein_matched_products,
                        data.protein_matched_product_intensity_sum,
                    ),
                ),
                ("Protein location",
                    format!(
                        "{} — {}",
                        data.peptide_start,
                        data.peptide_start + data.peptide.len(),
                    ),
                ),
                (
                    "Peptide score",
                    format!(
                        "<span class='colour-dot {curate}' title='{curate}'></span> score: {:.3} (raw: {:.3}) rank: {}",
                        data.peptide_score,
                        data.peptide_raw_score,
                        data.peptide_rank,
                        curate=data.peptide_auto_curate,
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
                        data.peptide_matched_products as f64 / data.peptide_matched_product_theoretical * 100.0,
                        data.peptide_matched_products,
                        data.peptide_matched_product_theoretical,
                        data.peptide_matched_product_intensity,
                    ),
                ),
                (
                    "Peptide matched products",
                    data.peptide_matched_product_string.clone(),
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
                (
                    "Product",
                    {
                        let charge = data.product_z.map_or(String::new(), |c| format!("+{}",c.value));
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
                    },
                ),
                ("Product peak fwhm", data.product_fwhm.to_optional_string()),
                (
                    "Product RT (min)",
                    format!(
                        "{:.3}<span style='color:var(--color-red)'>◞</span> {} <span style='color:var(--color-red)'>◠</span> {:.3} <span style='color:var(--color-red)'>◟</span>{:.3}",
                        data.product_lift_off_rt
                            .map(|v| format!("{:.3}", v.value))
                            .to_optional_string(),
                        data.product_inf_up_rt.map(|v| format!("{:.3}", v.value)).to_optional_string(),
                        data.product_inf_down_rt
                            .map(|v| format!("{:.3}", v.value))
                            .to_optional_string(),
                        data.product_touch_down_rt
                            .map(|v| format!("{:.3}", v.value))
                            .to_optional_string(),
                    ),
                ),
            ],
            MetaData::DeepNovoFamily(_)
            | MetaData::InstaNovo(_)
            | MetaData::PowerNovo(_)
            | MetaData::PepNet(_) => Vec::new(),
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
    fn to_html(&self) -> HtmlElement {
        HtmlTag::span
            .new()
            .content(self.term.clone())
            .title(format!(
                "[{}, {}, {}, {}]",
                self.ontology, self.id, self.term, self.comment
            ))
            .clone()
    }
}
