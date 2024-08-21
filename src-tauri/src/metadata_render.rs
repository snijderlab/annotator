use itertools::Itertools;
use rustyms::identification::{
    FastaData, IdentifiedPeptide, MSFraggerData, MaxQuantData, MetaData, NovorData, OpairData,
    PeaksData, SageData,
};

use crate::{
    html_builder::{HtmlContent, HtmlElement, HtmlTag},
    render::{display_mass, display_masses},
};

pub trait RenderToHtml {
    fn to_html(&self) -> HtmlElement;
}

impl RenderToHtml for IdentifiedPeptide {
    fn to_html(&self) -> HtmlElement {
        // Render the peptide with its local confidence
        let mut peptide = HtmlElement::new(HtmlTag::div);
        peptide.class("original-sequence").style("--max-value:1");
        let lc = self
            .local_confidence
            .clone()
            .unwrap_or(vec![0.0; self.peptide.len()]);

        if let Some(n) = self.peptide.n_term.as_ref() {
            let mut modification = String::new();
            n.display(&mut modification, false).unwrap();
            peptide.content(
                HtmlElement::new(HtmlTag::div)
                    .style("--value:0")
                    .children([
                        HtmlElement::new(HtmlTag::p)
                            .content("⚬".to_string())
                            .clone(),
                        HtmlElement::new(HtmlTag::p)
                            .class("modification")
                            .content(modification)
                            .clone(),
                    ])
                    .clone(),
            );
        }

        for (aa, confidence) in self.peptide.sequence.iter().zip(lc) {
            peptide.content(
                HtmlElement::new(HtmlTag::div)
                    .style(format!("--value:{confidence}"))
                    .children([
                        HtmlElement::new(HtmlTag::p).content(aa.aminoacid.char().to_string()),
                        HtmlElement::new(HtmlTag::p).class("modification").content(
                            aa.modifications
                                .iter()
                                .map(|m| {
                                    let mut s = String::new();
                                    m.display(&mut s, false).unwrap();
                                    s
                                })
                                .join(","),
                        ),
                    ]),
            );
        }

        if let Some(c) = self.peptide.c_term.as_ref() {
            let mut modification = String::new();
            c.display(&mut modification, false).unwrap();
            peptide.content(
                HtmlElement::new(HtmlTag::div).style("--value:0").children([
                    HtmlElement::new(HtmlTag::p).content("⚬".to_string()),
                    HtmlElement::new(HtmlTag::p)
                        .class("modification")
                        .content(modification),
                ]),
            );
        }

        let formula = &self.peptide.formulas()[0];
        HtmlElement::new(HtmlTag::div)
            .children([
                HtmlElement::new(HtmlTag::p)
                    .content(format!(
                        "Score: {}, Length: {} AA, Mass: {}, Charge: {}, m/z: {} Th, Mode: {}",
                        self.score.map_or(String::from("-"), |s| format!("{s:.3}")),
                        self.peptide.sequence.len(),
                        display_masses(formula),
                        self.metadata
                            .charge()
                            .map_or("-".to_string(), |c| c.value.to_string()),
                        self.metadata.charge().map_or("-".to_string(), |c| format!(
                            "{:.3}",
                            formula.monoisotopic_mass().value / c.value as f64
                        )),
                        self.metadata
                            .mode()
                            .map_or("-".to_string(), |c| c.to_string())
                    ))
                    .clone(),
                peptide,
                self.metadata.to_html(),
            ])
            .clone()
    }
}

impl RenderToHtml for MetaData {
    fn to_html(&self) -> HtmlElement {
        match self {
            MetaData::Novor(n) => n.to_html(),
            MetaData::Peaks(p) => p.to_html(),
            MetaData::Opair(o) => o.to_html(),
            MetaData::Fasta(f) => f.to_html(),
            MetaData::MaxQuant(m) => m.to_html(),
            MetaData::MSFragger(m) => m.to_html(),
            MetaData::Sage(s) => s.to_html(),
            MetaData::None => HtmlElement::new(HtmlTag::i).content("No metadata").clone(),
        }
    }
}

impl RenderToHtml for PeaksData {
    fn to_html(&self) -> HtmlElement {
        HtmlElement::new(HtmlTag::div)
            .children([
                HtmlElement::new(HtmlTag::p)
                    .content(format!(
                        "Additional MetaData Peaks {}",
                        self.scan.iter().map(|i| i.to_string()).join(";")
                    ))
                    .clone(),
                HtmlElement::table::<HtmlContent, _>(
                    None,
                    &[
                        &["Fraction".to_string(), self.fraction.to_optional_string()],
                        &[
                            "Source file".to_string(),
                            self.source_file.clone().to_optional_string(),
                        ],
                        &[
                            "Feature".to_string(),
                            self.feature.as_ref().to_optional_string(),
                        ],
                        &[
                            "De novo score".to_string(),
                            self.de_novo_score.to_optional_string(),
                        ],
                        &["ALC".to_string(), self.alc.to_string()],
                        &["RT".to_string(), self.rt.value.to_string()],
                        &[
                            "Predicted RT".to_string(),
                            self.predicted_rt.map(|v| v.value).to_optional_string(),
                        ],
                        &["Area".to_string(), self.area.to_optional_string()],
                        &["ppm".to_string(), self.ppm.to_string()],
                        &[
                            "Post Translational Modifications".to_string(),
                            self.ptm.to_string(),
                        ],
                        &[
                            "Accession".to_string(),
                            self.accession.as_ref().to_optional_string(),
                        ],
                        &["Version".to_string(), self.version.to_string()],
                    ],
                )
                .clone(),
            ])
            .clone()
    }
}

impl RenderToHtml for NovorData {
    fn to_html(&self) -> HtmlElement {
        HtmlElement::new(HtmlTag::div)
            .children([
                HtmlElement::new(HtmlTag::p)
                    .content(format!(
                        "Additional MetaData Novor {}",
                        self.id.unwrap_or(self.scan)
                    ))
                    .clone(),
                HtmlElement::table::<HtmlContent, _>(
                    None,
                    &[
                        &["Scan".to_string(), self.scan.to_string()],
                        &["ppm".to_string(), self.ppm.value.to_string()],
                        &["Score".to_string(), self.score.to_string()],
                        &["ID".to_string(), self.id.to_optional_string()],
                        &[
                            "Spectra ID".to_string(),
                            self.spectra_id.to_optional_string(),
                        ],
                        &["Fraction".to_string(), self.fraction.to_optional_string()],
                        &[
                            "RT".to_string(),
                            self.rt.map(|v| v.value).to_optional_string(),
                        ],
                        &[
                            "Mass error".to_string(),
                            self.mass_err.map(|v| v.value).to_optional_string(),
                        ],
                        &["Protein".to_string(), self.protein.to_optional_string()],
                        &[
                            "Protein Start".to_string(),
                            self.protein_start.to_optional_string(),
                        ],
                        &[
                            "Protein Origin".to_string(),
                            self.protein_origin.as_ref().to_optional_string(),
                        ],
                        &[
                            "Protein All".to_string(),
                            self.protein_all.as_ref().to_optional_string(),
                        ],
                        &[
                            "Database Sequence".to_string(),
                            self.database_sequence.as_ref().to_optional_string(),
                        ],
                        &["Version".to_string(), self.version.to_string()],
                    ],
                )
                .clone(),
            ])
            .clone()
    }
}

impl RenderToHtml for OpairData {
    fn to_html(&self) -> HtmlElement {
        HtmlElement::new(HtmlTag::div)
            .children([
                HtmlElement::new(HtmlTag::p)
                    .content(format!("Additional MetaData Opair {}", self.scan))
                    .clone(),
                HtmlElement::table::<HtmlContent, _>(
                    None,
                    &[
                        &["File name".to_string(), self.file_name.to_string()],
                        &["RT".to_string(), self.rt.value.to_string()],
                        &[
                            "Precursor scan number".to_string(),
                            self.precursor_scan_number.to_string(),
                        ],
                        &[
                            "Theoretical Mass".to_string(),
                            self.theoretical_mass.value.to_string(),
                        ],
                        &[
                            "Glycan Mass".to_string(),
                            self.glycan_mass.value.to_string(),
                        ],
                        &["Accession".to_string(), self.accession.to_string()],
                        &["Organism".to_string(), self.organism.to_string()],
                        &["Protein name".to_string(), self.protein_name.to_string()],
                        &[
                            "Protein location".to_string(),
                            format!("{} to {}", self.protein_location.0, self.protein_location.1),
                        ],
                        &[
                            "Flanking residues".to_string(),
                            format!(
                                "{}_(seq)_{}",
                                self.flanking_residues.0.char(),
                                self.flanking_residues.1.char()
                            ),
                        ],
                        &["Rank".to_string(), self.rank.to_string()],
                        &[
                            "Matched ion series".to_string(),
                            self.matched_ion_series.to_string(),
                        ],
                        &[
                            "Matched ion mz ratios".to_string(),
                            self.matched_ion_mz_ratios.to_string(),
                        ],
                        &[
                            "Matched ion mass error".to_string(),
                            self.matched_ion_mass_error.to_string(),
                        ],
                        &[
                            "Matched ion mass error (ppm)".to_string(),
                            self.matched_ion_ppm.to_string(),
                        ],
                        &[
                            "Matched ion intensities".to_string(),
                            self.matched_ion_intensities.to_string(),
                        ],
                        &[
                            "Matched ion counts".to_string(),
                            self.matched_ion_counts.to_string(),
                        ],
                        &["Kind".to_string(), self.kind.to_string()],
                        &["Q value".to_string(), self.q_value.to_string()],
                        &["PEP".to_string(), self.pep.to_string()],
                        &["PEP Q value".to_string(), self.pep_q_value.to_string()],
                        &[
                            "Localisation score".to_string(),
                            self.localisation_score.to_string(),
                        ],
                        &["Yion score".to_string(), self.yion_score.to_string()],
                        &[
                            "Diagnostic ion score".to_string(),
                            self.diagnostic_ion_score.to_string(),
                        ],
                        &[
                            "Plausible glycan number".to_string(),
                            self.plausible_glycan_number.to_string(),
                        ],
                        &[
                            "Total glycosylation sites".to_string(),
                            self.total_glycosylation_sites.to_string(),
                        ],
                        &[
                            "Plausible glycan composition".to_string(),
                            self.plausible_glycan_composition.to_string(),
                        ],
                        &[
                            "Plausible glycan structure".to_string(),
                            self.plausible_glycan_structure.to_string(),
                        ],
                        &[
                            "N glycan motif check passed".to_string(),
                            self.n_glycan_motif.to_string(),
                        ],
                        &["R138/144".to_string(), self.r138_144.to_string()],
                        &[
                            "Glycan localisation level".to_string(),
                            self.glycan_localisation_level.to_string(),
                        ],
                        &[
                            "Glycan peptide site specificity".to_string(),
                            self.glycan_peptide_site_specificity.to_string(),
                        ],
                        &[
                            "Glycan protein site specificity".to_string(),
                            self.glycan_protein_site_specificity.to_string(),
                        ],
                        &[
                            "All potential glycan localisations".to_string(),
                            self.all_potential_glycan_localisations.to_string(),
                        ],
                        &[
                            "All site specific localisation probabilities".to_string(),
                            self.all_site_specific_localisation_probabilities
                                .to_string(),
                        ],
                        &["Version".to_string(), self.version.to_string()],
                    ],
                )
                .clone(),
            ])
            .clone()
    }
}

impl RenderToHtml for MaxQuantData {
    fn to_html(&self) -> HtmlElement {
        HtmlElement::new(HtmlTag::div)
            .children([
                HtmlElement::new(HtmlTag::p)
                    .content(format!(
                        "Additional MetaData MaxQuant {}",
                        self.id
                            .unwrap_or(self.scan_number.first().copied().unwrap_or(0))
                    ))
                    .clone(),
                HtmlElement::table::<HtmlContent, _>(
                    None,
                    &[
                        &["Raw file".to_string(), self.raw_file.to_string()],
                        &["Scan number".to_string(), self.scan_number.iter().join(";")],
                        &[
                            "Scan index".to_string(),
                            self.scan_index.to_optional_string(),
                        ],
                        &["Modifications".to_string(), self.modifications.to_string()],
                        &["Proteins".to_string(), self.proteins.to_string()],
                        &[
                            "Fragmentation".to_string(),
                            self.fragmentation.as_ref().to_optional_string(),
                        ],
                        &[
                            "Mass analyser".to_string(),
                            self.mass_analyser.as_ref().to_optional_string(),
                        ],
                        &["Type".to_string(), self.ty.to_string()],
                        &[
                            "Scan event number".to_string(),
                            self.scan_event_number.to_optional_string(),
                        ],
                        &["Pep".to_string(), self.pep.to_string()],
                        &["Score".to_string(), self.score.to_string()],
                        &[
                            "Precursor full scan number".to_string(),
                            self.precursor.to_optional_string(),
                        ],
                        &[
                            "Precursor intensity".to_string(),
                            self.precursor_intensity.to_optional_string(),
                        ],
                        &[
                            "Precursor apex".to_string(),
                            self.precursor_apex_function
                                .map(|function| {
                                    format!(
                                        "function: {function} offset: {} offset time: {}",
                                        self.precursor_apex_offset.to_optional_string(),
                                        self.precursor_apex_offset_time.to_optional_string(),
                                    )
                                })
                                .to_optional_string(),
                        ],
                        &[
                            "Missed cleavages".to_string(),
                            self.missed_cleavages.to_optional_string(),
                        ],
                        &[
                            "Isotope index".to_string(),
                            self.isotope_index.to_optional_string(),
                        ],
                        &[
                            "Simple mass error [ppm]".to_string(),
                            self.simple_mass_error_ppm.to_optional_string(),
                        ],
                        &[
                            "Retention time".to_string(),
                            self.retention_time.map(|v| v.value).to_optional_string(),
                        ],
                        &[
                            "Number of matches".to_string(),
                            self.number_of_matches.to_optional_string(),
                        ],
                        &[
                            "Intensity coverage".to_string(),
                            self.intensity_coverage.to_optional_string(),
                        ],
                        &[
                            "Peak coverage".to_string(),
                            self.peak_coverage.to_optional_string(),
                        ],
                        &[
                            "Delta score".to_string(),
                            self.delta_score.to_optional_string(),
                        ],
                        &[
                            "Score diff".to_string(),
                            self.score_diff.to_optional_string(),
                        ],
                        &[
                            "Localization probability".to_string(),
                            self.localisation_probability.to_optional_string(),
                        ],
                        &[
                            "All modified sequences".to_string(),
                            self.all_modified_sequences
                                .as_ref()
                                .map(|v| v.iter().map(ToString::to_string).join(","))
                                .to_optional_string(),
                        ],
                        &[
                            "Protein group ids".to_string(),
                            self.protein_group_ids
                                .as_ref()
                                .map(|v| v.iter().map(ToString::to_string).join(","))
                                .to_optional_string(),
                        ],
                        &[
                            "IDs".to_string(),
                            format!(
                                "id: {} peptide: {} mod. peptide: {} evidence: {}",
                                self.id.to_optional_string(),
                                self.peptide_id.to_optional_string(),
                                self.modified_peptide_id.to_optional_string(),
                                self.evidence_id.to_optional_string()
                            ),
                        ],
                        &[
                            "Base peak intensity".to_string(),
                            self.base_peak_intensity.to_optional_string(),
                        ],
                        &[
                            "Total ion current".to_string(),
                            self.total_ion_current.to_optional_string(),
                        ],
                        &[
                            "Collision energy".to_string(),
                            self.collision_energy.to_optional_string(),
                        ],
                        &[
                            "DN sequence".to_string(),
                            self.dn_sequence.as_ref().to_optional_string(),
                        ],
                        &[
                            "DN combined score".to_string(),
                            self.dn_combined_score.to_optional_string(),
                        ],
                        &[
                            "DN mass left".to_string(),
                            self.dn_missing_mass
                                .map(|v| {
                                    format!(
                                        "N: {} other: {} C: {}",
                                        self.dn_n_mass
                                            .map(|v| display_mass(v, None),)
                                            .to_optional_string(),
                                        display_mass(v, None),
                                        self.dn_c_mass
                                            .map(|v| display_mass(v, None),)
                                            .to_optional_string(),
                                    )
                                })
                                .to_optional_string(),
                        ],
                        &[
                            "Ratio H/L".to_string(),
                            self.ration_h_l
                                .map(|v| {
                                    format!(
                                        "{v} normalised: {}",
                                        self.ration_h_l_normalised.to_optional_string()
                                    )
                                })
                                .to_optional_string(),
                        ],
                        &[
                            "Intensity".to_string(),
                            self.intensity
                                .map(|i| {
                                    format!(
                                        "{i} H: {} L: {}",
                                        self.intensity_h.to_optional_string(),
                                        self.intensity_l.to_optional_string()
                                    )
                                })
                                .to_optional_string(),
                        ],
                        &[
                            "Experiment".to_string(),
                            self.experiment.as_ref().to_optional_string(),
                        ],
                        &[
                            "Labeling state".to_string(),
                            self.labeling_state.to_optional_string(),
                        ],
                        &["Version".to_string(), self.version.to_string()],
                    ],
                )
                .clone(),
            ])
            .clone()
    }
}

impl RenderToHtml for SageData {
    fn to_html(&self) -> HtmlElement {
        HtmlElement::new(HtmlTag::div)
            .children([
                HtmlElement::new(HtmlTag::p)
                    .content(format!("Additional MetaData Sage {}", self.psm_id))
                    .clone(),
                HtmlElement::table::<HtmlContent, _>(
                    None,
                    &[
                        &["Raw file".to_string(), self.filename.to_string()],
                        &[
                            "Scan number".to_string(),
                            format!("{},{},{}", self.scan_nr.0, self.scan_nr.1, self.scan_nr.2),
                        ],
                        &["Proteins".to_string(), self.proteins.join(";")],
                        &["Rank".to_string(), self.rank.to_string()],
                        &["Decoy".to_string(), self.decoy.to_string()],
                        &[
                            "Experimental mass".to_string(),
                            display_mass(self.experimental_mass, None).to_string(),
                        ],
                        &[
                            "Missed cleavages".to_string(),
                            self.missed_cleavages.to_string(),
                        ],
                        &[
                            "Semi enzymatic".to_string(),
                            self.semi_enzymatic.to_string(),
                        ],
                        &["Isotope error".to_string(), self.isotope_error.to_string()],
                        &[
                            "Precursor error (ppm)".to_string(),
                            self.precursor_ppm.value.to_string(),
                        ],
                        &[
                            "Fragment error (average ppm)".to_string(),
                            self.precursor_ppm.value.to_string(),
                        ],
                        &["Hyperscore".to_string(), self.hyperscore.to_string()],
                        &[
                            "Hyperscore delta next".to_string(),
                            self.delta_next.to_string(),
                        ],
                        &[
                            "Hyperscore delta best".to_string(),
                            self.delta_best.to_string(),
                        ],
                        &[
                            "Retention time (min)".to_string(),
                            self.rt.value.to_string(),
                        ],
                        &[
                            "Retention time aligned".to_string(),
                            self.aligned_rt.value.to_string(),
                        ],
                        &[
                            "Retention time predicted".to_string(),
                            self.predicted_rt.value.to_string(),
                        ],
                        &[
                            "Retention time delta".to_string(),
                            self.delta_rt_model.to_string(),
                        ],
                        &["Ion mobility".to_string(), self.ion_mobility.to_string()],
                        &[
                            "Ion mobility predicted".to_string(),
                            self.predicted_mobility.to_string(),
                        ],
                        &[
                            "Ion mobility delta".to_string(),
                            self.delta_mobility.to_string(),
                        ],
                        &["Matched peaks".to_string(), self.matched_peaks.to_string()],
                        &["Longest b".to_string(), self.longest_b.to_string()],
                        &["Longest y".to_string(), self.longest_y.to_string()],
                        &[
                            "Matched intensity".to_string(),
                            self.matched_intensity_pct.value.to_string(),
                        ],
                        &[
                            "Scored candidates".to_string(),
                            self.scored_candidates.to_string(),
                        ],
                        &["Poisson".to_string(), self.poisson.to_string()],
                        &[
                            "Sage discriminant score".to_string(),
                            self.sage_discriminant_score.to_string(),
                        ],
                        &[
                            "Posterior error".to_string(),
                            self.posterior_error.to_string(),
                        ],
                        &["Spectrum q".to_string(), self.spectrum_q.to_string()],
                        &["Peptide q".to_string(), self.peptide_q.to_string()],
                        &["Protein q".to_string(), self.protein_q.to_string()],
                        &["MS2 intensity".to_string(), self.ms2_intensity.to_string()],
                    ],
                )
                .clone(),
            ])
            .clone()
    }
}

impl RenderToHtml for MSFraggerData {
    fn to_html(&self) -> HtmlElement {
        HtmlElement::new(HtmlTag::div)
            .children([
                HtmlElement::new(HtmlTag::p)
                    .content(format!("Additional MetaData MSFragger {}", self.spectrum))
                    .clone(),
                HtmlElement::table::<HtmlContent, _>(
                    None,
                    &[
                        &["Raw file".to_string(), self.spectrum.file.to_string()],
                        &["Spectrum file".to_string(), self.spectrum_file.to_string()],
                        &[
                            "Scan number".to_string(),
                            format!(
                                "{}.{}.{}",
                                self.spectrum.scan.0, self.spectrum.scan.1, self.spectrum.scan.2
                            ),
                        ],
                        &[
                            "Extended Peptide".to_string(),
                            self.extended_peptide.to_string(),
                        ],
                        &[
                            "Assigned modifications".to_string(),
                            self.assigned_modifications.to_string(),
                        ],
                        &[
                            "Experimental mass".to_string(),
                            display_mass(self.experimental_mass, None).to_string(),
                        ],
                        &[
                            "Calibrated experimental mass".to_string(),
                            display_mass(self.calibrated_experimental_mass, None).to_string(),
                        ],
                        &[
                            "Experimental m/z".to_string(),
                            self.experimental_mz.value.to_string(),
                        ],
                        &[
                            "Calibrated experimental m/z".to_string(),
                            self.calibrated_experimental_mz.value.to_string(),
                        ],
                        &[
                            "Retention time (min)".to_string(),
                            self.rt.value.to_string(),
                        ],
                        &["Expectation".to_string(), self.expectation.to_string()],
                        &["Hyperscore".to_string(), self.hyperscore.to_string()],
                        &["Next score".to_string(), self.next_score.to_string()],
                        &[
                            "Peptide Prophet probability".to_string(),
                            self.peptide_prophet_probability.to_string(),
                        ],
                        &[
                            "Missed cleavages".to_string(),
                            self.missed_cleavages.to_string(),
                        ],
                        &[
                            "Number of enzymatic termini".to_string(),
                            self.enzymatic_termini.to_string(),
                        ],
                        &["Protein start".to_string(), self.protein_start.to_string()],
                        &["Protein end".to_string(), self.protein_end.to_string()],
                        &["Intensity".to_string(), format!("{:e}", self.intensity)],
                        &["Purity".to_string(), self.purity.to_string()],
                        &["Is unique".to_string(), self.is_unique.to_string()],
                        &["Protein".to_string(), self.protein.to_string()],
                        &["Protein ID".to_string(), self.protein_id.to_string()],
                        &["Entry name".to_string(), self.entry_name.to_string()],
                        &["Gene".to_string(), self.gene.to_string()],
                        &[
                            "Protein description".to_string(),
                            self.protein_description.to_string(),
                        ],
                        &["Mapped genes".to_string(), self.mapped_genes.join(",")],
                        &[
                            "Mapped proteins".to_string(),
                            self.mapped_proteins.join(","),
                        ],
                        &[
                            "Condition".to_string(),
                            self.condition.as_ref().to_optional_string(),
                        ],
                        &[
                            "Group".to_string(),
                            self.group.as_ref().to_optional_string(),
                        ],
                    ],
                )
                .clone(),
            ])
            .clone()
    }
}

impl RenderToHtml for FastaData {
    fn to_html(&self) -> HtmlElement {
        HtmlElement::new(HtmlTag::div)
            .children([
                HtmlElement::new(HtmlTag::p)
                    .content(format!("Additional MetaData Fasta {}", self.id))
                    .clone(),
                HtmlElement::table::<HtmlContent, _>(
                    None,
                    &[&["Header".to_string(), self.full_header.to_string()]],
                )
                .clone(),
            ])
            .clone()
    }
}

trait OptionalString {
    fn to_optional_string(self) -> String;
}

impl<T: ToString> OptionalString for Option<T> {
    fn to_optional_string(self) -> String {
        self.map_or("-".to_string(), |v| v.to_string())
    }
}
