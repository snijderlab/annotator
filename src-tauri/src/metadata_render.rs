use itertools::Itertools;
use rustyms::identifications::{FastaData, MetaData, NovorData, OpairData, PeaksData};

use crate::html_builder::{HtmlContent, HtmlElement, HtmlTag};

pub trait RenderToHtml {
    fn to_html(&self) -> HtmlElement;
}

impl RenderToHtml for MetaData {
    fn to_html(&self) -> HtmlElement {
        match self {
            MetaData::Novor(n) => n.to_html(),
            MetaData::Peaks(p) => p.to_html(),
            MetaData::Opair(o) => o.to_html(),
            MetaData::Fasta(f) => f.to_html(),
        }
    }
}

impl RenderToHtml for PeaksData {
    fn to_html(&self) -> HtmlElement {
        HtmlElement::new(HtmlTag::details).children([
            HtmlElement::new(HtmlTag::summary).content(format!(
                "MetaData Peaks read {}",
                self.scan.iter().map(|i| i.to_string()).join(";")
            )),
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
                    &["Length".to_string(), self.length.to_string()],
                    &["mz".to_string(), self.mz.to_string()],
                    &["z".to_string(), self.z.to_string()],
                    &["RT".to_string(), self.rt.to_string()],
                    &[
                        "Predicted RT".to_string(),
                        self.predicted_rt.to_optional_string(),
                    ],
                    &["Area".to_string(), self.area.to_optional_string()],
                    &["Mass".to_string(), self.mass.to_string()],
                    &["ppm".to_string(), self.ppm.to_string()],
                    &[
                        format!("Tag (length {})", self.tag_length),
                        self.tag.to_string(),
                    ],
                    &[
                        "Post Translational Modifications".to_string(),
                        self.ptm.to_string(),
                    ],
                    &["Mode".to_string(), self.mode.to_string()],
                    &[
                        "Accession".to_string(),
                        self.accession.as_ref().to_optional_string(),
                    ],
                    &["Version".to_string(), self.version.to_string()],
                ],
            ),
        ])
    }
}

impl RenderToHtml for NovorData {
    fn to_html(&self) -> HtmlElement {
        HtmlElement::new(HtmlTag::details).children([
            HtmlElement::new(HtmlTag::summary).content(format!(
                "MetaData Novor read {}",
                self.id.unwrap_or(self.scan)
            )),
            HtmlElement::table::<HtmlContent, _>(
                None,
                &[
                    &["Scan".to_string(), self.scan.to_string()],
                    &["mz".to_string(), self.mz.to_string()],
                    &["z".to_string(), self.z.to_string()],
                    &["Mass".to_string(), self.mass.to_string()],
                    &["ppm".to_string(), self.ppm.to_string()],
                    &["Score".to_string(), self.score.to_string()],
                    &["ID".to_string(), self.id.to_optional_string()],
                    &[
                        "Spectra ID".to_string(),
                        self.spectra_id.to_optional_string(),
                    ],
                    &["Fraction".to_string(), self.fraction.to_optional_string()],
                    &["RT".to_string(), self.rt.to_optional_string()],
                    &["Mass error".to_string(), self.mass_err.to_optional_string()],
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
            ),
        ])
    }
}

impl RenderToHtml for OpairData {
    fn to_html(&self) -> HtmlElement {
        HtmlElement::new(HtmlTag::details).children([
            HtmlElement::new(HtmlTag::summary)
                .content(format!("MetaData Opair read {}", self.scan)),
            HtmlElement::table::<HtmlContent, _>(
                None,
                &[
                    &["File name".to_string(), self.file_name.to_string()],
                    &["RT".to_string(), self.rt.to_string()],
                    &[
                        "Precursor scan number".to_string(),
                        self.precursor_scan_number.to_string(),
                    ],
                    &["mz".to_string(), self.mz.to_string()],
                    &["z".to_string(), self.z.to_string()],
                    &["Mass".to_string(), self.mass.to_string()],
                    &[
                        "Theoretical Mass".to_string(),
                        self.theoretical_mass.to_string(),
                    ],
                    &["Glycan Mass".to_string(), self.glycan_mass.to_string()],
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
            ),
        ])
    }
}

impl RenderToHtml for FastaData {
    fn to_html(&self) -> HtmlElement {
        HtmlElement::new(HtmlTag::details).children([
            HtmlElement::new(HtmlTag::summary).content(format!("MetaData Fasta read {}", self.id)),
            HtmlElement::table::<HtmlContent, _>(
                None,
                &[&["Header".to_string(), self.full_header.to_string()]],
            ),
        ])
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
