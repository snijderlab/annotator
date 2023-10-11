use itertools::Itertools;
use rustyms::identifications::{MetaData, NovorData, PeaksData};

use crate::html_builder::{HtmlContent, HtmlElement, HtmlTag};

pub trait RenderToHtml {
    fn to_html(&self) -> HtmlElement;
}

impl RenderToHtml for MetaData {
    fn to_html(&self) -> HtmlElement {
        match self {
            MetaData::Novor(n) => n.to_html(),
            MetaData::Peaks(p) => p.to_html(),
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

trait OptionalString {
    fn to_optional_string(self) -> String;
}

impl<T: ToString> OptionalString for Option<T> {
    fn to_optional_string(self) -> String {
        self.map_or("-".to_string(), |v| v.to_string())
    }
}
