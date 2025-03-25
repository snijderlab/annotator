use itertools::Itertools;
use serde::{Deserialize, Serialize};

use rustyms::{
    error::*,
    model::*,
    system::{mz, MassOverCharge},
    *,
};

use crate::validate::{parse_aa_neutral_loss, parse_amino_acid, parse_satellite_ion};

#[derive(Debug, PartialEq, PartialOrd, Serialize, Deserialize)]
pub struct ModelParameters {
    pub a: PrimarySeriesParameters,
    pub b: PrimarySeriesParameters,
    pub c: PrimarySeriesParameters,
    pub d: SatelliteSeriesParameters,
    pub v: SatelliteSeriesParameters,
    pub w: SatelliteSeriesParameters,
    pub x: PrimarySeriesParameters,
    pub y: PrimarySeriesParameters,
    pub z: PrimarySeriesParameters,
    pub precursor: (PrecursorLosses, ChargeRange),
    pub immonium: (bool, ChargeRange),
    pub modification_diagnostic: (bool, ChargeRange),
    pub modification_neutral: bool,
    pub cleave_cross_links: bool,
    pub glycan: (bool, (usize, usize), Vec<String>, ChargeRange, ChargeRange),
}

#[derive(Debug, PartialEq, PartialOrd, Serialize, Deserialize)]
pub struct PrecursorLosses {
    neutral_losses: Vec<String>,
    amino_acid_neutral_losses: Vec<String>,
    amino_acid_side_chain_losses: u8,
    amino_acid_side_chain_losses_selection: Vec<String>,
}

#[derive(Debug, PartialEq, PartialOrd, Serialize, Deserialize)]
pub struct PrimarySeriesParameters {
    location: Location,
    neutral_losses: Vec<String>,
    amino_acid_neutral_losses: Vec<String>,
    amino_acid_side_chain_losses: u8,
    amino_acid_side_chain_losses_selection: Vec<String>,
    charge_range: ChargeRange,
    variants: Vec<i8>,
}

impl PrimarySeriesParameters {
    fn into_final(self) -> Result<PrimaryIonSeries, CustomError> {
        let selection = self
            .amino_acid_side_chain_losses_selection
            .iter()
            .map(|aa| parse_amino_acid(aa))
            .collect::<Result<Vec<_>, _>>()?;
        Ok(PrimaryIonSeries {
            location: self.location,
            neutral_losses: parse_neutral_losses(&self.neutral_losses)?,
            amino_acid_neutral_losses: self
                .amino_acid_neutral_losses
                .iter()
                .map(|rule| parse_aa_neutral_loss(rule))
                .collect::<Result<_, _>>()?,
            amino_acid_side_chain_losses: (
                self.amino_acid_side_chain_losses,
                (!selection.is_empty()).then_some(selection),
            ),
            charge_range: self.charge_range,
            allowed_variants: self.variants,
        })
    }
}

#[derive(Debug, PartialEq, PartialOrd, Serialize, Deserialize)]
pub struct SatelliteSeriesParameters {
    location_rules: Vec<String>,
    location_base: Option<u8>,
    neutral_losses: Vec<String>,
    amino_acid_neutral_losses: Vec<String>,
    amino_acid_side_chain_losses: u8,
    amino_acid_side_chain_losses_selection: Vec<String>,
    charge_range: ChargeRange,
    variants: Vec<i8>,
}

impl SatelliteSeriesParameters {
    fn into_final(self) -> Result<SatelliteIonSeries, CustomError> {
        let selection = self
            .amino_acid_side_chain_losses_selection
            .iter()
            .map(|aa| parse_amino_acid(aa))
            .collect::<Result<Vec<_>, _>>()?;
        Ok(SatelliteIonSeries {
            location: SatelliteLocation {
                rules: self
                    .location_rules
                    .into_iter()
                    .map(|rule| parse_satellite_ion(&rule))
                    .collect::<Result<Vec<_>, _>>()?,
                base: self.location_base,
            },
            neutral_losses: parse_neutral_losses(&self.neutral_losses)?,
            amino_acid_neutral_losses: self
                .amino_acid_neutral_losses
                .iter()
                .map(|rule| parse_aa_neutral_loss(rule))
                .collect::<Result<_, _>>()?,
            amino_acid_side_chain_losses: (
                self.amino_acid_side_chain_losses,
                (!selection.is_empty()).then_some(selection),
            ),
            charge_range: self.charge_range,
            allowed_variants: self.variants,
        })
    }
}

fn parse_neutral_losses(neutral_losses: &[String]) -> Result<Vec<NeutralLoss>, CustomError> {
    neutral_losses
        .iter()
        .filter(|n| !n.is_empty())
        .map(|n| n.parse::<NeutralLoss>())
        .collect()
}

impl ModelParameters {
    pub fn create_model(
        self,
        name: &str,
        tolerance: (f64, &str),
        mz_range: (Option<f64>, Option<f64>),
    ) -> Result<Model, CustomError> {
        let mut model = match name {
            "all" => Model::all(),
            "ethcd" => Model::ethcd(),
            "hot_eacid" => Model::hot_eacid(),
            "ead" => Model::ead(),
            "cidhcd" => Model::cid_hcd(),
            "etd" => Model::etd(),
            "td_etd" => Model::td_etd(),
            "none" => Model::none(),
            "custom" => Model::none()
                .a(self.a.into_final()?)
                .b(self.b.into_final()?)
                .c(self.c.into_final()?)
                .d(self.d.into_final()?)
                .v(self.v.into_final()?)
                .w(self.w.into_final()?)
                .x(self.x.into_final()?)
                .y(self.y.into_final()?)
                .z(self.z.into_final()?)
                .precursor(
                    parse_neutral_losses(&self.precursor.0.neutral_losses)?,
                    self.precursor
                        .0
                        .amino_acid_neutral_losses
                        .iter()
                        .map(|rule| parse_aa_neutral_loss(rule))
                        .collect::<Result<_, _>>()?,
                    (self.precursor.0.amino_acid_side_chain_losses, {
                        let selection = self
                            .precursor
                            .0
                            .amino_acid_side_chain_losses_selection
                            .iter()
                            .map(|aa| parse_amino_acid(aa))
                            .collect::<Result<Vec<_>, _>>()?;
                        (!selection.is_empty()).then_some(selection)
                    }),
                    self.precursor.1,
                )
                .immonium(self.immonium.0.then_some(self.immonium.1))
                .modification_specific_diagnostic_ions(
                    self.modification_diagnostic
                        .0
                        .then_some(self.modification_diagnostic.1),
                )
                .modification_specific_neutral_losses(self.modification_neutral)
                .allow_cross_link_cleavage(self.cleave_cross_links)
                .glycan(GlycanModel {
                    allow_structural: self.glycan.0,
                    compositional_range: self.glycan.1 .0..=self.glycan.1 .1,
                    neutral_losses: parse_neutral_losses(&self.glycan.2)?,
                    oxonium_charge_range: self.glycan.3,
                    other_charge_range: self.glycan.4,
                }),
            _ => Model::all(),
        };
        if tolerance.1 == "ppm" {
            model.tolerance = Tolerance::new_ppm(tolerance.0);
        } else if tolerance.1 == "th" {
            model.tolerance =
                Tolerance::new_absolute(MassOverCharge::new::<rustyms::system::mz>(tolerance.0));
        } else {
            return Err(CustomError::error(
                "Invalid tolerance unit",
                "",
                Context::None,
            ));
        }
        let min = mz_range.0.unwrap_or(0.0);
        let max = mz_range.1.unwrap_or(f64::MAX);
        if min > max {
            return Err(CustomError::error(
                "m/z range invalid",
                "The minimal value is less then the maximal value",
                Context::None,
            ));
        }
        model.mz_range = MassOverCharge::new::<mz>(min)..=MassOverCharge::new::<mz>(max);
        Ok(model)
    }
}
