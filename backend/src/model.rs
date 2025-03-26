use serde::{Deserialize, Serialize};

use rustyms::{
    error::*,
    model::*,
    system::{mz, MassOverCharge},
    *,
};

use crate::validate::{
    display_aa_neutral_loss, display_satellite_ion, parse_aa_neutral_loss, parse_amino_acid,
    parse_satellite_ion,
};

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
    pub glycan: GlycanParameters,
}

impl TryFrom<ModelParameters> for Model {
    type Error = CustomError;
    fn try_from(value: ModelParameters) -> Result<Self, Self::Error> {
        Ok(Model::none()
            .a(value.a.try_into()?)
            .b(value.b.try_into()?)
            .c(value.c.try_into()?)
            .d(value.d.try_into()?)
            .v(value.v.try_into()?)
            .w(value.w.try_into()?)
            .x(value.x.try_into()?)
            .y(value.y.try_into()?)
            .z(value.z.try_into()?)
            .precursor(
                parse_neutral_losses(&value.precursor.0.neutral_losses)?,
                value
                    .precursor
                    .0
                    .amino_acid_neutral_losses
                    .iter()
                    .map(|rule| parse_aa_neutral_loss(rule))
                    .collect::<Result<_, _>>()?,
                (value.precursor.0.amino_acid_side_chain_losses, {
                    let selection = value
                        .precursor
                        .0
                        .amino_acid_side_chain_losses_selection
                        .iter()
                        .map(|aa| parse_amino_acid(aa))
                        .collect::<Result<Vec<_>, _>>()?;
                    (!selection.is_empty()).then_some(selection)
                }),
                value.precursor.1,
            )
            .immonium(value.immonium.0.then_some(value.immonium.1))
            .modification_specific_diagnostic_ions(
                value
                    .modification_diagnostic
                    .0
                    .then_some(value.modification_diagnostic.1),
            )
            .modification_specific_neutral_losses(value.modification_neutral)
            .allow_cross_link_cleavage(value.cleave_cross_links)
            .glycan(value.glycan.try_into()?))
    }
}

impl From<Model> for ModelParameters {
    fn from(value: Model) -> Self {
        ModelParameters {
            a: value.a.into(),
            b: value.b.into(),
            c: value.c.into(),
            d: value.d.into(),
            v: value.v.into(),
            w: value.w.into(),
            x: value.x.into(),
            y: value.y.into(),
            z: value.z.into(),
            precursor: ((&value.precursor).into(), value.precursor.3),
            immonium: value
                .immonium
                .map_or((false, ChargeRange::ONE), |c| (true, c)),
            modification_diagnostic: value
                .modification_specific_diagnostic_ions
                .map_or((false, ChargeRange::ONE), |c| (true, c)),
            modification_neutral: value.modification_specific_neutral_losses,
            cleave_cross_links: value.allow_cross_link_cleavage,
            glycan: value.glycan.into(),
        }
    }
}

#[derive(Debug, PartialEq, PartialOrd, Serialize, Deserialize)]
pub struct GlycanParameters {
    allow_structural: bool,
    compositional_range: (usize, usize),
    neutral_losses: Vec<String>,
    oxonium_charge_range: ChargeRange,
    other_charge_range: ChargeRange,
}

impl TryFrom<GlycanParameters> for GlycanModel {
    type Error = CustomError;
    fn try_from(value: GlycanParameters) -> Result<Self, Self::Error> {
        Ok(Self {
            allow_structural: value.allow_structural,
            compositional_range: value.compositional_range.0..=value.compositional_range.1,
            neutral_losses: parse_neutral_losses(&value.neutral_losses)?,
            oxonium_charge_range: value.oxonium_charge_range,
            other_charge_range: value.other_charge_range,
        })
    }
}

impl From<GlycanModel> for GlycanParameters {
    fn from(value: GlycanModel) -> Self {
        Self {
            allow_structural: value.allow_structural,
            compositional_range: (
                *value.compositional_range.start(),
                *value.compositional_range.end(),
            ),
            neutral_losses: value.neutral_losses.iter().map(|l| l.to_string()).collect(),
            oxonium_charge_range: value.oxonium_charge_range,
            other_charge_range: value.other_charge_range,
        }
    }
}

#[derive(Debug, PartialEq, PartialOrd, Serialize, Deserialize)]
pub struct PrecursorLosses {
    neutral_losses: Vec<String>,
    amino_acid_neutral_losses: Vec<String>,
    amino_acid_side_chain_losses: u8,
    amino_acid_side_chain_losses_selection: Vec<String>,
}

impl
    From<&(
        Vec<NeutralLoss>,
        Vec<(Vec<AminoAcid>, Vec<NeutralLoss>)>,
        (u8, Option<Vec<AminoAcid>>),
        ChargeRange,
    )> for PrecursorLosses
{
    fn from(
        value: &(
            Vec<NeutralLoss>,
            Vec<(Vec<AminoAcid>, Vec<NeutralLoss>)>,
            (u8, Option<Vec<AminoAcid>>),
            ChargeRange,
        ),
    ) -> Self {
        Self {
            neutral_losses: value.0.iter().map(|l| l.to_string()).collect(),
            amino_acid_neutral_losses: value
                .1
                .iter()
                .map(|l| display_aa_neutral_loss(&l.0, &l.1))
                .collect(),
            amino_acid_side_chain_losses: value.2 .0,
            amino_acid_side_chain_losses_selection: value
                .2
                 .1
                .iter()
                .flat_map(|selection| selection.iter().map(|a| a.to_string()))
                .collect(),
        }
    }
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

impl TryFrom<PrimarySeriesParameters> for PrimaryIonSeries {
    type Error = CustomError;
    fn try_from(value: PrimarySeriesParameters) -> Result<Self, Self::Error> {
        let selection = value
            .amino_acid_side_chain_losses_selection
            .iter()
            .map(|aa| parse_amino_acid(aa))
            .collect::<Result<Vec<_>, _>>()?;
        Ok(Self {
            location: value.location,
            neutral_losses: parse_neutral_losses(&value.neutral_losses)?,
            amino_acid_neutral_losses: value
                .amino_acid_neutral_losses
                .iter()
                .map(|rule| parse_aa_neutral_loss(rule))
                .collect::<Result<_, _>>()?,
            amino_acid_side_chain_losses: (
                value.amino_acid_side_chain_losses,
                (!selection.is_empty()).then_some(selection),
            ),
            charge_range: value.charge_range,
            allowed_variants: value.variants,
        })
    }
}

impl From<PrimaryIonSeries> for PrimarySeriesParameters {
    fn from(value: PrimaryIonSeries) -> Self {
        PrimarySeriesParameters {
            location: value.location,
            neutral_losses: value.neutral_losses.iter().map(|l| l.to_string()).collect(),
            amino_acid_neutral_losses: value
                .amino_acid_neutral_losses
                .iter()
                .map(|l| display_aa_neutral_loss(&l.0, &l.1))
                .collect(),
            amino_acid_side_chain_losses: value.amino_acid_side_chain_losses.0,
            amino_acid_side_chain_losses_selection: value
                .amino_acid_side_chain_losses
                .1
                .iter()
                .flat_map(|selection| selection.iter().map(|a| a.to_string()))
                .collect(),
            charge_range: value.charge_range,
            variants: value.allowed_variants,
        }
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

impl TryFrom<SatelliteSeriesParameters> for SatelliteIonSeries {
    type Error = CustomError;
    fn try_from(value: SatelliteSeriesParameters) -> Result<Self, Self::Error> {
        let selection = value
            .amino_acid_side_chain_losses_selection
            .iter()
            .map(|aa| parse_amino_acid(aa))
            .collect::<Result<Vec<_>, _>>()?;
        Ok(SatelliteIonSeries {
            location: SatelliteLocation {
                rules: value
                    .location_rules
                    .into_iter()
                    .map(|rule| parse_satellite_ion(&rule))
                    .collect::<Result<Vec<_>, _>>()?,
                base: value.location_base,
            },
            neutral_losses: parse_neutral_losses(&value.neutral_losses)?,
            amino_acid_neutral_losses: value
                .amino_acid_neutral_losses
                .iter()
                .map(|rule| parse_aa_neutral_loss(rule))
                .collect::<Result<_, _>>()?,
            amino_acid_side_chain_losses: (
                value.amino_acid_side_chain_losses,
                (!selection.is_empty()).then_some(selection),
            ),
            charge_range: value.charge_range,
            allowed_variants: value.variants,
        })
    }
}

impl From<SatelliteIonSeries> for SatelliteSeriesParameters {
    fn from(value: SatelliteIonSeries) -> Self {
        SatelliteSeriesParameters {
            location_base: value.location.base,
            location_rules: value
                .location
                .rules
                .iter()
                .map(|rule| display_satellite_ion(&rule.0, rule.1))
                .collect(),
            neutral_losses: value.neutral_losses.iter().map(|l| l.to_string()).collect(),
            amino_acid_neutral_losses: value
                .amino_acid_neutral_losses
                .iter()
                .map(|l| display_aa_neutral_loss(&l.0, &l.1))
                .collect(),
            amino_acid_side_chain_losses: value.amino_acid_side_chain_losses.0,
            amino_acid_side_chain_losses_selection: value
                .amino_acid_side_chain_losses
                .1
                .iter()
                .flat_map(|selection| selection.iter().map(|a| a.to_string()))
                .collect(),
            charge_range: value.charge_range,
            variants: value.allowed_variants,
        }
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
    pub fn get_model(self, name: &str) -> Result<Model, CustomError> {
        Ok(match name {
            "all" => Model::all(),
            "ethcd" => Model::ethcd(),
            "hot_eacid" => Model::hot_eacid(),
            "ead" => Model::ead(),
            "cidhcd" => Model::cid_hcd(),
            "etd" => Model::etd(),
            "td_etd" => Model::td_etd(),
            "none" => Model::none(),
            "custom" => self.try_into()?,
            _ => Model::all(),
        })
    }

    pub fn create_model(
        self,
        name: &str,
        tolerance: (f64, &str),
        mz_range: (Option<f64>, Option<f64>),
    ) -> Result<Model, CustomError> {
        let mut model = self.get_model(name)?;
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
