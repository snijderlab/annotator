use std::{io::BufWriter, sync::Mutex};

use serde::{Deserialize, Serialize};

use rustyms::{
    error::*,
    fragment::FragmentKind,
    model::*,
    system::{MassOverCharge, mz},
    *,
};
use tauri::Manager;

use crate::{
    ModifiableState,
    state::State,
    validate::{
        display_aa_neutral_loss, display_monosaccharide_neutral_loss, display_satellite_ion,
        parse_aa_neutral_loss, parse_amino_acid, parse_monosaccharide_neutral_loss,
        parse_satellite_ion,
    },
};

pub fn get_models(state: &State) -> (usize, Vec<(bool, &str, &FragmentationModel)>) {
    let mut output = vec![
        (true, "All", FragmentationModel::all()),
        (true, "EtHCD/EtCAD", FragmentationModel::ethcd()),
        (true, "EAciD", FragmentationModel::eacid()),
        (true, "EAD", FragmentationModel::ead()),
        (true, "CID/HCD", FragmentationModel::cid_hcd()),
        (true, "ETD", FragmentationModel::etd()),
        (true, "Top-down ETD", FragmentationModel::td_etd()),
        (true, "UVPD", FragmentationModel::uvpd()),
        (true, "None", FragmentationModel::none()),
    ];
    let built_in_length = output.len();
    output.extend(state.models.iter().map(|(n, m)| (false, n.as_str(), m)));
    (built_in_length, output)
}

#[tauri::command]
pub fn get_custom_models(
    state: ModifiableState,
) -> Result<Vec<(bool, usize, String)>, &'static str> {
    let state = state.lock().map_err(|_| "Could not lock mutex")?;
    Ok(get_models(&state)
        .1
        .iter()
        .enumerate()
        .map(|(index, (built_in, name, _))| (*built_in, index, name.to_string()))
        .collect())
}

#[tauri::command]
pub fn duplicate_custom_model(
    id: usize,
    state: ModifiableState,
) -> Result<(usize, bool, String, ModelParameters), &'static str> {
    let mut locked_state = state.lock().map_err(|_| "Could not lock mutex")?;
    let models = get_models(&locked_state).1;
    if let Some((built_in, name, model)) = models.get(id).cloned() {
        let model = (*model).clone();
        let name = name.to_string();
        let models_len = models.len();
        locked_state.models.push((name.clone(), model.clone()));
        drop(locked_state);
        Ok((models_len, built_in, name, model.into()))
    } else {
        Err("Could not find specified id")
    }
}

#[tauri::command]
pub fn delete_custom_model(id: usize, state: ModifiableState) -> Result<(), &'static str> {
    let mut state = state.lock().map_err(|_| "Could not lock mutex")?;
    let (offset, models) = get_models(&state);
    if id < models.len() {
        let (built_in, _, _) = models[id];
        drop(models);
        if built_in {
            Err("Can not delete built in model")
        } else {
            state.models.remove(id - offset);
            Ok(())
        }
    } else {
        Err("Could not find specified id")
    }
}

#[tauri::command]
pub fn get_custom_model(
    id: usize,
    state: ModifiableState,
) -> Result<(usize, String, ModelParameters), &'static str> {
    let state = state.lock().map_err(|_| "Could not lock mutex")?;
    let models = get_models(&state);
    if let Some((_, name, model)) = models.1.get(id) {
        Ok((id, name.to_string(), (*model).clone().into()))
    } else {
        Err("Given index does not exist")
    }
}

#[tauri::command]
pub async fn update_model(
    id: usize,
    name: String,
    model: ModelParameters,
    app: tauri::AppHandle,
) -> Result<(), CustomError> {
    let model = (name, model.try_into()?);

    if let Ok(mut state) = app.state::<Mutex<State>>().lock() {
        let (built_in_models_length, models) = get_models(&state);
        let index = id.saturating_sub(built_in_models_length);
        if index < state.models.len() {
            let (built_in, _, _) = models[id];
            if built_in {
                return Err(CustomError::error(
                    "Can not update built in model",
                    "A built in model cannot be edited",
                    Context::None,
                ));
            } else {
                state.models[index] = model;
            }
        } else {
            state.models.push(model);
        }

        // Store mods config file
        let path = app
            .path()
            .app_config_dir()
            .map(|dir| dir.join(crate::CUSTOM_MODELS_FILE))
            .map_err(|e| {
                CustomError::error(
                    "Cannot find app data directory",
                    e.to_string(),
                    Context::None,
                )
            })?;
        let parent = path.parent().ok_or_else(|| {
            CustomError::error(
                "Custom models configuration does not have a valid directory",
                "Please report",
                Context::show(path.to_string_lossy()),
            )
        })?;
        std::fs::create_dir_all(parent).map_err(|err| {
            CustomError::error(
                "Could not create parent directories for custom models configuration file",
                err,
                Context::show(parent.to_string_lossy()),
            )
        })?;
        let file = BufWriter::new(std::fs::File::create(&path).map_err(|err| {
            CustomError::error(
                "Could not open custom models configuration file",
                err,
                Context::show(path.to_string_lossy()),
            )
        })?);
        serde_json::to_writer_pretty(file, &state.models).map_err(|err| {
            CustomError::error(
                "Could not write custom models to configuration file",
                err,
                Context::None,
            )
        })?;

        Ok(())
    } else {
        Err(CustomError::error(
            "State locked",
            "Cannot unlock the mutable state, are you doing many things in parallel?",
            Context::None,
        ))
    }
}

#[derive(Debug, PartialEq, Serialize, Deserialize)]
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

impl TryFrom<ModelParameters> for FragmentationModel {
    type Error = CustomError;
    fn try_from(value: ModelParameters) -> Result<Self, Self::Error> {
        Ok(FragmentationModel::none()
            .clone()
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
            .immonium(value.immonium.0.then_some((value.immonium.1, Vec::new())))
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

impl From<FragmentationModel> for ModelParameters {
    fn from(value: FragmentationModel) -> Self {
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
                .map_or((false, ChargeRange::ONE), |(c, _)| (true, c)),
            modification_diagnostic: value
                .modification_specific_diagnostic_ions
                .map_or((false, ChargeRange::ONE), |c| (true, c)),
            modification_neutral: value.modification_specific_neutral_losses,
            cleave_cross_links: value.allow_cross_link_cleavage,
            glycan: value.glycan.into(),
        }
    }
}

#[derive(Debug, PartialEq, Serialize, Deserialize)]
pub struct GlycanParameters {
    allow_structural: bool,
    compositional_range: (usize, usize),
    neutral_losses: Vec<String>,
    diagnostic_neutral_losses: Vec<String>,
    default_glycan_peptide_fragment: GlycanPeptideFragment,
    specific_glycan_peptide_fragment: Vec<(Vec<char>, Vec<FragmentKind>, GlycanPeptideFragment)>,
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
            specific_neutral_losses: value
                .diagnostic_neutral_losses
                .into_iter()
                .map(|v| parse_monosaccharide_neutral_loss(&v))
                .collect::<Result<Vec<_>, _>>()?,
            default_peptide_fragment: value.default_glycan_peptide_fragment,
            peptide_fragment_rules: value
                .specific_glycan_peptide_fragment
                .into_iter()
                .map(|(aas, fragments, settings)| {
                    (
                        aas.iter()
                            .filter_map(|aa| AminoAcid::try_from(*aa).ok())
                            .collect(),
                        fragments,
                        settings,
                    )
                })
                .collect(),
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
            diagnostic_neutral_losses: value
                .specific_neutral_losses
                .into_iter()
                .map(display_monosaccharide_neutral_loss)
                .collect(),
            default_glycan_peptide_fragment: value.default_peptide_fragment,
            specific_glycan_peptide_fragment: value
                .peptide_fragment_rules
                .into_iter()
                .map(|(aas, fragments, settings)| {
                    (
                        aas.iter().filter_map(|aa| aa.one_letter_code()).collect(),
                        fragments,
                        settings,
                    )
                })
                .collect(),
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
            amino_acid_side_chain_losses: value.2.0,
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

pub fn parameters(
    tolerance: (f64, &str),
    mz_range: (Option<f64>, Option<f64>),
) -> Result<MatchingParameters, CustomError> {
    let mut parameters = MatchingParameters::default();
    if tolerance.1 == "ppm" {
        parameters.tolerance = Tolerance::new_ppm(tolerance.0);
    } else if tolerance.1 == "th" {
        parameters.tolerance =
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
    parameters.mz_range = MassOverCharge::new::<mz>(min)..=MassOverCharge::new::<mz>(max);
    Ok(parameters)
}
