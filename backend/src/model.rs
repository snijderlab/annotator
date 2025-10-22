use std::{io::BufWriter, sync::Mutex};

use context_error::{BasicKind, BoxedError, Context, CreateError, FullErrorContent};
use mzannotate::{
    annotation::model::{
        BuiltInFragmentationModel, GlycanModel, Location, PrimaryIonSeries, SatelliteIonSeries,
        SatelliteLocation,
    },
    prelude::*,
};
use mzcore::{
    chemistry::{ChargeRange, NeutralLoss},
    glycan::{BackboneFragmentKind, GlycanPeptideFragment},
    prelude::*,
    quantities::Tolerance,
    system::{MassOverCharge, thomson},
};
use serde::{Deserialize, Serialize};
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

pub fn get_models(
    state: &State,
) -> (
    usize,
    Vec<(Option<BuiltInFragmentationModel>, &str, &FragmentationModel)>,
) {
    let mut output = vec![
        (
            Some(BuiltInFragmentationModel::All),
            "All",
            FragmentationModel::all(),
        ),
        (
            Some(BuiltInFragmentationModel::CID),
            "CID",
            FragmentationModel::cid(),
        ),
        (
            Some(BuiltInFragmentationModel::ETD),
            "ETD",
            FragmentationModel::etd(),
        ),
        (
            Some(BuiltInFragmentationModel::ETD_TD),
            "Top-down ETD",
            FragmentationModel::td_etd(),
        ),
        (
            Some(BuiltInFragmentationModel::ETciD),
            "ETciD",
            FragmentationModel::etcid(),
        ),
        (
            Some(BuiltInFragmentationModel::EAD),
            "EAD",
            FragmentationModel::ead(),
        ),
        (
            Some(BuiltInFragmentationModel::EAciD),
            "EAciD",
            FragmentationModel::eacid(),
        ),
        (
            Some(BuiltInFragmentationModel::UVPD),
            "UVPD",
            FragmentationModel::uvpd(),
        ),
        (
            Some(BuiltInFragmentationModel::None),
            "None",
            FragmentationModel::none(),
        ),
    ];
    let built_in_length = output.len();
    output.extend(
        state
            .custom_models
            .iter()
            .map(|(n, m)| (None, n.as_str(), m)),
    );
    (built_in_length, output)
}

const BUILT_IN_MODELS: &[&[&str]] = &[
    &["all"],
    &["cid", "hcd"],
    &["etd"],
    &["td_etd"],
    &["ethcd", "etcid", "etcad"],
    &["ead"],
    &["eacid"],
    &["uvpd"],
    &["none"],
];

pub fn get_built_in_index(built_in: BuiltInFragmentationModel) -> usize {
    match built_in {
        BuiltInFragmentationModel::All => 0,
        BuiltInFragmentationModel::CID => 1,
        BuiltInFragmentationModel::ETD => 2,
        BuiltInFragmentationModel::ETD_TD => 3,
        BuiltInFragmentationModel::ETciD => 4,
        BuiltInFragmentationModel::EAD => 5,
        BuiltInFragmentationModel::EAciD => 6,
        BuiltInFragmentationModel::UVPD => 7,
        BuiltInFragmentationModel::None => 8,
    }
}

pub fn get_model_index(
    custom_models: &[(String, FragmentationModel)],
    name: &str,
) -> Option<usize> {
    for (index, options) in BUILT_IN_MODELS.iter().enumerate() {
        if options.iter().any(|o| o.eq_ignore_ascii_case(name)) {
            return Some(index);
        }
    }
    for (index, option) in custom_models.iter().enumerate() {
        if option.0.eq_ignore_ascii_case(name) {
            return Some(index);
        }
    }
    None
}

#[tauri::command]
pub fn get_custom_models(
    state: ModifiableState,
) -> Result<
    (
        Vec<(Option<BuiltInFragmentationModel>, usize, String)>,
        Option<(String, Vec<String>)>,
    ),
    &'static str,
> {
    let state = state.lock().map_err(|_| "Could not lock mutex")?;
    Ok((
        get_models(&state)
            .1
            .iter()
            .enumerate()
            .map(|(index, (built_in, name, _))| (*built_in, index, name.to_string()))
            .collect(),
        state.custom_models_error.clone(),
    ))
}

#[tauri::command]
pub fn duplicate_custom_model(
    id: usize,
    state: ModifiableState,
) -> Result<
    (
        usize,
        Option<BuiltInFragmentationModel>,
        String,
        ModelParameters,
    ),
    &'static str,
> {
    let mut locked_state = state.lock().map_err(|_| "Could not lock mutex")?;
    let models = get_models(&locked_state).1;
    if let Some((built_in, name, model)) = models.get(id).cloned() {
        let model = (*model).clone();
        let name = name.to_string();
        let models_len = models.len();
        locked_state
            .custom_models
            .push((name.clone(), model.clone()));
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
        if built_in.is_some() {
            Err("Can not delete built in model")
        } else {
            state.custom_models.remove(id - offset);
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
) -> Result<(), String> {
    let model = (
        name,
        model
            .try_into()
            .map_err(|err: BoxedError<'static, BasicKind>| err.to_html(false))?,
    );

    if let Ok(mut state) = app.state::<Mutex<State>>().lock() {
        let (built_in_models_length, models) = get_models(&state);
        let index = id.saturating_sub(built_in_models_length);
        if index < state.custom_models.len() {
            let (built_in, _, _) = models[id];
            if built_in.is_some() {
                return Err(BoxedError::small(
                    BasicKind::Error,
                    "Can not update built in model",
                    "A built in model cannot be edited",
                )
                .to_html(false));
            } else {
                state.custom_models[index] = model;
            }
        } else {
            state.custom_models.push(model);
        }

        // Store mods config file
        let path = app
            .path()
            .app_config_dir()
            .map(|dir| dir.join(crate::CUSTOM_MODELS_FILE))
            .map_err(|e| {
                BoxedError::small(
                    BasicKind::Error,
                    "Cannot find app data directory",
                    e.to_string(),
                )
                .to_html(false)
            })?;
        let parent = path.parent().ok_or_else(|| {
            BoxedError::new(
                BasicKind::Error,
                "Custom models configuration does not have a valid directory",
                "Please report",
                Context::show(path.to_string_lossy()).to_owned(),
            )
            .to_html(false)
        })?;
        std::fs::create_dir_all(parent).map_err(|err| {
            BoxedError::new(
                BasicKind::Error,
                "Could not create parent directories for custom models configuration file",
                err.to_string(),
                Context::show(parent.to_string_lossy()).to_owned(),
            )
            .to_html(false)
        })?;
        let file = BufWriter::new(std::fs::File::create(&path).map_err(|err| {
            BoxedError::new(
                BasicKind::Error,
                "Could not open custom models configuration file",
                err.to_string(),
                Context::show(path.to_string_lossy()).to_owned(),
            )
            .to_html(false)
        })?);
        serde_json::to_writer_pretty(file, &state.custom_models).map_err(|err| {
            BoxedError::small(
                BasicKind::Error,
                "Could not write custom models to configuration file",
                err.to_string(),
            )
            .to_html(false)
        })?;

        Ok(())
    } else {
        Err(BoxedError::small(
            BasicKind::Error,
            "State locked",
            "Cannot unlock the mutable state, are you doing many things in parallel?",
        )
        .to_html(false))
    }
}

#[derive(Debug, Deserialize, PartialEq, Serialize)]
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
    pub immonium: (bool, ChargeRange, Vec<String>),
    pub modification_diagnostic: (bool, ChargeRange),
    pub modification_neutral: bool,
    pub cleave_cross_links: bool,
    pub glycan: GlycanParameters,
}

impl TryFrom<ModelParameters> for FragmentationModel {
    type Error = BoxedError<'static, BasicKind>;
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
                parse_neutral_losses(&value.precursor.0.neutral_losses)
                    .map_err(BoxedError::to_owned)?,
                value
                    .precursor
                    .0
                    .amino_acid_neutral_losses
                    .iter()
                    .map(|rule| parse_aa_neutral_loss(rule).map_err(BoxedError::to_owned))
                    .collect::<Result<_, _>>()?,
                (value.precursor.0.amino_acid_side_chain_losses, {
                    let selection = value
                        .precursor
                        .0
                        .amino_acid_side_chain_losses_selection
                        .iter()
                        .map(|aa| parse_amino_acid(aa).map_err(BoxedError::to_owned))
                        .collect::<Result<Vec<_>, _>>()?;
                    (!selection.is_empty()).then_some(selection)
                }),
                value.precursor.1,
            )
            .immonium(
                value.immonium.0.then_some((
                    value.immonium.1,
                    value
                        .immonium
                        .2
                        .iter()
                        .map(|rule| parse_aa_neutral_loss(rule).map_err(BoxedError::to_owned))
                        .collect::<Result<Vec<_>, _>>()?,
                )),
            )
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
            immonium: value.immonium.map_or(
                (false, ChargeRange::ONE, Vec::new()),
                |(c, losses)| {
                    (
                        true,
                        c,
                        losses
                            .iter()
                            .map(|(aa, loss)| display_aa_neutral_loss(aa, loss))
                            .collect(),
                    )
                },
            ),
            modification_diagnostic: value
                .modification_specific_diagnostic_ions
                .map_or((false, ChargeRange::ONE), |c| (true, c)),
            modification_neutral: value.modification_specific_neutral_losses,
            cleave_cross_links: value.allow_cross_link_cleavage,
            glycan: value.glycan.into(),
        }
    }
}

#[derive(Debug, Deserialize, PartialEq, Serialize)]
pub struct GlycanParameters {
    allow_structural: bool,
    compositional_range: (usize, usize),
    neutral_losses: Vec<String>,
    diagnostic_neutral_losses: Vec<String>,
    default_glycan_peptide_fragment: GlycanPeptideFragment,
    specific_glycan_peptide_fragment:
        Vec<(Vec<char>, Vec<BackboneFragmentKind>, GlycanPeptideFragment)>,
    oxonium_charge_range: ChargeRange,
    other_charge_range: ChargeRange,
}

impl TryFrom<GlycanParameters> for GlycanModel {
    type Error = BoxedError<'static, BasicKind>;
    fn try_from(value: GlycanParameters) -> Result<Self, Self::Error> {
        Ok(Self {
            allow_structural: value.allow_structural,
            compositional_range: (
                Some(value.compositional_range.0),
                Some(value.compositional_range.1),
            ),
            neutral_losses: parse_neutral_losses(&value.neutral_losses)
                .map_err(BoxedError::to_owned)?,
            specific_neutral_losses: value
                .diagnostic_neutral_losses
                .into_iter()
                .map(|v| parse_monosaccharide_neutral_loss(&v).map_err(BoxedError::to_owned))
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
                value.compositional_range.0.unwrap_or(0),
                value.compositional_range.1.unwrap_or(usize::MAX),
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

#[derive(Debug, Deserialize, PartialEq, PartialOrd, Serialize)]
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

#[derive(Debug, Deserialize, PartialEq, PartialOrd, Serialize)]
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
    type Error = BoxedError<'static, BasicKind>;
    fn try_from(value: PrimarySeriesParameters) -> Result<Self, Self::Error> {
        let selection = value
            .amino_acid_side_chain_losses_selection
            .iter()
            .map(|aa| parse_amino_acid(aa).map_err(BoxedError::to_owned))
            .collect::<Result<Vec<_>, _>>()?;
        Ok(Self {
            location: value.location,
            neutral_losses: parse_neutral_losses(&value.neutral_losses)
                .map_err(BoxedError::to_owned)?,
            amino_acid_neutral_losses: value
                .amino_acid_neutral_losses
                .iter()
                .map(|rule| parse_aa_neutral_loss(rule).map_err(BoxedError::to_owned))
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

#[derive(Debug, Deserialize, PartialEq, PartialOrd, Serialize)]
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
    type Error = BoxedError<'static, BasicKind>;
    fn try_from(value: SatelliteSeriesParameters) -> Result<Self, Self::Error> {
        let selection = value
            .amino_acid_side_chain_losses_selection
            .iter()
            .map(|aa| parse_amino_acid(aa).map_err(BoxedError::to_owned))
            .collect::<Result<Vec<_>, _>>()?;
        Ok(SatelliteIonSeries {
            location: SatelliteLocation {
                rules: value
                    .location_rules
                    .into_iter()
                    .map(|rule| parse_satellite_ion(&rule).map_err(BoxedError::to_owned))
                    .collect::<Result<Vec<_>, _>>()?,
                base: value.location_base,
            },
            neutral_losses: parse_neutral_losses(&value.neutral_losses)
                .map_err(BoxedError::to_owned)?,
            amino_acid_neutral_losses: value
                .amino_acid_neutral_losses
                .iter()
                .map(|rule| parse_aa_neutral_loss(rule).map_err(BoxedError::to_owned))
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

fn parse_neutral_losses(
    neutral_losses: &[String],
) -> Result<Vec<NeutralLoss>, BoxedError<'static, BasicKind>> {
    neutral_losses
        .iter()
        .filter(|n| !n.is_empty())
        .map(|n| n.parse::<NeutralLoss>())
        .collect()
}

pub fn parameters(
    tolerance: (f64, &str),
    mz_range: (Option<f64>, Option<f64>),
) -> Result<MatchingParameters, BoxedError<'static, BasicKind>> {
    let mut parameters = MatchingParameters::default();
    if tolerance.1 == "ppm" {
        parameters.tolerance = Tolerance::new_ppm(tolerance.0);
    } else if tolerance.1 == "th" {
        parameters.tolerance =
            Tolerance::new_absolute(MassOverCharge::new::<mzcore::system::thomson>(tolerance.0));
    } else {
        return Err(BoxedError::new(
            BasicKind::Error,
            "Invalid tolerance unit",
            "",
            Context::none(),
        ));
    }
    let min = mz_range.0.unwrap_or(0.0);
    let max = mz_range.1.unwrap_or(f64::MAX);
    if min > max {
        return Err(BoxedError::new(
            BasicKind::Error,
            "m/z range invalid",
            "The minimal value is less then the maximal value",
            Context::none(),
        ));
    }
    parameters.mz_range = MassOverCharge::new::<thomson>(min)..=MassOverCharge::new::<thomson>(max);
    Ok(parameters)
}
