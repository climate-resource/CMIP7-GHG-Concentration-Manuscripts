"""
Checks of the values behind the scenario GHG manuscript's statements

See [local.value_checks][] for how these are used.
"""

from __future__ import annotations

import importlib.util
import json
import re
from functools import cache, partial
from pathlib import Path

import numpy as np
import openscm_units
import pandas as pd
import xarray as xr
import yaml

from local.historical_ghg_forcing_for_cmip7.comparison_data import (
    TIME_COLUMN,
    VALUE_COLUMN,
    get_igcc_comparison,
    get_radiative_effect_per_unit,
)
from local.scenario_ghg_forcing_for_cmip7.bundle import BUNDLE_DIR
from local.scenario_ghg_forcing_for_cmip7.esgf_data import (
    get_history_yearly_global_mean,
    get_scenario_yearly_global_mean,
)
from local.scenario_ghg_forcing_for_cmip7.magicc_jump import (
    EXPERIMENTS,
    MAGICC_DIR,
    MAGICC_PROB_DISTRIBUTION,
    get_complete_emissions_file,
)
from local.scenario_ghg_forcing_for_cmip7.magicc_jump_summary import (
    JUMP_PERIODS,
    JUMP_VARIABLES,
    MAGICC_JUMP_SUMMARY,
    load_published,
)
from local.scenario_ghg_forcing_for_cmip7.methods_data import (
    CO2_SEASONALITY_REGRESSION,
    ESGF_READY_DIRS,
    HISTORICAL_ESGF_READY_DIR,
    HISTORICAL_INTERIM_DIR,
    INTERIM_DIR,
    LAST_HISTORICAL_YEAR,
    get_max_radiative_effect_difference,
    get_projected_pcs,
    get_regression,
    load_historical_yearly_global_mean,
    load_published_yearly_global_mean,
)
from local.scenario_ghg_forcing_for_cmip7.rerun_projections import (
    MARKERS,
    SINGLE_PROJECTION_GASES,
    TWO_EOF_GASES,
    get_regression_info_file,
)
from local.scenario_ghg_forcing_for_cmip7.rerun_vl_cf import (
    VL_CF_HARMONISATION_YEAR,
    VL_CF_MARKER,
)
from local.scenario_ghg_forcing_for_cmip7.scenarios import (
    CMIP6_SCENARIOS,
    CMIP7_SCENARIOS,
    SCENARIOS,
)
from local.value_checks import CheckValue, ValueCheck

MAGICC_SUMMARY_FILE = Path("manuscript-outputs") / "magicc-summary.feather"
"""MAGICC summary, relative to the bundle's root"""


def get_single_projection_magicc_max_erf_difference_sum(
    summary_file: Path, scenario: str
) -> CheckValue:
    """
    Get the sum of the single-projection gases' largest radiative effect differences

    For each gas, the largest (over time) difference in approximate radiative effect
    between MAGICC's concentrations and ours, summed over the gases.

    Parameters
    ----------
    summary_file
        File which holds the MAGICC summary

    scenario
        Scenario of interest

    Returns
    -------
    :
        The sum
    """
    total = sum(
        get_max_radiative_effect_difference(summary_file, gas, scenario)
        for gas in SINGLE_PROJECTION_GASES
    )

    return openscm_units.unit_registry.Quantity(total, "W / m^2")


def get_single_projection_magicc_max_erf_difference_sum_range(
    summary_file: Path,
) -> CheckValue:
    """
    Get the range across all experiments of the sum of the largest differences

    All experiments means the ScenarioMIP scenarios and VL-CF.

    See [get_single_projection_magicc_max_erf_difference_sum][].

    Returns a quantity with the lowest and highest values,
    so the check passes only if both are within the statement's bounds.
    """
    values = [
        get_single_projection_magicc_max_erf_difference_sum(summary_file, name).m
        for _, _, name in (*MARKERS, VL_CF_MARKER)
    ]

    return openscm_units.unit_registry.Quantity([min(values), max(values)], "W / m^2")


IGCC_GASES = ("co2", "ch4", "n2o")
"""Gases whose results we compare with IGCC

IGCC's equivalent species are defined differently to ours, so aren't compared.
"""


def to_yearly_series(da: xr.DataArray) -> pd.Series[float]:
    """
    Get yearly data as a series, indexed by (integer) year

    Parameters
    ----------
    da
        Yearly data, with `time` as decimal years and nothing else

    Returns
    -------
    :
        `da`'s values, indexed by year
    """
    return pd.Series(
        da.values, index=np.floor(da["time"].values).astype(int), dtype=float
    )


def get_igcc_yearly(gas: str, units: str) -> pd.Series[float]:
    """
    Get IGCC's global-, annual-mean, indexed by (integer) year

    Parameters
    ----------
    gas
        Gas of interest

    units
        Units to return the data in

    Returns
    -------
    :
        IGCC's record of `gas`, in `units`
    """
    comparison = get_igcc_comparison(gas)
    factor = openscm_units.unit_registry.Quantity(1.0, comparison.units).to(units).m

    return pd.Series(
        comparison.data[VALUE_COLUMN].to_numpy() * factor,
        index=np.floor(comparison.data[TIME_COLUMN].to_numpy()).astype(int),
    )


def get_igcc_last_year() -> CheckValue:
    """
    Get the last year of IGCC's record

    The same for every gas we compare with IGCC.

    Returns
    -------
    :
        The last year
    """
    (last_year,) = {
        int(get_igcc_yearly(gas, get_igcc_comparison(gas).units).index.max())
        for gas in IGCC_GASES
    }

    return openscm_units.unit_registry.Quantity(last_year, "yr")


def get_scenariomip_minus_igcc_last_igcc_year(gas: str) -> CheckValue:
    """
    Get the range across the ScenarioMIP scenarios of their difference from IGCC

    In the last year of IGCC's record.
    VL-CF is left out: it diverges from history in 2016,
    so isn't meant to follow the observations.

    Parameters
    ----------
    gas
        Gas of interest

    Returns
    -------
    :
        The lowest and highest difference (scenario minus IGCC)
    """
    last_year = int(get_igcc_last_year().m)
    differences = []
    for scenario in CMIP7_SCENARIOS:
        if scenario.name == VL_CF_MARKER[-1]:
            continue

        da = get_scenario_yearly_global_mean(gas, scenario)
        igcc = get_igcc_yearly(gas, da.attrs["units"])
        differences.append(to_yearly_series(da)[last_year] - igcc[last_year])

    return openscm_units.unit_registry.Quantity(
        [min(differences), max(differences)], da.attrs["units"]
    )


def get_history_minus_igcc(gas: str, year: int) -> CheckValue:
    """
    Get the difference between the CMIP7 history and IGCC in a year

    Parameters
    ----------
    gas
        Gas of interest

    year
        Year of interest

    Returns
    -------
    :
        The difference (history minus IGCC)
    """
    da = get_history_yearly_global_mean(gas, "CMIP7")
    igcc = get_igcc_yearly(gas, da.attrs["units"])

    return openscm_units.unit_registry.Quantity(
        to_yearly_series(da)[year] - igcc[year], da.attrs["units"]
    )


HISTORICAL_CH4_LAT_GRADIENT_REGRESSION_FILE = (
    HISTORICAL_INTERIM_DIR / "ch4" / "ch4_pc0-ch4-fossil-emissions-regression.yaml"
)
"""The historical dataset's regression of CH4's first latitudinal gradient PC

Against PRIMAP-hist's fossil CH4 emissions.
"""


def get_ch4_lat_gradient_regression_slope_ratio_to_historical(
    bundle_dir: Path,
) -> CheckValue:
    """
    Get the ratio of our CH4 latitudinal gradient regression's slope to the historical

    Parameters
    ----------
    bundle_dir
        Directory which holds the original run's bundle (including our re-runs)

    Returns
    -------
    :
        Our slope divided by the historical dataset's slope
    """
    _, info = get_regression("ch4", bundle_dir=bundle_dir)
    ours = openscm_units.unit_registry.Quantity(info["slope"], info["slope_unit"])

    historical_info = yaml.safe_load(
        HISTORICAL_CH4_LAT_GRADIENT_REGRESSION_FILE.read_text()
    )
    historical = openscm_units.unit_registry.Quantity(*historical_info["m"])

    return (ours / historical).to("dimensionless")


def get_ch4_lat_gradient_regression_year(bundle_dir: Path, *, last: bool) -> CheckValue:
    """
    Get the first or last year our CH4 latitudinal gradient regression is fitted over

    Parameters
    ----------
    bundle_dir
        Directory which holds the original run's bundle (including our re-runs)

    last
        Whether to get the last year (otherwise the first)

    Returns
    -------
    :
        The year
    """
    data, _ = get_regression("ch4", bundle_dir=bundle_dir)
    years = data["year"]

    return openscm_units.unit_registry.Quantity(
        int(years.max() if last else years.min()), "yr"
    )


VL_CF_REGRESSION_KEYS = (CO2_SEASONALITY_REGRESSION, "co2", "ch4", "n2o")
"""Keys of the regressions which were re-calculated for VL-CF, see [get_regression][]"""


def get_regression_slope(key: str, bundle_dir: Path, last_year: int | None) -> float:
    """
    Get the slope of one of our regressions, fitted up to a given year

    Parameters
    ----------
    key
        Key of the regression, see [get_regression][]

    bundle_dir
        Directory which holds the original run's bundle (including our re-runs)

    last_year
        Last year to fit over

        If `None`, all the years the regression uses (as for the ScenarioMIP scenarios).

    Returns
    -------
    :
        Slope, in the regression's units
    """
    data, _ = get_regression(key, bundle_dir=bundle_dir)
    if last_year is not None:
        data = data[data["year"] <= last_year]

    return float(np.polyfit(data["x"], data["y"], 1)[0])


def get_vl_cf_regression_slope_ratio(key: str, bundle_dir: Path) -> CheckValue:
    """
    Get the ratio of VL-CF's regression slope to the ScenarioMIP scenarios'

    VL-CF's regressions only use the years up to its harmonisation year,
    see [get_vl_cf_regression_reproduced][].

    Parameters
    ----------
    key
        Key of the regression, see [get_regression][]

    bundle_dir
        Directory which holds the original run's bundle (including our re-runs)

    Returns
    -------
    :
        VL-CF's slope divided by the ScenarioMIP scenarios' slope
    """
    vl_cf = get_regression_slope(key, bundle_dir, VL_CF_HARMONISATION_YEAR)
    scenariomip = get_regression_slope(key, bundle_dir, None)

    return openscm_units.unit_registry.Quantity(vl_cf / scenariomip, "dimensionless")


def get_vl_cf_regression_reproduced(bundle_dir: Path) -> CheckValue:
    """
    Get whether fitting up to VL-CF's harmonisation year reproduces its regression

    We only have VL-CF's own regression result for the CO2 seasonality change PC
    (from our re-run of the VL-CF notebook).
    Reproducing it this way confirms how VL-CF's regressions were calculated,
    which [get_vl_cf_regression_slope_ratio][] relies on for the other regressions.

    Parameters
    ----------
    bundle_dir
        Directory which holds the original run's bundle (including our re-runs)

    Returns
    -------
    :
        Whether the slopes agree
    """
    info = json.loads(
        (
            bundle_dir / get_regression_info_file(f"{CO2_SEASONALITY_REGRESSION}-vl-cf")
        ).read_text()
    )
    refit = get_regression_slope(
        CO2_SEASONALITY_REGRESSION, bundle_dir, VL_CF_HARMONISATION_YEAR
    )

    return bool(np.isclose(refit, info["slope"], rtol=1e-8))


def get_regression_year(key: str, bundle_dir: Path, *, last: bool) -> CheckValue:
    """
    Get the first or last year one of our regressions is fitted over

    For the ScenarioMIP scenarios.

    Parameters
    ----------
    key
        Key of the regression, see [get_regression][]

    bundle_dir
        Directory which holds the original run's bundle (including our re-runs)

    last
        Whether to get the last year (otherwise the first)

    Returns
    -------
    :
        The year
    """
    data, _ = get_regression(key, bundle_dir=bundle_dir)
    years = data["year"]

    return openscm_units.unit_registry.Quantity(
        int(years.max() if last else years.min()), "yr"
    )


def get_pc0_in_year(gas: str, scenario: str, year: int, bundle_dir: Path) -> CheckValue:
    """
    Get the projected first latitudinal gradient PC in a year

    Parameters
    ----------
    gas
        Gas of interest

    scenario
        Scenario (or VL-CF) of interest

    year
        Year of interest

    bundle_dir
        Directory which holds the original run's bundle (including our re-runs)

    Returns
    -------
    :
        The PC's annual value
    """
    pcs = get_projected_pcs(gas, bundle_dir=bundle_dir, include_vl_cf=True)
    in_year = pcs.loc[(pcs.index >= year) & (pcs.index < year + 1), (0, scenario)]
    (value,) = in_year.to_numpy()

    return openscm_units.unit_registry.Quantity(float(value), "dimensionless")


COMPARISON_YEARS = (2022, 2500)
"""First and last year over which we compare experiments (and CMIP phases)

CMIP7's experiments start in 2022 (VL-CF's data starts earlier, but follows history),
CMIP6's in 2015, so we compare over the years both cover.
"""

ABSTRACT_GASES = ("co2", "ch4", "n2o")
"""Gases whose ranges the abstract gives"""


@cache
def get_yearly_global_means(gas: str, cmip_era: str) -> tuple[pd.DataFrame, str]:
    """
    Get every experiment's yearly global-mean, over [COMPARISON_YEARS][]

    CMIP7 comes from the bundle's copy of the published files,
    CMIP6 from ESGF.

    Parameters
    ----------
    gas
        Gas of interest

    cmip_era
        CMIP era of interest

    Returns
    -------
    :
        Yearly global-means (one column per experiment, indexed by integer year)
        and their units
    """
    scenarios = CMIP7_SCENARIOS if cmip_era == "CMIP7" else CMIP6_SCENARIOS
    columns = {}
    units = set()
    for scenario in scenarios:
        if cmip_era == "CMIP7":
            da = load_published_yearly_global_mean(gas, scenario.name)
        else:
            da = get_scenario_yearly_global_mean(gas, scenario)

        columns[scenario.name] = to_yearly_series(da)
        units.add(da.attrs["units"])

    (unit,) = units
    res = pd.DataFrame(columns).loc[COMPARISON_YEARS[0] : COMPARISON_YEARS[1]]

    return res, unit


def get_extreme(gas: str, cmip_era: str, kind: str) -> tuple[float, int, str, str]:
    """
    Get the lowest or highest yearly global-mean across experiments and years

    Parameters
    ----------
    gas
        Gas of interest

    cmip_era
        CMIP era of interest

    kind
        "min" or "max"

    Returns
    -------
    :
        Value, year, experiment and units
    """
    df, unit = get_yearly_global_means(gas, cmip_era)
    stacked = df.stack()
    year, experiment = stacked.idxmin() if kind == "min" else stacked.idxmax()

    return float(stacked[(year, experiment)]), int(year), str(experiment), unit


def get_extreme_value(gas: str, cmip_era: str, kind: str) -> CheckValue:
    """Get the value of [get_extreme][]"""
    value, _, _, unit = get_extreme(gas, cmip_era, kind)

    return openscm_units.unit_registry.Quantity(value, unit)


def get_extreme_year(gas: str, cmip_era: str, kind: str) -> CheckValue:
    """Get the year of [get_extreme][]"""
    _, year, _, _ = get_extreme(gas, cmip_era, kind)

    return openscm_units.unit_registry.Quantity(year, "yr")


def get_extreme_is_in(
    gas: str, cmip_era: str, kind: str, experiment: str
) -> CheckValue:
    """Get whether [get_extreme][] is in a given experiment"""
    _, _, extreme_experiment, _ = get_extreme(gas, cmip_era, kind)

    return extreme_experiment == experiment


def get_range_width_relative_to_cmip6(gas: str) -> CheckValue:
    """
    Get the width of CMIP7's range (over experiments and years) relative to CMIP6's

    Parameters
    ----------
    gas
        Gas of interest

    Returns
    -------
    :
        CMIP7's width divided by CMIP6's
    """
    widths = []
    for cmip_era in ("CMIP7", "CMIP6"):
        lowest, _, _, unit = get_extreme(gas, cmip_era, "min")
        highest, _, _, _ = get_extreme(gas, cmip_era, "max")
        widths.append(
            openscm_units.unit_registry.Quantity(highest - lowest, unit).to("ppt")
        )

    return (widths[0] / widths[1]).to("dimensionless")


def get_cfc12_cmip7_above_all_cmip6() -> CheckValue:
    """Get whether CMIP7's CFC-12 is above every CMIP6 scenario's in every year"""
    cmip7, unit_cmip7 = get_yearly_global_means("cfc12", "CMIP7")
    cmip6, unit_cmip6 = get_yearly_global_means("cfc12", "CMIP6")
    if unit_cmip7 != unit_cmip6:
        raise NotImplementedError((unit_cmip7, unit_cmip6))

    return bool((cmip6.lt(cmip7.iloc[:, 0], axis="index")).all().all())


def get_cfc12_relative_to_cmip6_mean(last_year: int) -> CheckValue:
    """
    Get how much higher CMIP7's CFC-12 is than the mean of the CMIP6 scenarios

    Parameters
    ----------
    last_year
        Last year to average over (from the start of [COMPARISON_YEARS][])

    Returns
    -------
    :
        Mean over the years of CMIP7 / mean of CMIP6 - 1
    """
    cmip7, _ = get_yearly_global_means("cfc12", "CMIP7")
    cmip6, _ = get_yearly_global_means("cfc12", "CMIP6")
    # The same in every CMIP7 experiment
    ratio = cmip7.iloc[:, 0] / cmip6.mean(axis="columns") - 1.0

    return openscm_units.unit_registry.Quantity(
        float(ratio.loc[:last_year].mean()), "dimensionless"
    )


def get_difference_radiative_effect(
    gas: str, minuend: str, subtrahend: str, year: int
) -> CheckValue:
    """
    Get the approximate radiative effect of the difference between two experiments

    Parameters
    ----------
    gas
        Gas of interest

    minuend
        Experiment to subtract from (CMIP6 or CMIP7)

    subtrahend
        Experiment to subtract (CMIP6 or CMIP7)

    year
        Year of interest

    Returns
    -------
    :
        Approximate radiative effect of `minuend` minus `subtrahend`
    """
    values = []
    for name in (minuend, subtrahend):
        cmip_era = SCENARIOS[name].cmip_era
        df, unit = get_yearly_global_means(gas, cmip_era)
        values.append(df.loc[year, name] * get_radiative_effect_per_unit(gas, unit))

    return openscm_units.unit_registry.Quantity(values[0] - values[1], "W / m^2")


def get_range_in_year(gas: str, experiments: tuple[str, ...], year: int) -> CheckValue:
    """
    Get the lowest and highest CMIP7 yearly global-mean across experiments in a year

    Parameters
    ----------
    gas
        Gas of interest

    experiments
        CMIP7 experiments to consider

    year
        Year of interest

    Returns
    -------
    :
        The lowest and highest value
    """
    df, unit = get_yearly_global_means(gas, "CMIP7")
    values = df.loc[year, list(experiments)]

    return openscm_units.unit_registry.Quantity(
        [float(values.min()), float(values.max())], unit
    )


def get_vl_minus_vl_cf(gas: str) -> tuple[pd.Series[float], str]:
    """
    Get VL's yearly global-mean minus VL-CF's

    Parameters
    ----------
    gas
        Gas of interest

    Returns
    -------
    :
        Difference, indexed by year, and its units
    """
    df, unit = get_yearly_global_means(gas, "CMIP7")

    return df["vl"] - df[VL_CF_MARKER[-1]], unit


def get_vl_minus_vl_cf_max(gas: str) -> CheckValue:
    """Get the largest difference in [get_vl_minus_vl_cf][]"""
    diff, unit = get_vl_minus_vl_cf(gas)

    return openscm_units.unit_registry.Quantity(float(diff.max()), unit)


def get_vl_minus_vl_cf_max_year(gas: str) -> CheckValue:
    """Get the year of the largest difference in [get_vl_minus_vl_cf][]"""
    diff, _ = get_vl_minus_vl_cf(gas)

    return openscm_units.unit_registry.Quantity(int(diff.idxmax()), "yr")


def get_single_projection_final_year_extrapolated() -> CheckValue:
    """
    Get whether every single projection gas' last year is linearly extrapolated

    I.e. whether the value in 2500 is the linear extrapolation
    of the values in 2498 and 2499.
    The same in every experiment, so we check L.
    """
    res = True
    for gas in SINGLE_PROJECTION_GASES:
        values = to_yearly_series(load_published_yearly_global_mean(gas, "l"))
        extrapolated = 2 * values[2499] - values[2498]
        res = res and bool(np.isclose(values[2500], extrapolated, rtol=1e-5, atol=1e-6))

    return res


WESTERN_ET_AL_2024_GASES = ("hcfc22", "hcfc141b", "hcfc142b")
"""Gases which use the Western et al. (2024) projections"""


def get_western_et_al_2024_2100_fraction_of_peak() -> CheckValue:
    """
    Get the largest fraction of its peak a Western et al. (2024) gas is at in 2100

    The peak is over the historical dataset and the projection.

    Returns
    -------
    :
        Largest (over the gases) concentration in 2100 divided by the gas' peak
    """
    fractions = []
    for gas in WESTERN_ET_AL_2024_GASES:
        history = to_yearly_series(get_history_yearly_global_mean(gas, "CMIP7"))
        projection = to_yearly_series(load_published_yearly_global_mean(gas, "l"))
        peak = max(history.max(), projection.max())
        fractions.append(projection[2100] / peak)

    return openscm_units.unit_registry.Quantity(max(fractions), "dimensionless")


def get_western_et_al_2024_2100_radiative_effect() -> CheckValue:
    """
    Get the approximate radiative effect of the Western et al. (2024) gases in 2100

    Summed over the gases, relative to their (zero) pre-industrial concentrations.
    """
    total = 0.0
    for gas in WESTERN_ET_AL_2024_GASES:
        da = load_published_yearly_global_mean(gas, "l")
        total += to_yearly_series(da)[2100] * get_radiative_effect_per_unit(
            gas, da.attrs["units"]
        )

    return openscm_units.unit_registry.Quantity(total, "W / m^2")


def load_native_resolution_in_year(
    directory: Path, gas: str, year: int
) -> xr.DataArray:
    """
    Load a published monthly native resolution file's values in a year

    Parameters
    ----------
    directory
        Directory of the dataset (i.e. of its source ID)

    gas
        Gas of interest

    year
        Year of interest

    Returns
    -------
    :
        Values in `year` (time and latitude)
    """
    for file in sorted((directory / "atmos" / "mon" / gas / "gnz").glob("*/*.nc")):
        da = xr.load_dataset(
            file, decode_times=xr.coders.CFDatetimeCoder(use_cftime=True)
        )[gas]
        in_year = da.sel(time=da["time"].dt.year == year)
        if in_year.size:
            return in_year

    msg = f"No {gas} data in {year} in {directory}"
    raise FileNotFoundError(msg)


def get_max_monthly_difference_from_history_radiative_effect(
    bundle_dir: Path, year: int = 2022
) -> CheckValue:
    """
    Get the largest difference from the historical dataset in its last year

    At native resolution, over every gas, ScenarioMIP scenario, latitude and month,
    in approximate radiative effect terms.

    Parameters
    ----------
    bundle_dir
        Directory which holds the original run's bundle

    year
        Year of interest (the last year of the historical dataset)

    Returns
    -------
    :
        The largest absolute difference
    """
    largest = 0.0
    for gas_dir in sorted((HISTORICAL_ESGF_READY_DIR / "atmos" / "mon").iterdir()):
        gas = gas_dir.name
        history = load_native_resolution_in_year(HISTORICAL_ESGF_READY_DIR, gas, year)
        per_unit = get_radiative_effect_per_unit(gas, history.attrs["units"])
        for scenario in CMIP7_SCENARIOS:
            if scenario.name == VL_CF_MARKER[-1]:
                continue

            projection = load_native_resolution_in_year(
                bundle_dir / ESGF_READY_DIRS[0] / scenario.source_ids[0], gas, year
            )
            if projection.attrs["units"] != history.attrs["units"]:
                raise NotImplementedError(gas)

            diff = np.abs(projection.values - history.values).max() * per_unit
            largest = max(largest, float(diff))

    return openscm_units.unit_registry.Quantity(largest, "W / m^2")


def load_magicc_jump_summary(bundle_dir: Path) -> dict:
    """
    Load the summary of the MAGICC runs with and without the CMIP7 history

    See [local.scenario_ghg_forcing_for_cmip7.magicc_jump_summary][].
    """
    return json.loads((bundle_dir / MAGICC_JUMP_SUMMARY).read_text())


def get_magicc_jump_reproduces_published(bundle_dir: Path) -> CheckValue:
    """Get whether our `published` MAGICC runs reproduce the published output exactly"""
    summary = load_magicc_jump_summary(bundle_dir)

    return all(
        diff == 0.0
        for by_variable in summary["reproduction"].values()
        for diff in by_variable.values()
    )


def get_magicc_jump_n_runs(bundle_dir: Path) -> CheckValue:
    """Get the number of ensemble members in the MAGICC runs"""
    return openscm_units.unit_registry.Quantity(
        load_magicc_jump_summary(bundle_dir)["n_runs"], "dimensionless"
    )


def get_magicc_jump_max_difference(
    short_name: str, period_name: str, bundle_dir: Path
) -> CheckValue:
    """
    Get the largest difference between the MAGICC runs with and without the jump

    Parameters
    ----------
    short_name
        Short name of the variable, see
        [local.scenario_ghg_forcing_for_cmip7.magicc_jump_summary.JUMP_VARIABLES][]

    period_name
        Name of the period, see
        [local.scenario_ghg_forcing_for_cmip7.magicc_jump_summary.JUMP_PERIODS][]

    bundle_dir
        Directory which holds the bundle

    Returns
    -------
    :
        Largest absolute difference between the ensemble medians
    """
    units = {"gsat": "K", "co2": "ppm", "ch4": "ppb", "n2o": "ppb"}.get(
        short_name, "W / m^2"
    )
    value = load_magicc_jump_summary(bundle_dir)["max_median_difference"][short_name][
        period_name
    ]

    return openscm_units.unit_registry.Quantity(value["value"], units)


SCENARIOMIP_EXPERIMENTS = tuple(name for _, _, name in MARKERS)
"""The ScenarioMIP experiments (i.e. every experiment except VL-CF)"""

HARMONISATION_YEARS = {
    **dict.fromkeys(SCENARIOMIP_EXPERIMENTS, LAST_HISTORICAL_YEAR),
    VL_CF_MARKER[-1]: VL_CF_HARMONISATION_YEAR,
}
"""Year in which each experiment is harmonised to (and leaves) history"""


def to_check_value(values: list[float] | list[int], units: str) -> CheckValue:
    """
    Convert values which should all be the same to a check value

    If they are all the same, a single quantity.
    Otherwise, a quantity with the lowest and highest values,
    so the check passes only if both are within the statement's bounds.
    """
    if len(set(values)) == 1:
        return openscm_units.unit_registry.Quantity(values[0], units)

    return openscm_units.unit_registry.Quantity([min(values), max(values)], units)


def get_published_source_id_dirs(bundle_dir: Path) -> dict[str, Path]:
    """
    Get the directory of every published source ID

    Parameters
    ----------
    bundle_dir
        Directory which holds the original run's bundle

    Returns
    -------
    :
        Directory of each source ID, keyed by source ID
    """
    return {
        source_id_dir.name: source_id_dir
        for esgf_ready_dir in ESGF_READY_DIRS
        for source_id_dir in sorted((bundle_dir / esgf_ready_dir).iterdir())
    }


def get_published_variables(bundle_dir: Path) -> tuple[str, ...]:
    """
    Get the variables (gases and equivalent species) we publish

    Taken from L (every experiment has the same, see [get_n_output_variables][]).
    """
    (l_dir,) = (
        source_id_dir
        for source_id, source_id_dir in get_published_source_id_dirs(bundle_dir).items()
        if source_id == SCENARIOS["l"].source_ids[0]
    )

    return tuple(sorted(p.name for p in (l_dir / "atmos" / "mon").iterdir()))


def get_native_grid_lat_band_width(bundle_dir: Path) -> CheckValue:
    """
    Get the width of the latitude bands of our native grid

    From the latitude bounds of every published source ID's CO2 native grid files.
    """
    widths = set()
    for source_id_dir in get_published_source_id_dirs(bundle_dir).values():
        for file in (source_id_dir / "atmos" / "mon" / "co2" / "gnz").glob("*/*.nc"):
            lat_bnds = xr.open_dataset(file, decode_times=False)["lat_bnds"].values
            widths.update(np.unique(lat_bnds[:, 1] - lat_bnds[:, 0]).tolist())

    return to_check_value(sorted(widths), "degree")


def get_file_time_range(file: Path) -> tuple[str, str]:
    """
    Get the time range of a published file from its name

    E.g. `("2022", "2100")` for yearly files
    and `("202201", "210012")` for monthly files.
    """
    return tuple(file.stem.rsplit("_", 1)[1].split("-"))  # type: ignore[return-value]


def get_output_year(
    experiments: tuple[str, ...], bundle_dir: Path, *, last: bool
) -> CheckValue:
    """
    Get the first or last year of our published data

    From the names of the monthly files, so it includes every month
    (e.g. a last year of 2500 means the data runs to December 2500).

    Parameters
    ----------
    experiments
        Experiments of interest

    bundle_dir
        Directory which holds the original run's bundle

    last
        Whether to get the last year (otherwise the first)

    Returns
    -------
    :
        The year (if not the same in all the experiments,
        the range across them, see [to_check_value][])
    """
    source_id_dirs = get_published_source_id_dirs(bundle_dir)
    years = []
    for experiment in experiments:
        months = [
            month
            for source_id in SCENARIOS[experiment].source_ids
            for file in (source_id_dirs[source_id] / "atmos" / "mon").glob("*/*/*/*.nc")
            for month in get_file_time_range(file)
        ]
        month = max(months) if last else min(months)
        if month.endswith("12" if last else "01"):
            years.append(int(month[:4]))
        else:
            msg = f"{experiment} doesn't cover whole years ({month=})"
            raise AssertionError(msg)

    return to_check_value(years, "yr")


def get_n_output_variables(bundle_dir: Path) -> CheckValue:
    """
    Get the number of variables we publish (gases and equivalent species)

    Counted for every published source ID, frequency and grid.
    """
    counts = []
    for source_id_dir in get_published_source_id_dirs(bundle_dir).values():
        for frequency_dir in (source_id_dir / "atmos").iterdir():
            variables = {p.parent.name for p in frequency_dir.glob("*/*")}
            for grid in {p.name for p in frequency_dir.glob("*/*")}:
                counts.append(
                    sum((frequency_dir / v / grid).is_dir() for v in variables)
                )

    return to_check_value(counts, "dimensionless")


ESGF_VERSION = "1-1-0"
"""Version of our published data, as it appears in the source IDs"""


def get_source_ids_as_described(bundle_dir: Path) -> CheckValue:
    """
    Get whether the published source IDs have the form we describe

    `CR-[experiment]-[version]` up to 2100
    and `CR-[experiment]-ext-[version]` for the extensions,
    with the experiment in lowercase.
    """
    expected = {
        source_id
        for scenario in CMIP7_SCENARIOS
        for source_id in (
            f"CR-{scenario.name.lower()}-{ESGF_VERSION}",
            f"CR-{scenario.name.lower()}-ext-{ESGF_VERSION}",
        )
    }
    known = {
        source_id for scenario in CMIP7_SCENARIOS for source_id in scenario.source_ids
    }

    return set(get_published_source_id_dirs(bundle_dir)) == expected == known


EXTENSION_FILE_YEARS = ((2101, 2200), (2201, 2300), (2301, 2400), (2401, 2500))
"""Years covered by each of the extensions' files"""


def get_output_time_chunks_as_described(bundle_dir: Path) -> CheckValue:
    """
    Get whether the published files are split in time as we describe

    Up to 2100, a single file per variable, frequency and grid
    (from 2022 for the ScenarioMIP experiments, 2016 for VL-CF).
    The extensions (2101-2500), one file per 100 years.
    """
    source_id_dirs = get_published_source_id_dirs(bundle_dir)
    for scenario in CMIP7_SCENARIOS:
        main, extension = scenario.source_ids
        first_year = HARMONISATION_YEARS[scenario.name] + (
            1 if scenario.name == VL_CF_MARKER[-1] else 0
        )
        expected_years = {
            main: {(first_year, 2100)},
            extension: set(EXTENSION_FILE_YEARS),
        }
        for source_id, years in expected_years.items():
            for frequency_dir in (source_id_dirs[source_id] / "atmos").iterdir():
                for grid_dir in frequency_dir.glob("*/*"):
                    ranges = {
                        get_file_time_range(file) for file in grid_dir.glob("*/*.nc")
                    }
                    if frequency_dir.name == "mon":
                        expected = {(f"{s}01", f"{e}12") for s, e in years}
                    else:
                        expected = {(str(s), str(e)) for s, e in years}

                    if ranges != expected:
                        return False

    return True


def get_magicc_n_ensemble_members(bundle_dir: Path) -> CheckValue:
    """
    Get the number of ensemble members in the published MAGICC output

    For every experiment.
    """
    return to_check_value(
        [
            int(
                load_published(name, ("Surface Air Temperature Change",), bundle_dir)
                .index.get_level_values("run_id")
                .nunique()
            )
            for _, _, name in EXPERIMENTS
        ],
        "dimensionless",
    )


def get_magicc_switch_from_conc_to_emissions_year(bundle_dir: Path) -> CheckValue:
    """
    Get the year in which MAGICC switches from concentrations to emissions

    From MAGICC's default configuration
    (the bundle's MAGICC, extracted by the MAGICC jump runs),
    after checking that the probabilistic distribution doesn't override it.

    Returns
    -------
    :
        The year (if not the same for every group of gases,
        the range across them, see [to_check_value][])
    """
    cfg = (bundle_dir / MAGICC_DIR / "run" / "MAGCFG_DEFAULTALL.CFG").read_text()
    years = [
        int(year)
        for year in re.findall(
            r"^\s*\w+_SWITCHFROMCONC2EMIS_YEAR\s*=\s*(\d+)", cfg, flags=re.M | re.I
        )
    ]
    if not years:
        raise AssertionError(cfg)

    prob_distribution = json.loads((bundle_dir / MAGICC_PROB_DISTRIBUTION).read_text())
    overridden = {
        key
        for configuration in prob_distribution["configurations"]
        for key in configuration["nml_allcfgs"]
        if "switchfromconc2emis" in key.lower()
    }
    if overridden:
        raise AssertionError(overridden)

    return to_check_value(years, "yr")


CONVERGED_RTOL = 5e-7
"""Relative difference below which our data equals MAGICC's median

Above the numerical noise in the published files (up to around 6e-8),
below the smallest difference the year before convergence (around 2e-6).
"""


def get_harmonisation_convergence_time(gas: str, bundle_dir: Path) -> CheckValue:
    """
    Get the convergence time of the harmonisation of MAGICC's output

    From the data:
    the time from the harmonisation year
    until our global-, annual-mean equals MAGICC's median (from then on).

    Parameters
    ----------
    gas
        Gas of interest

    bundle_dir
        Directory which holds the original run's bundle

    Returns
    -------
    :
        The convergence time (if not the same in every experiment,
        the range across them, see [to_check_value][])
    """
    times = []
    for _, _, name in EXPERIMENTS:
        magicc = load_published(
            name, (f"Atmospheric Concentrations|{gas.upper()}",), bundle_dir
        ).median(axis="index")
        ours = to_yearly_series(
            load_published_yearly_global_mean(gas, name, bundle_dir)
        )
        years = ours.index.intersection(magicc.index)
        relative_difference = (ours[years] - magicc[years].to_numpy()).abs() / ours[
            years
        ].abs()
        not_converged = relative_difference[relative_difference > CONVERGED_RTOL]
        times.append(int(not_converged.index.max()) + 1 - HARMONISATION_YEARS[name])

    return to_check_value(times, "yr")


def get_scenariomip_max_relative_difference_monthly_global_mean(
    year: int, bundle_dir: Path
) -> CheckValue:
    """
    Get the largest difference between the ScenarioMIP experiments' monthly global-means

    Relative to the first experiment,
    over every published variable and month in `year`.
    """
    source_id_dirs = get_published_source_id_dirs(bundle_dir)
    largest = 0.0
    for variable in get_published_variables(bundle_dir):
        values = []
        for experiment in SCENARIOMIP_EXPERIMENTS:
            source_id_dir = source_id_dirs[SCENARIOS[experiment].source_ids[0]]
            (file,) = (source_id_dir / "atmos" / "mon" / variable / "gm").glob("*/*.nc")
            da = xr.load_dataset(
                file, decode_times=xr.coders.CFDatetimeCoder(use_cftime=True)
            )[variable]
            values.append(da.sel(time=da["time"].dt.year == year).to_numpy())

        stacked = np.stack(values)
        largest = max(
            largest,
            float((np.abs(stacked - stacked[0]) / np.abs(stacked[0])).max()),
        )

    return openscm_units.unit_registry.Quantity(largest, "dimensionless")


ANNUAL_MEANS_DIR = INTERIM_DIR / "annual-means"
"""The original run's global-, annual-means, relative to the bundle's root"""

PROJECTION_METHODS = (
    "gradient-aware-harmonisation",
    "one-box",
    "single-concentration-projection",
)
"""The ways we project global-, annual-means (and the prefix of their interim files)"""

PROJECTION_METHOD_RTOL = 1e-5
"""Difference (relative to the largest value) below which our data matches a method's"""


@cache
def get_projection_method(gas: str, bundle_dir: Path) -> str:
    """
    Get the method used to project a gas' global-, annual-mean

    From the data: the method whose interim global-, annual-mean
    matches our published data.
    We compare from 2023,
    because in the harmonisation year our data is the historical dataset's.
    The interim files hold different experiments for different methods
    (VL-CF for MAGICC plus harmonisation, see `rerun_projections`),
    so we compare the experiment the file holds
    (L for the single concentration projections, which are the same in all experiments).

    Parameters
    ----------
    gas
        Gas of interest

    bundle_dir
        Directory which holds the original run's bundle

    Returns
    -------
    :
        The method, one of [PROJECTION_METHODS][]
    """
    matches = []
    for method in PROJECTION_METHODS:
        file = bundle_dir / ANNUAL_MEANS_DIR / f"{method}_{gas}_annual-mean.feather"
        if not file.exists():
            continue

        interim_df = pd.read_feather(file)
        experiment = interim_df.index.get_level_values("scenario")[0]
        interim = interim_df.iloc[0]
        interim.index = interim.index.astype(int)
        published = to_yearly_series(
            load_published_yearly_global_mean(
                gas, "l" if experiment == "all" else experiment, bundle_dir
            )
        )
        years = published.index.intersection(interim.index)
        years = years[years > LAST_HISTORICAL_YEAR]
        difference = np.abs(published[years].to_numpy() - interim[years].to_numpy())
        if difference.max() / np.abs(interim[years].to_numpy()).max() < (
            PROJECTION_METHOD_RTOL
        ):
            matches.append(method)

    if len(matches) != 1:
        raise AssertionError((gas, matches))

    return matches[0]


def get_gases_with_projection_method(method: str, bundle_dir: Path) -> set[str]:
    """
    Get the published gases whose global-, annual-mean uses a projection method

    See [get_projection_method][].
    Equivalent species are left out (they are calculated from the other gases).
    """
    return {
        gas
        for gas in get_published_variables(bundle_dir)
        if not gas.endswith("eq") and get_projection_method(gas, bundle_dir) == method
    }


CFC12_LIKE_GASES_WMO_2022 = (
    "ccl4",
    "cfc11",
    "cfc113",
    "cfc114",
    "cfc115",
    "ch3br",
    "ch3ccl3",
    "ch3cl",
    "halon1211",
    "halon1301",
    "halon2402",
)
"""Gases processed like CFC-12 with projections from WMO (2022), as listed"""

CFC12_LIKE_GASES_WESTERN_ET_AL_2024 = ("hcfc22", "hcfc141b", "hcfc142b")
"""Gases processed like CFC-12 with projections from Western et al. (2024), as listed"""

HFC125_LIKE_GASES = (
    "c2f6",
    "c3f8",
    "c4f10",
    "c5f12",
    "c6f14",
    "c7f16",
    "c8f18",
    "cc4f8",
    "cf4",
    "ch2cl2",
    "chcl3",
    "hfc143a",
    "hfc227ea",
    "hfc236fa",
    "hfc4310mee",
    "nf3",
    "sf6",
    "so2f2",
)
"""Gases processed like HFC-125, as listed in the text"""


def get_cfc12_like_gases_as_described(bundle_dir: Path) -> CheckValue:
    """
    Get whether the gases processed like CFC-12 are those we list

    From the data (see [get_projection_method][]),
    with the Western et al. (2024) gases being those
    whose projections were extended (`data/interim/western-et-al-2024`).
    """
    single_projection = get_gases_with_projection_method(
        "single-concentration-projection", bundle_dir
    )
    western = {
        file.stem.split("_")[1]
        for file in (bundle_dir / INTERIM_DIR / "western-et-al-2024").glob(
            "western-et-al-2024_*_extended.feather"
        )
    }

    return single_projection == {
        "cfc12",
        *CFC12_LIKE_GASES_WMO_2022,
        *CFC12_LIKE_GASES_WESTERN_ET_AL_2024,
    } and western == set(CFC12_LIKE_GASES_WESTERN_ET_AL_2024)


def get_hfc125_like_gases_as_described(bundle_dir: Path) -> CheckValue:
    """
    Get whether the gases processed like HFC-125 are those we list

    From the data (see [get_projection_method][]).
    """
    return get_gases_with_projection_method("one-box", bundle_dir) == {
        "hfc125",
        *HFC125_LIKE_GASES,
    }


def get_single_projection_max_relative_difference_from_history(
    year: int, bundle_dir: Path
) -> CheckValue:
    """
    Get the largest difference between the single concentration projections and history

    Relative to the historical dataset,
    over the gases processed like CFC-12 (in our output), in `year`.
    The projections as used, i.e. before they are joined to the history
    (from the interim files).
    """
    largest = 0.0
    for gas in get_gases_with_projection_method(
        "single-concentration-projection", bundle_dir
    ):
        projection_df = pd.read_feather(
            bundle_dir
            / ANNUAL_MEANS_DIR
            / f"single-concentration-projection_{gas}_annual-mean.feather"
        )
        projection = projection_df.iloc[0]
        projection.index = projection.index.astype(int)
        history_da = load_historical_yearly_global_mean(gas)
        if history_da.attrs["units"] != projection_df.index.get_level_values("unit")[0]:
            raise NotImplementedError(gas)

        history = to_yearly_series(history_da)
        largest = max(
            largest, float(abs(projection[year] - history[year]) / history[year])
        )

    return openscm_units.unit_registry.Quantity(largest, "dimensionless")


WMO_2022_TABLE_A5_LIFETIMES = {
    "c2f6": (10000.0, 10000.0),
    "c3f8": (2600.0, 2600.0),
    "c4f10": (2600.0, 2600.0),
    "c5f12": (4100.0, 4100.0),
    "c6f14": (3100.0, 3100.0),
    "c7f16": (3000.0, 3000.0),
    "c8f18": (3000.0, 3000.0),
    "cc4f8": (3200.0, 3200.0),
    "cf4": (50000.0, 50000.0),
    "ch2cl2": (176.0 / 365.0, 176.0 / 365.0),
    "chcl3": (178.0 / 365.0, 178.0 / 365.0),
    "hcfc141b": (8.81, 8.81),
    "hcfc142b": (17.1, 17.1),
    "hcfc22": (11.6, 11.6),
    "hfc125": (30.7, 30.7),
    "hfc143a": (51.8, 51.8),
    "hfc227ea": (35.8, 35.8),
    "hfc236fa": (213.0, 213.0),
    "hfc4310mee": (17.0, 17.0),
    "nf3": (569.0, 569.0),
    "sf6": (850.0, 1280.0),
    "so2f2": (36.0, 36.0),
}
"""Total lifetimes (years) from Table A-5 of WMO (2022), as (low, high)

Only for the gases whose lifetimes we use (the one-box gases
and the gases whose Western et al. (2024) projections we extend).
Copied from the 'WMO (2022) Total Lifetime' column of
https://csl.noaa.gov/assessments/ozone/2022/downloads/Annex_2022OzoneAssessment.pdf
(lifetimes given in days converted with 365 days per year, as in the original run).
The table gives a range for SF6.
"""

ORIGINAL_RUN_CONSTANTS = (
    Path("src") / "cmip7_scenariomip_ghg_generation" / "constants.py"
)
"""The original run's constants (including lifetimes), relative to the bundle's root"""


def load_original_run_lifetimes(bundle_dir: Path) -> dict[str, float]:
    """
    Load the lifetimes used by the original run, in years

    From the bundle's copy of the original run's code
    (which only needs `openscm_units`, so we can load it here).
    """
    spec = importlib.util.spec_from_file_location(
        "original_run_constants", bundle_dir / ORIGINAL_RUN_CONSTANTS
    )
    if spec is None or spec.loader is None:
        raise AssertionError(spec)

    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    return {gas: float(v.to("yr").m) for gas, v in module.GHG_LIFETIMES.items()}


def get_lifetimes_from_wmo_2022(bundle_dir: Path) -> CheckValue:
    """
    Get whether the lifetimes we use are those in Table A-5 of WMO (2022)

    For the one-box gases
    and the gases whose Western et al. (2024) projections we extend
    (where the table gives a range, the original run uses the middle of it).
    """
    gases = {
        *get_gases_with_projection_method("one-box", bundle_dir),
        *CFC12_LIKE_GASES_WESTERN_ET_AL_2024,
    }
    if gases != set(WMO_2022_TABLE_A5_LIFETIMES):
        raise AssertionError(gases ^ set(WMO_2022_TABLE_A5_LIFETIMES))

    used = load_original_run_lifetimes(bundle_dir)

    return all(
        np.isclose(used[gas], (low + high) / 2.0, rtol=1e-10)
        for gas, (low, high) in WMO_2022_TABLE_A5_LIFETIMES.items()
    )


def get_historical_radiative_effect(gas: str, year: int) -> float:
    """
    Get a gas' approximate radiative effect in the historical dataset

    Relative to 1750, in W / m^2.
    """
    da = load_historical_yearly_global_mean(gas)
    values = to_yearly_series(da)

    return float(
        (values[year] - values[1750])
        * get_radiative_effect_per_unit(gas, da.attrs["units"])
    )


def get_hfc125_largest_radiative_effect_in_group(
    year: int, bundle_dir: Path
) -> CheckValue:
    """
    Get whether HFC-125 has the largest radiative effect of the one-box gases

    In the historical dataset, see [get_historical_radiative_effect][].
    """
    effects = {
        gas: get_historical_radiative_effect(gas, year)
        for gas in get_gases_with_projection_method("one-box", bundle_dir)
    }

    return max(effects, key=effects.__getitem__) == "hfc125"


def get_hfc125_radiative_effect(year: int) -> CheckValue:
    """Get HFC-125's radiative effect, see [get_historical_radiative_effect][]"""
    return openscm_units.unit_registry.Quantity(
        get_historical_radiative_effect("hfc125", year), "W / m^2"
    )


HISTORICAL_PCS_FILE_PATTERN = "{gas}_allyears-lat-gradient-eofs-pcs.nc"
"""Historical dataset's latitudinal gradient PCs (in its interim data)"""


@cache
def load_historical_pc(gas: str, eof: int) -> pd.Series[float]:
    """
    Load one of the historical dataset's annual latitudinal gradient PCs

    Indexed by year.
    """
    (file,) = HISTORICAL_INTERIM_DIR.rglob(HISTORICAL_PCS_FILE_PATTERN.format(gas=gas))
    pc = xr.load_dataset(file)["principal-components"].sel(eof=eof)

    return pd.Series(pc.to_numpy(), index=pc["year"].to_numpy().astype(int))


def get_second_pc_follows_history_then_constant(bundle_dir: Path) -> CheckValue:
    """
    Get whether the projected second PCs are as described

    For CO2, CH4 and N2O, in every experiment:
    the historical dataset's PC up to its last year (2023),
    then constant at its value in the experiment's harmonisation year.
    """
    for gas in TWO_EOF_GASES:
        history = load_historical_pc(gas, 1)
        pcs = get_projected_pcs(gas, bundle_dir=bundle_dir, include_vl_cf=True)[1]
        for experiment in pcs.columns:
            projected = pcs[experiment].dropna()
            projected.index = np.floor(projected.index).astype(int)
            expected = pd.Series(
                [
                    history[year]
                    if year <= history.index.max()
                    else history[HARMONISATION_YEARS[experiment]]
                    for year in projected.index
                ],
                index=projected.index,
            )
            if not np.allclose(projected, expected, rtol=1e-5, atol=1e-5):
                return False

    return True


def get_n_constant_pc0_regression_years(gas: str, bundle_dir: Path) -> CheckValue:
    """
    Get the number of years in a PC regression in which the first PC is constant

    I.e. years in which the historical dataset's first PC
    has the same value as in the year before
    (so isn't derived from observations or emissions).
    """
    data, _ = get_regression(gas, bundle_dir=bundle_dir)
    history = load_historical_pc(gas, 0)
    constant = history[history.diff() == 0.0].index

    return openscm_units.unit_registry.Quantity(
        int(np.isin(data["year"].astype(int), constant).sum()), "dimensionless"
    )


@cache
def load_complete_emissions(name: str, bundle_dir: Path) -> pd.DataFrame:
    """
    Load the complete emissions an experiment's MAGICC run used (1750 to 2500)

    Indexed by variable and unit, one column per year.
    """
    ((model, scenario),) = (
        (model, scenario) for model, scenario, n in EXPERIMENTS if n == name
    )
    res = pd.read_feather(bundle_dir / get_complete_emissions_file(model, scenario))
    res.columns = res.columns.astype(int)

    return res.reset_index(
        [n for n in res.index.names if n not in ("variable", "unit")], drop=True
    )


EMISSIONS_RTOL = 1e-4
"""Relative difference above which emissions differ (well above numerical noise)"""


def get_last_year_emissions_identical(
    names: tuple[str, ...], bundle_dir: Path
) -> CheckValue:
    """
    Get the last year in which experiments' emissions are all the same

    Before this year, all the experiments follow history
    (they are harmonised to it), so this is the last year they follow history.

    Parameters
    ----------
    names
        Experiments to compare

    bundle_dir
        Directory which holds the original run's bundle

    Returns
    -------
    :
        The last year, over every species, before which the emissions are all the same
    """
    reference = load_complete_emissions(names[0], bundle_dir)
    differs = pd.Series(False, index=reference.columns)
    for name in names[1:]:
        emissions = load_complete_emissions(name, bundle_dir).loc[reference.index]
        relative_difference = (emissions - reference).abs() / reference.abs().where(
            reference != 0.0
        )
        differs |= (relative_difference > EMISSIONS_RTOL).any(axis="index")

    return openscm_units.unit_registry.Quantity(int(differs.idxmax()) - 1, "yr")


def get_value_checks(bundle_dir: Path = BUNDLE_DIR) -> tuple[ValueCheck, ...]:
    """
    Get the value checks for the scenario GHG manuscript

    Parameters
    ----------
    bundle_dir
        Directory which holds the original run's bundle

    Returns
    -------
    :
        The value checks
    """
    summary_file = bundle_dir / MAGICC_SUMMARY_FILE

    return (
        ValueCheck(
            tag="single-projection-magicc-max-erf-difference-sum-l",
            calculate=lambda: get_single_projection_magicc_max_erf_difference_sum(
                summary_file, "l"
            ),
            description=(
                "Sum over the single-projection gases of the largest (over time) "
                "difference in approximate radiative effect "
                "between MAGICC's median and our concentrations, in L"
            ),
        ),
        ValueCheck(
            tag="single-projection-magicc-max-erf-difference-sum-range",
            calculate=lambda: get_single_projection_magicc_max_erf_difference_sum_range(
                summary_file
            ),
            description=(
                "Range across the scenarios of the sum "
                "over the single-projection gases "
                "of the largest (over time) difference in approximate radiative effect "
                "between MAGICC's median and our concentrations"
            ),
        ),
        ValueCheck(
            tag="igcc-last-year",
            calculate=get_igcc_last_year,
            description="Last year of IGCC's record",
        ),
        *(
            ValueCheck(
                tag=f"{gas}-scenariomip-minus-igcc-last-igcc-year",
                calculate=partial(get_scenariomip_minus_igcc_last_igcc_year, gas),
                description=(
                    "Range across the ScenarioMIP scenarios of their global-, "
                    "annual-mean minus IGCC's, in the last year of IGCC's record"
                ),
            )
            for gas in IGCC_GASES
        ),
        ValueCheck(
            tag="ch4-history-minus-igcc-2022",
            calculate=partial(get_history_minus_igcc, "ch4", 2022),
            description=(
                "CMIP7 historical CH4 global-, annual-mean minus IGCC's, in 2022"
            ),
        ),
        ValueCheck(
            tag="ch4-lat-gradient-regression-slope-ratio-to-historical",
            calculate=partial(
                get_ch4_lat_gradient_regression_slope_ratio_to_historical, bundle_dir
            ),
            description=(
                "Slope of our CH4 latitudinal gradient PC0 regression "
                "divided by the historical dataset's"
            ),
        ),
        *(
            ValueCheck(
                tag=f"ch4-lat-gradient-regression-{end}-year",
                calculate=partial(
                    get_ch4_lat_gradient_regression_year,
                    bundle_dir,
                    last=end == "last",
                ),
                description=(
                    f"{end.capitalize()} year our CH4 latitudinal gradient PC0 "
                    "regression is fitted over"
                ),
            )
            for end in ("first", "last")
        ),
        *(
            ValueCheck(
                tag=f"{key}-vl-cf-regression-slope-ratio",
                calculate=partial(get_vl_cf_regression_slope_ratio, key, bundle_dir),
                description=(
                    f"Slope of VL-CF's {key} regression "
                    "divided by the ScenarioMIP scenarios'"
                ),
            )
            for key in VL_CF_REGRESSION_KEYS
        ),
        ValueCheck(
            tag="vl-cf-regression-is-up-to-harmonisation-year",
            calculate=partial(get_vl_cf_regression_reproduced, bundle_dir),
            description=(
                "Fitting up to VL-CF's harmonisation year "
                "reproduces VL-CF's CO2 seasonality change regression"
            ),
        ),
        *(
            ValueCheck(
                tag=f"{CO2_SEASONALITY_REGRESSION}-regression-{end}-year",
                calculate=partial(
                    get_regression_year,
                    CO2_SEASONALITY_REGRESSION,
                    bundle_dir,
                    last=end == "last",
                ),
                description=(
                    f"{end.capitalize()} year the CO2 seasonality change PC "
                    "regression is fitted over (ScenarioMIP scenarios)"
                ),
            )
            for end in ("first", "last")
        ),
        *(
            ValueCheck(
                tag=f"{gas}-{scenario}-pc0-2300",
                calculate=partial(get_pc0_in_year, gas, scenario, 2300, bundle_dir),
                description=f"{gas} first latitudinal gradient PC in 2300, {scenario}",
            )
            for gas in ("n2o", "hfc134a")
            for scenario in ("vl", VL_CF_MARKER[-1])
        ),
        *(
            ValueCheck(
                tag=f"{gas}-{cmip_era.lower()}-{kind}",
                calculate=partial(get_extreme_value, gas, cmip_era, kind),
                description=(
                    f"{kind} {gas} yearly global-mean over {cmip_era}'s experiments "
                    f"and {COMPARISON_YEARS[0]}-{COMPARISON_YEARS[1]}"
                ),
            )
            for gas in ABSTRACT_GASES
            for cmip_era in ("CMIP7", "CMIP6")
            for kind in ("min", "max")
        ),
        *(
            ValueCheck(
                tag=f"{gas}-cmip7-{kind}-year",
                calculate=partial(get_extreme_year, gas, "CMIP7", kind),
                description=f"Year of the {kind} {gas} CMIP7 yearly global-mean",
            )
            for gas in ABSTRACT_GASES
            for kind in ("min", "max")
        ),
        *(
            ValueCheck(
                tag=f"{gas}-cmip7-{kind}-is-{scenario.name}",
                calculate=partial(get_extreme_is_in, gas, "CMIP7", kind, scenario.name),
                description=(
                    f"The {kind} {gas} CMIP7 yearly global-mean is in {scenario.name}"
                ),
            )
            for gas in ABSTRACT_GASES
            for kind in ("min", "max")
            for scenario in CMIP7_SCENARIOS
        ),
        *(
            ValueCheck(
                tag=f"{gas}-cmip7-range-width-relative-to-cmip6",
                calculate=partial(get_range_width_relative_to_cmip6, gas),
                description=f"Width of CMIP7's {gas} range divided by CMIP6's",
            )
            for gas in ABSTRACT_GASES
        ),
        ValueCheck(
            tag="cfc12-cmip7-above-all-cmip6",
            calculate=get_cfc12_cmip7_above_all_cmip6,
            description="CMIP7's CFC-12 is above every CMIP6 scenario's in every year",
        ),
        ValueCheck(
            tag="cfc12-cmip7-relative-to-cmip6-mean-to-2100",
            calculate=partial(get_cfc12_relative_to_cmip6_mean, 2100),
            description=(
                "Mean over 2022-2100 of CMIP7's CFC-12 divided by "
                "the mean of the CMIP6 scenarios', minus one"
            ),
        ),
        ValueCheck(
            tag="n2o-ssp370-minus-h-2500-radiative-effect",
            calculate=partial(
                get_difference_radiative_effect, "n2o", "ssp370", "h", 2500
            ),
            description=(
                "Approximate radiative effect of N2O in SSP3-7.0 minus H, in 2500"
            ),
        ),
        ValueCheck(
            tag="co2-settled-2500-range",
            calculate=partial(
                get_range_in_year, "co2", ("vl-cf", "vl", "l", "ml", "hl"), 2500
            ),
            description=(
                "Range of CO2 in 2500 across the experiments "
                "other than LN (the exception) and M and H (still stabilising)"
            ),
        ),
        ValueCheck(
            tag="ch4-stabilised-2500-range",
            calculate=partial(
                get_range_in_year, "ch4", ("vl-cf", "vl", "ln", "l", "ml", "hl"), 2500
            ),
            description=(
                "Range of CH4 in 2500 across the experiments other than M and H"
            ),
        ),
        ValueCheck(
            tag="ch4-vl-minus-vl-cf-max",
            calculate=partial(get_vl_minus_vl_cf_max, "ch4"),
            description="Largest difference between VL's and VL-CF's CH4",
        ),
        ValueCheck(
            tag="ch4-vl-minus-vl-cf-max-year",
            calculate=partial(get_vl_minus_vl_cf_max_year, "ch4"),
            description="Year of the largest difference between VL's and VL-CF's CH4",
        ),
        ValueCheck(
            tag="single-projection-final-year-linearly-extrapolated",
            calculate=get_single_projection_final_year_extrapolated,
            description=(
                "Every single projection gas' 2500 value is the linear extrapolation "
                "of its 2498 and 2499 values"
            ),
        ),
        ValueCheck(
            tag="western-et-al-2024-2100-max-fraction-of-peak",
            calculate=get_western_et_al_2024_2100_fraction_of_peak,
            description=(
                "Largest fraction of its peak a Western et al. (2024) gas is at in 2100"
            ),
        ),
        ValueCheck(
            tag="western-et-al-2024-2100-radiative-effect",
            calculate=get_western_et_al_2024_2100_radiative_effect,
            description=(
                "Approximate radiative effect of the Western et al. (2024) gases "
                "in 2100"
            ),
        ),
        ValueCheck(
            tag="max-2022-monthly-difference-from-history-radiative-effect",
            calculate=partial(
                get_max_monthly_difference_from_history_radiative_effect, bundle_dir
            ),
            description=(
                "Largest 2022 difference between the ScenarioMIP scenarios "
                "and the historical dataset at native resolution, "
                "over every gas, latitude and month, in approximate radiative effect"
            ),
        ),
        ValueCheck(
            tag="magicc-jump-runs-reproduce-published",
            calculate=partial(get_magicc_jump_reproduces_published, bundle_dir),
            description=(
                "Our re-run of the published MAGICC set-up reproduces "
                "the published output exactly"
            ),
        ),
        ValueCheck(
            tag="magicc-jump-n-runs",
            calculate=partial(get_magicc_jump_n_runs, bundle_dir),
            description="Ensemble members in the MAGICC jump runs",
        ),
        *(
            ValueCheck(
                tag=f"magicc-jump-max-{short_name}-difference-{period_name}",
                calculate=partial(
                    get_magicc_jump_max_difference,
                    short_name,
                    period_name,
                    bundle_dir,
                ),
                description=(
                    f"Largest difference in the {short_name} ensemble median "
                    "between MAGICC with the CMIP7 history and as published, "
                    f"over {period[0]}-{period[1]}"
                ),
            )
            for short_name in JUMP_VARIABLES.values()
            for period_name, period in JUMP_PERIODS.items()
        ),
        ValueCheck(
            tag="native-grid-lat-band-width",
            calculate=partial(get_native_grid_lat_band_width, bundle_dir),
            description="Width of the latitude bands of our native grid",
        ),
        ValueCheck(
            tag="output-first-year",
            calculate=partial(
                get_output_year, SCENARIOMIP_EXPERIMENTS, bundle_dir, last=False
            ),
            description="First year of our published ScenarioMIP data",
        ),
        ValueCheck(
            tag="output-first-year-vl-cf",
            calculate=partial(
                get_output_year, (VL_CF_MARKER[-1],), bundle_dir, last=False
            ),
            description="First year of our published VL-CF data",
        ),
        ValueCheck(
            tag="output-last-year",
            calculate=partial(
                get_output_year, SCENARIOMIP_EXPERIMENTS, bundle_dir, last=True
            ),
            description="Last year of our published ScenarioMIP data (to December)",
        ),
        ValueCheck(
            tag="n-output-variables",
            calculate=partial(get_n_output_variables, bundle_dir),
            description="Number of variables we publish",
        ),
        ValueCheck(
            tag="source-ids",
            calculate=partial(get_source_ids_as_described, bundle_dir),
            description="The published source IDs have the form we describe",
        ),
        ValueCheck(
            tag="output-time-chunks",
            calculate=partial(get_output_time_chunks_as_described, bundle_dir),
            description="The published files are split in time as we describe",
        ),
        ValueCheck(
            tag="magicc-n-ensemble-members",
            calculate=partial(get_magicc_n_ensemble_members, bundle_dir),
            description="Ensemble members in the published MAGICC output",
        ),
        ValueCheck(
            tag="magicc-switch-from-conc-to-emissions-year",
            calculate=partial(
                get_magicc_switch_from_conc_to_emissions_year, bundle_dir
            ),
            description="Year MAGICC switches from concentrations to emissions",
        ),
        *(
            ValueCheck(
                tag=f"{gas}-harmonisation-convergence-time",
                calculate=partial(get_harmonisation_convergence_time, gas, bundle_dir),
                description=f"Convergence time of the harmonisation of {gas}",
            )
            for gas in TWO_EOF_GASES
        ),
        ValueCheck(
            tag="scenariomip-max-relative-difference-2022-monthly-global-mean",
            calculate=partial(
                get_scenariomip_max_relative_difference_monthly_global_mean,
                LAST_HISTORICAL_YEAR,
                bundle_dir,
            ),
            description=(
                "Largest relative difference between the ScenarioMIP experiments' "
                "monthly global-means in 2022"
            ),
        ),
        ValueCheck(
            tag="cfc12-like-gases",
            calculate=partial(get_cfc12_like_gases_as_described, bundle_dir),
            description="The gases processed like CFC-12 are those we list",
        ),
        ValueCheck(
            tag="hfc125-like-gases",
            calculate=partial(get_hfc125_like_gases_as_described, bundle_dir),
            description="The gases processed like HFC-125 are those we list",
        ),
        ValueCheck(
            tag="cfc12-like-wmo-historical-agreement-2022",
            calculate=partial(
                get_single_projection_max_relative_difference_from_history,
                LAST_HISTORICAL_YEAR,
                bundle_dir,
            ),
            description=(
                "Largest relative difference between the single concentration "
                "projections and the historical dataset in 2022"
            ),
        ),
        ValueCheck(
            tag="one-box-lifetimes",
            calculate=partial(get_lifetimes_from_wmo_2022, bundle_dir),
            description="The lifetimes we use are those in WMO (2022) Table A-5",
        ),
        ValueCheck(
            tag="hfc125-largest-erf-one-box-group-2022",
            calculate=partial(
                get_hfc125_largest_radiative_effect_in_group,
                LAST_HISTORICAL_YEAR,
                bundle_dir,
            ),
            description="HFC-125 has the largest radiative effect of the one-box gases",
        ),
        ValueCheck(
            tag="hfc125-erf-2022",
            calculate=partial(get_hfc125_radiative_effect, LAST_HISTORICAL_YEAR),
            description="HFC-125's approximate radiative effect in 2022",
        ),
        ValueCheck(
            tag="second-pc-follows-history-then-constant",
            calculate=partial(get_second_pc_follows_history_then_constant, bundle_dir),
            description=(
                "The second PCs follow history to 2023, "
                "then are constant at their harmonisation year value"
            ),
        ),
        *(
            ValueCheck(
                tag=f"{gas}-lat-gradient-regression-{which}-year",
                calculate=partial(
                    get_regression_year, gas, bundle_dir, last=which == "last"
                ),
                description=f"{which.capitalize()} year of {gas}'s PC regression",
            )
            for gas in ("co2", "n2o")
            for which in ("first", "last")
        ),
        *(
            ValueCheck(
                tag=f"{gas}-lat-gradient-regression-n-constant-pc0-years",
                calculate=partial(get_n_constant_pc0_regression_years, gas, bundle_dir),
                description=(
                    f"Years in {gas}'s PC regression in which the first PC is constant"
                ),
            )
            for gas in TWO_EOF_GASES
        ),
        ValueCheck(
            tag="scenariomip-emissions-last-historical-year",
            calculate=partial(
                get_last_year_emissions_identical, SCENARIOMIP_EXPERIMENTS, bundle_dir
            ),
            description="Last year the ScenarioMIP experiments' emissions are the same",
        ),
        ValueCheck(
            tag="vl-cf-emissions-last-historical-year",
            calculate=partial(
                get_last_year_emissions_identical, ("vl", VL_CF_MARKER[-1]), bundle_dir
            ),
            description="Last year VL-CF's and VL's emissions are the same",
        ),
    )
