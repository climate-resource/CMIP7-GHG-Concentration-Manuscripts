"""
Checks of the values behind the scenario GHG manuscript's statements

See [local.value_checks][] for how these are used.
"""

from __future__ import annotations

from functools import partial
from pathlib import Path

import numpy as np
import openscm_units
import pandas as pd
import xarray as xr

from local.historical_ghg_forcing_for_cmip7.comparison_data import (
    TIME_COLUMN,
    VALUE_COLUMN,
    get_igcc_comparison,
)
from local.scenario_ghg_forcing_for_cmip7.bundle import BUNDLE_DIR
from local.scenario_ghg_forcing_for_cmip7.esgf_data import (
    get_history_yearly_global_mean,
    get_scenario_yearly_global_mean,
)
from local.scenario_ghg_forcing_for_cmip7.methods_data import (
    get_max_radiative_effect_difference,
)
from local.scenario_ghg_forcing_for_cmip7.rerun_projections import (
    MARKERS,
    SINGLE_PROJECTION_GASES,
)
from local.scenario_ghg_forcing_for_cmip7.rerun_vl_cf import VL_CF_MARKER
from local.scenario_ghg_forcing_for_cmip7.scenarios import CMIP7_SCENARIOS
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
    )
