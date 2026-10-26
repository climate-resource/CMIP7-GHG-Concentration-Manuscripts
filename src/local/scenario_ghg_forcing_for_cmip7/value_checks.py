"""
Checks of the values behind the scenario GHG manuscript's statements

See [local.value_checks][] for how these are used.
"""

from __future__ import annotations

from pathlib import Path

import openscm_units

from local.scenario_ghg_forcing_for_cmip7.bundle import BUNDLE_DIR
from local.scenario_ghg_forcing_for_cmip7.methods_data import (
    get_max_radiative_effect_difference,
)
from local.scenario_ghg_forcing_for_cmip7.rerun_projections import (
    MARKERS,
    SINGLE_PROJECTION_GASES,
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
    Get the range across the scenarios of the sum of the largest differences

    See [get_single_projection_magicc_max_erf_difference_sum][].

    Returns a quantity with the lowest and highest values,
    so the check passes only if both are within the statement's bounds.
    """
    values = [
        get_single_projection_magicc_max_erf_difference_sum(summary_file, name).m
        for _, _, name in MARKERS
    ]

    return openscm_units.unit_registry.Quantity([min(values), max(values)], "W / m^2")


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
    )
