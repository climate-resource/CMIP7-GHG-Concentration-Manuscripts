"""
Checks of the values behind the scenario GHG manuscript's statements

See [local.value_checks][] for how these are used.
"""

from __future__ import annotations

from pathlib import Path

import openscm_units

from local.scenario_ghg_forcing_for_cmip7.bundle import (
    BUNDLE_DIR,
    INVERSE_EMISSIONS_FILE,
    load_inverse_emissions,
)
from local.value_checks import CheckValue, ValueCheck


def get_ccl4_inverse_emissions(bundle_dir: Path, year: int) -> CheckValue:
    """
    Get the CCl4 inverse emissions in a given year

    Parameters
    ----------
    bundle_dir
        Directory which holds the original run's bundle

    year
        Year of interest

    Returns
    -------
    :
        The inverse emissions in `year`
    """
    unit = "kt CCl4 / yr"
    inverse_emissions = load_inverse_emissions(
        bundle_dir / INVERSE_EMISSIONS_FILE, unit=unit
    )

    return openscm_units.unit_registry.Quantity(inverse_emissions.loc[year], unit)


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
    return (
        ValueCheck(
            tag="ccl4-inverse-emissions-2100",
            calculate=lambda: get_ccl4_inverse_emissions(bundle_dir, 2100),
            description="CCl4 inverse emissions in 2100",
        ),
    )
