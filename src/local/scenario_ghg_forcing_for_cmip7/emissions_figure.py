"""
Generation of the scenario emissions figures

One figure per species, laid out like the yearly panels of the results figures
([local.scenario_ghg_forcing_for_cmip7.results_figure][]),
so the emissions can be read against the concentrations they lead to:

- every CMIP7 scenario's global emissions
- the CMIP7 scenarios' global emissions next to CMIP6's,
  in groups of higher, medium and lower peak scenarios

The CMIP7 emissions are the complete emissions MAGICC was run with
in the original run: the history up to the harmonisation year,
then the scenario and its extension beyond 2100.
The CMIP6 emissions are those collected in RCMIP v5.1.0
(Nicholls et al., 2020, https://doi.org/10.5194/gmd-13-5175-2020),
the emissions behind the CMIP6 concentrations
(Meinshausen et al., 2020, https://doi.org/10.5194/gmd-13-3571-2020):
the history to 2014, then the scenarios and their extensions beyond 2100.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cache
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import openscm_units
import pandas as pd
from loguru import logger
from pandas_openscm.unit_conversion import convert_unit

from local.cmip_ghg_generation import DEFAULT_BUNDLE_DIR as HISTORICAL_BUNDLE_DIR
from local.historical_ghg_forcing_for_cmip7.comparison_data import (
    TIME_COLUMN,
    VALUE_COLUMN,
)
from local.historical_ghg_forcing_for_cmip7.equivalent_species import (
    CMIP7_DEFINITION,
)
from local.historical_ghg_forcing_for_cmip7.layout import (
    create_figure,
    label_panels,
    lay_out_figure,
    manuscript_style,
)
from local.scenario_ghg_forcing_for_cmip7.bundle import BUNDLE_DIR
from local.scenario_ghg_forcing_for_cmip7.magicc_output import MAGICC_GAS_NAMES
from local.scenario_ghg_forcing_for_cmip7.methods_data import (
    MID_YEAR,
    get_complete_emissions,
)
from local.scenario_ghg_forcing_for_cmip7.rerun_projections import (
    HARMONISATION_YEAR,
)
from local.scenario_ghg_forcing_for_cmip7.rerun_vl_cf import VL_CF_MARKER
from local.scenario_ghg_forcing_for_cmip7.results_figure import (
    LAYOUT_SETTINGS,
    LEGEND_PANELS,
    YEARLY_ROWS,
    branch_off_history,
    get_axes_of_panel,
    plot_yearly_panels,
)
from local.scenario_ghg_forcing_for_cmip7.results_figure import (
    TITLES as RESULTS_TITLES,
)
from local.scenario_ghg_forcing_for_cmip7.scenarios import (
    CMIP6_HISTORY_LAST_YEAR,
    CMIP6_SCENARIOS,
    CMIP7_SCENARIOS,
    SCENARIO_GROUPS,
    Scenario,
)

CMIP6_EMISSIONS_FILE = (
    HISTORICAL_BUNDLE_DIR
    / "data"
    / "processed"
    / "historical_emissions"
    / "4e0139cb5adff04ce6ccd33237419e53-rcmip-emissions-annual-means-v5-1-0.csv"
)
"""RCMIP v5.1.0 emissions, which hold the CMIP6 emissions

Part of the historical run's bundle (its `data--processed` tarball).
"""

CMIP6_HISTORY_SCENARIO = "historical"
"""Name of the CMIP6 history in [CMIP6_EMISSIONS_FILE][]"""

GWP_CONTEXT = "AR6GWP100"
"""Context in which emissions of different gases are converted into each other

AR6's 100-year global warming potentials.
"""

UR = openscm_units.unit_registry
"""Unit registry to use for the conversions"""


@dataclass(frozen=True)
class EmissionsSpecies:
    """
    A species whose emissions we show
    """

    components: tuple[str, ...]
    """
    Gases whose emissions make up the species' emissions

    Converted into [units][(c).] with [GWP_CONTEXT][], then summed.
    """

    units: str
    """Units to show the emissions in"""

    units_label: str
    """Units, as written on the vertical axis"""


EMISSIONS_SPECIES = {
    "co2": EmissionsSpecies(("co2",), "Gt CO2 / yr", "Gt CO$_2$ / yr"),
    "ch4": EmissionsSpecies(("ch4",), "Mt CH4 / yr", "Mt CH$_4$ / yr"),
    "n2o": EmissionsSpecies(("n2o",), "Mt N2O / yr", "Mt N$_2$O / yr"),
    # The same gases as the HFC-134a-eq concentrations.
    # The concentrations are combined with radiative efficiencies,
    # the emissions with GWP100s.
    "hfc134aeq": EmissionsSpecies(
        CMIP7_DEFINITION.components["hfc134aeq"],
        "kt HFC134a / yr",
        "kt HFC-134a-eq / yr",
    ),
}
"""Species whose emissions we show, by name"""

TITLES = {
    "cmip7": "CMIP7 global emissions",
    **{group: RESULTS_TITLES[group] for group in SCENARIO_GROUPS},
}
"""Title of each panel"""


def get_cmip7_variables(gas: str) -> tuple[str, ...]:
    """
    Get the variables which hold a gas' emissions in the CMIP7 emissions

    Parameters
    ----------
    gas
        Gas of interest

    Returns
    -------
    :
        Variables whose sum is the gas' emissions
    """
    if gas == "co2":
        return ("Emissions|CO2|Fossil", "Emissions|CO2|Biosphere")

    return (f"Emissions|{MAGICC_GAS_NAMES[gas]}",)


def get_cmip6_variable(gas: str, variables: pd.Index) -> str:
    """
    Get the variable which holds a gas' emissions in the CMIP6 emissions

    Parameters
    ----------
    gas
        Gas of interest

    variables
        Variables in the CMIP6 emissions

    Returns
    -------
    :
        The variable, e.g. "Emissions|CO2" or "Emissions|F-Gases|HFC|HFC134a"
    """
    matches = [v for v in variables if v.split("|")[-1] == MAGICC_GAS_NAMES[gas]]
    if len(matches) != 1:
        msg = f"Expected exactly one CMIP6 variable for {gas}, found {matches}"
        raise AssertionError(msg)

    return matches[0]


def sum_in_units(
    emissions: pd.DataFrame, variables: list[str], units: str
) -> pd.Series:
    """
    Sum emissions of several variables, in common units

    Parameters
    ----------
    emissions
        Emissions, with `variable` and `unit` index levels and years as columns

    variables
        Variables to sum

    units
        Units to convert to before summing, with [GWP_CONTEXT][]

    Returns
    -------
    :
        Summed emissions, indexed by year

        Years in which any variable has no value have no value.
    """
    selected = emissions[emissions.index.get_level_values("variable").isin(variables)]
    if sorted(selected.index.get_level_values("variable")) != sorted(variables):
        msg = f"Expected one row for each of {variables}, found {selected.index}"
        raise AssertionError(msg)

    with UR.context(GWP_CONTEXT):
        converted = convert_unit(selected, units, ur=UR)

    res = converted.sum(axis="index", min_count=len(variables))
    res.index = res.index.astype(int)

    return res


def get_cmip7_emissions(
    species: EmissionsSpecies, scenario: str, bundle_dir: Path = BUNDLE_DIR
) -> pd.Series:
    """
    Get a species' emissions in a CMIP7 scenario, including the history

    Parameters
    ----------
    species
        Species of interest

    scenario
        Scenario of interest

    bundle_dir
        Directory which holds the bundle

    Returns
    -------
    :
        Emissions, indexed by year
    """
    variables = [v for gas in species.components for v in get_cmip7_variables(gas)]

    return sum_in_units(
        get_complete_emissions(scenario, bundle_dir), variables, species.units
    )


def get_cmip7_history(
    species: EmissionsSpecies, bundle_dir: Path = BUNDLE_DIR
) -> pd.Series:
    """
    Get the history behind a species' CMIP7 emissions

    Parameters
    ----------
    species
        Species of interest

    bundle_dir
        Directory which holds the bundle

    Returns
    -------
    :
        Emissions up to the harmonisation year, indexed by year
    """
    histories = [
        get_cmip7_emissions(species, s.name, bundle_dir).loc[:HARMONISATION_YEAR]
        # VL-CF's fossil CO2 emissions leave the history in 2016
        for s in CMIP7_SCENARIOS
        if s.name != VL_CF_MARKER[2]
    ]
    for history in histories[1:]:
        if not np.allclose(history, histories[0]):
            msg = "The CMIP7 scenarios don't share one history"
            raise AssertionError(msg)

    return histories[0]


def get_cmip7_scenario(
    species: EmissionsSpecies,
    scenario: Scenario,
    history: pd.Series,
    bundle_dir: Path = BUNDLE_DIR,
) -> pd.Series:
    """
    Get a species' emissions in a CMIP7 scenario, from where it leaves the history

    Parameters
    ----------
    species
        Species of interest

    scenario
        Scenario of interest

    history
        History, as returned by [get_cmip7_history][]

    bundle_dir
        Directory which holds the bundle

    Returns
    -------
    :
        Emissions from the first year in which they differ from the history
        (the year after the harmonisation year for most scenarios,
        earlier for VL-CF's CO2), indexed by year
    """
    emissions = get_cmip7_emissions(species, scenario.name, bundle_dir)
    overlap = emissions.loc[history.index]
    differs = ~np.isclose(overlap, history)
    first_year = (
        overlap.index[differs].min() if differs.any() else history.index.max() + 1
    )

    return emissions.loc[first_year:]


@cache
def load_cmip6_emissions(file: Path) -> pd.DataFrame:
    """
    Load the global CMIP6 emissions (cached, they are used repeatedly)

    Parameters
    ----------
    file
        File in which the emissions are saved (RCMIP format)

    Returns
    -------
    :
        Emissions, with `scenario`, `variable` and `unit` index levels
        and years as columns
    """
    raw = pd.read_csv(file)
    res = (
        raw[raw["Region"] == "World"]
        .drop(columns=["Model", "Region", "Mip_Era", "Activity_Id"])
        .rename(
            columns={"Scenario": "scenario", "Variable": "variable", "Unit": "unit"}
        )
        .set_index(["scenario", "variable", "unit"])
    )
    res.columns = res.columns.astype(int)

    return res


def get_cmip6_emissions(
    species: EmissionsSpecies, scenario: str, file: Path = CMIP6_EMISSIONS_FILE
) -> pd.Series:
    """
    Get a species' emissions in a CMIP6 scenario (or the CMIP6 history)

    Parameters
    ----------
    species
        Species of interest

    scenario
        Scenario of interest, or [CMIP6_HISTORY_SCENARIO][]

    file
        File in which the CMIP6 emissions are saved

    Returns
    -------
    :
        Emissions, indexed by year

        After 2100, the emissions are only given every ten years,
        so we interpolate linearly between them to get yearly values.
    """
    emissions = load_cmip6_emissions(file).loc[scenario]
    variables = [
        get_cmip6_variable(gas, emissions.index.get_level_values("variable"))
        for gas in species.components
    ]

    return (
        sum_in_units(emissions, variables, species.units)
        .interpolate(method="index", limit_area="inside")
        .dropna()
    )


def to_frame(emissions: pd.Series) -> pd.DataFrame:
    """
    Get emissions as a frame, placed at the middle of each year

    The same way the yearly concentrations are placed,
    see [local.scenario_ghg_forcing_for_cmip7.esgf_data.to_decimal_year][].

    Parameters
    ----------
    emissions
        Emissions, indexed by year

    Returns
    -------
    :
        Emissions with a [TIME_COLUMN][] and a [VALUE_COLUMN][]
    """
    return pd.DataFrame(
        {
            TIME_COLUMN: emissions.index.astype(float) + MID_YEAR,
            VALUE_COLUMN: emissions.to_numpy(),
        }
    )


@manuscript_style
def generate_emissions_figure(
    species_name: str,
    outfile: Path,
    bundle_dir: Path = BUNDLE_DIR,
    cmip6_emissions_file: Path = CMIP6_EMISSIONS_FILE,
) -> Path:
    """
    Generate a scenario emissions figure

    Parameters
    ----------
    species_name
        Species of interest, a key of [EMISSIONS_SPECIES][]

    outfile
        File in which to write the figure

    bundle_dir
        Directory which holds the bundle

    cmip6_emissions_file
        File in which the CMIP6 emissions are saved

    Returns
    -------
    :
        `outfile`
    """
    species = EMISSIONS_SPECIES[species_name]
    fig, axes = create_figure(YEARLY_ROWS, settings=LAYOUT_SETTINGS)

    cmip7_history = get_cmip7_history(species, bundle_dir)
    cmip6_history = get_cmip6_emissions(
        species, CMIP6_HISTORY_SCENARIO, cmip6_emissions_file
    ).loc[:CMIP6_HISTORY_LAST_YEAR]
    histories = {"CMIP7": to_frame(cmip7_history), "CMIP6": to_frame(cmip6_history)}

    scenarios = {
        **{
            s: branch_off_history(
                s,
                to_frame(get_cmip7_scenario(species, s, cmip7_history, bundle_dir)),
                histories["CMIP7"],
            )
            for s in CMIP7_SCENARIOS
        },
        **{
            s: branch_off_history(
                s,
                to_frame(
                    get_cmip6_emissions(species, s.name, cmip6_emissions_file).loc[
                        CMIP6_HISTORY_LAST_YEAR + 1 :
                    ]
                ),
                histories["CMIP6"],
            )
            for s in CMIP6_SCENARIOS
        },
    }
    plot_yearly_panels(YEARLY_ROWS, axes, histories, scenarios, species.units_label)

    # Some species' emissions go negative, so show where zero is
    for panel in ("cmip7", *SCENARIO_GROUPS):
        for ax in get_axes_of_panel(YEARLY_ROWS, axes, panel):
            if ax.get_ylim()[0] < 0.0:
                ax.axhline(0.0, color="0.7", linewidth=0.5, zorder=0)

    label_panels(YEARLY_ROWS, axes, TITLES, unlabelled=tuple(LEGEND_PANELS.values()))

    # Last, because it needs to know how much room everything takes up
    lay_out_figure(fig, axes, YEARLY_ROWS, settings=LAYOUT_SETTINGS)

    outfile.parent.mkdir(exist_ok=True, parents=True)
    logger.info(f"Writing {outfile}")
    fig.savefig(outfile)
    plt.close(fig)

    return outfile
