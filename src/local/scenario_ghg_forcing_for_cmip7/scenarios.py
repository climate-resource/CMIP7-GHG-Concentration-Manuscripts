"""
The scenarios we show, and how we show them

The single source of truth for each scenario's name, colour, group
and the ESGF datasets which hold it,
so that every figure draws a scenario the same way.
"""

from __future__ import annotations

from dataclasses import dataclass

CMIP7_HISTORY_SOURCE_ID = "CR-CMIP-1-0-0"
"""Source ID of the CMIP7 historical concentrations"""

CMIP6_HISTORY_SOURCE_ID = "UoM-CMIP-1-2-0"
"""Source ID of the CMIP6 historical concentrations"""

CMIP7_HISTORY_LAST_YEAR = 2021
"""Last year of the CMIP7 history we show

The CMIP7 historical files run to the end of 2022,
but the scenario files start at the beginning of 2022
(and differ slightly from the history in that year).
The scenarios are what the models run from 2022, so we show those.
VL-CF starts in 2015, so it is drawn over the end of the history we show.
"""

CMIP6_HISTORY_LAST_YEAR = 2014
"""Last year of the CMIP6 history (the scenarios start in 2015)"""

HISTORY_COLOUR = "#000000"
"""Colour of the history (CMIP6's and CMIP7's), told apart by line style"""

CMIP7_LINESTYLE = "-"
"""Line style for CMIP7 data, history and scenarios alike"""

CMIP6_LINESTYLE = "--"
"""Line style for CMIP6 data, history and scenarios alike"""


@dataclass(frozen=True)
class Scenario:
    """
    A scenario we show
    """

    name: str
    """Name of the scenario in the data, e.g. "vl" or "ssp119\""""

    label: str
    """Name of the scenario as shown to a reader"""

    cmip_era: str
    """CMIP era the scenario belongs to"""

    colour: str
    """Colour to draw the scenario in"""

    group: str
    """Group the scenario is shown in, one of [SCENARIO_GROUPS][]"""

    source_ids: tuple[str, ...]
    """
    Source IDs of the ESGF datasets which hold the scenario, in time order

    CMIP7 splits each scenario into the part to 2100 and its extension,
    CMIP6 has one dataset per scenario.
    """

    @property
    def linestyle(self) -> str:
        """Line style to draw the scenario in"""
        return CMIP7_LINESTYLE if self.cmip_era == "CMIP7" else CMIP6_LINESTYLE


SCENARIO_GROUPS = ("higher-peak", "medium-peak", "lower-peak")
"""Groups of scenarios, in the order they are shown (top to bottom)"""

# TODO: get official colour palettes.
# These are the colours we have been using so far
# (`notebooks/cci-2026-poster.py`, also the GHG viewer).
# The CMIP6 ones are the IPCC AR6 SSP colours.
CMIP7_SCENARIOS = (
    # PolMIP's counterfactual to VL (emissions diverge from history from 2016).
    # TODO: get an official colour for VL-CF.
    # Grey, as a counterfactual, and to stay clear of the other low scenarios' colours.
    Scenario(
        "vl-cf",
        "VL-CF",
        "CMIP7",
        "#8c8c8c",
        "lower-peak",
        ("CR-vl-cf-1-1-0", "CR-vl-cf-ext-1-1-0"),
    ),
    Scenario(
        "vl", "VL", "CMIP7", "#24a4ff", "lower-peak", ("CR-vl-1-1-0", "CR-vl-ext-1-1-0")
    ),
    Scenario(
        "ln", "LN", "CMIP7", "#4a0daf", "lower-peak", ("CR-ln-1-1-0", "CR-ln-ext-1-1-0")
    ),
    Scenario(
        "l", "L", "CMIP7", "#00cc69", "lower-peak", ("CR-l-1-1-0", "CR-l-ext-1-1-0")
    ),
    Scenario(
        "ml",
        "ML",
        "CMIP7",
        "#f5ac00",
        "medium-peak",
        ("CR-ml-1-1-0", "CR-ml-ext-1-1-0"),
    ),
    Scenario(
        "m", "M", "CMIP7", "#ffa9dc", "medium-peak", ("CR-m-1-1-0", "CR-m-ext-1-1-0")
    ),
    Scenario(
        "hl",
        "HL",
        "CMIP7",
        "#8f003b",
        "higher-peak",
        ("CR-hl-1-1-0", "CR-hl-ext-1-1-0"),
    ),
    Scenario(
        "h", "H", "CMIP7", "#700000", "higher-peak", ("CR-h-1-1-0", "CR-h-ext-1-1-0")
    ),
)
"""CMIP7 scenarios, lowest to highest"""

CMIP6_SCENARIOS = (
    Scenario(
        "ssp119",
        "SSP1-1.9",
        "CMIP6",
        "#00a9cf",
        "lower-peak",
        ("UoM-IMAGE-ssp119-1-2-1",),
    ),
    Scenario(
        "ssp126",
        "SSP1-2.6",
        "CMIP6",
        "#003466",
        "lower-peak",
        ("UoM-IMAGE-ssp126-1-2-1",),
    ),
    Scenario(
        "ssp434",
        "SSP4-3.4",
        "CMIP6",
        "#2274ae",
        "medium-peak",
        ("UoM-GCAM4-ssp434-1-2-1",),
    ),
    Scenario(
        "ssp245",
        "SSP2-4.5",
        "CMIP6",
        "#f69320",
        "medium-peak",
        ("UoM-MESSAGE-GLOBIOM-ssp245-1-2-1",),
    ),
    Scenario(
        "ssp460",
        "SSP4-6.0",
        "CMIP6",
        "#b0724e",
        "medium-peak",
        ("UoM-GCAM4-ssp460-1-2-1",),
    ),
    Scenario(
        "ssp534-over",
        "SSP5-3.4-over",
        "CMIP6",
        "#92397a",
        "higher-peak",
        ("UoM-REMIND-MAGPIE-ssp534-over-1-2-1",),
    ),
    Scenario(
        "ssp370",
        "SSP3-7.0",
        "CMIP6",
        "#df0000",
        "higher-peak",
        ("UoM-AIM-ssp370-1-2-1",),
    ),
    Scenario(
        "ssp585",
        "SSP5-8.5",
        "CMIP6",
        "#980002",
        "higher-peak",
        ("UoM-REMIND-MAGPIE-ssp585-1-2-1",),
    ),
)
"""CMIP6 scenarios, lowest to highest (roughly)"""

SCENARIOS = {
    scenario.name: scenario for scenario in (*CMIP7_SCENARIOS, *CMIP6_SCENARIOS)
}
"""Every scenario, by name"""

if len(SCENARIOS) != len(CMIP7_SCENARIOS) + len(CMIP6_SCENARIOS):
    msg = "Scenario names are not unique"
    raise AssertionError(msg)

if not {s.group for s in SCENARIOS.values()} <= set(SCENARIO_GROUPS):
    msg = "Unknown scenario group"
    raise AssertionError(msg)


def get_scenarios_in_group(
    group: str, scenarios: tuple[Scenario, ...]
) -> tuple[Scenario, ...]:
    """
    Get the scenarios in a group

    Parameters
    ----------
    group
        Group of interest, one of [SCENARIO_GROUPS][]

    scenarios
        Scenarios to choose from

    Returns
    -------
    :
        The scenarios in `scenarios` which are in `group`, in the same order
    """
    if group not in SCENARIO_GROUPS:
        raise ValueError(group)

    return tuple(s for s in scenarios if s.group == group)
