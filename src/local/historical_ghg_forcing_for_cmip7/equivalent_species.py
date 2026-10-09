r"""
Comparison of the equivalent species in different datasets

An equivalent species is a radiative efficiency weighted sum
of the concentrations of a group of gases, expressed as a concentration
of a reference gas (e.g. CFC-12 for CFC-12 equivalent, `cfc12eq`):

$$
C_{\text{eq}} = \sum_{g \in G} C_g \frac{R_g}{R_{\text{ref}}}
$$

where $G$ is the group of gases included,
$C_g$ is the concentration of gas $g$
and $R_g$ is its radiative efficiency.

Two datasets' equivalent species can therefore differ for three reasons:

1. they include different groups of gases ($G$)
1. they use different radiative efficiencies ($R$)
1. their concentrations of the individual gases differ ($C$)

[decompose_difference][] splits the difference into these three parts.
Each dataset's definition ([EquivalenceDefinition][])
is copied from the code which produced it,
and [check_reproduction][] checks that the definition
reproduces the dataset's published equivalent species,
so we know we have the definition right.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from pathlib import Path

import matplotlib.figure
import matplotlib.pyplot as plt
import openscm_units
import pandas as pd
import pint
import xarray as xr

from local.cmip_ghg_generation import DEFAULT_BUNDLE_DIR
from local.historical_ghg_forcing_for_cmip7.comparison_data import (
    IGCC_COLUMNS,
    IGCC_CONCENTRATIONS_FILE,
    IGCC_CONCENTRATIONS_URL,
    RADIATIVE_EFFICIENCIES,
    ensure_file_downloaded,
    get_cmip6_spatial_means,
    get_igcc_units,
)
from local.historical_ghg_forcing_for_cmip7.plotting import (
    OKABE_ITO,
    get_only_data_variable,
    label_name,
)

Q = openscm_units.unit_registry.Quantity

EQUIVALENT_SPECIES = ("cfc11eq", "cfc12eq", "hfc134aeq")
"""The equivalent species we produce"""

UNITS = "ppt"
"""Units every concentration is converted to before it is used here"""


@dataclass(frozen=True)
class EquivalenceDefinition:
    """
    How a dataset defines its equivalent species
    """

    label: str
    """Name of the dataset"""

    components: Mapping[str, tuple[str, ...]]
    """Gases included in each equivalent species, keyed by the equivalent species"""

    radiative_efficiencies: Mapping[str, pint.Quantity]
    """Radiative efficiency of each gas"""

    def get_reference_gas(self, equivalent_species: str) -> str:
        """
        Get the gas an equivalent species is expressed as a concentration of

        Parameters
        ----------
        equivalent_species
            Equivalent species of interest e.g. `cfc12eq`

        Returns
        -------
        :
            Reference gas e.g. `cfc12`
        """
        return equivalent_species.removesuffix("eq")


# Gas names follow ours throughout (e.g. `cfc12`, `hfc4310mee`),
# whatever the dataset itself calls them.

CMIP7_DEFINITION = EquivalenceDefinition(
    label="CMIP7",
    # Mirrors `EQUIVALENT_COMPONENTS` in the original run's
    # `src/local/config_creation/crunch_equivalent_species.py`
    components={
        "cfc11eq": (
            "c2f6",
            "c3f8",
            "c4f10",
            "c5f12",
            "c6f14",
            "c7f16",
            "c8f18",
            "cc4f8",
            "ccl4",
            "cf4",
            "cfc11",
            "cfc113",
            "cfc114",
            "cfc115",
            "ch2cl2",
            "ch3br",
            "ch3ccl3",
            "ch3cl",
            "chcl3",
            "halon1211",
            "halon1301",
            "halon2402",
            "hcfc141b",
            "hcfc142b",
            "hcfc22",
            "hfc125",
            "hfc134a",
            "hfc143a",
            "hfc152a",
            "hfc227ea",
            "hfc23",
            "hfc236fa",
            "hfc245fa",
            "hfc32",
            "hfc365mfc",
            "hfc4310mee",
            "nf3",
            "sf6",
            "so2f2",
        ),
        "cfc12eq": (
            "cfc11",
            "cfc113",
            "cfc114",
            "cfc115",
            "cfc12",
            "ccl4",
            "ch2cl2",
            "ch3br",
            "ch3ccl3",
            "ch3cl",
            "chcl3",
            "halon1211",
            "halon1301",
            "halon2402",
            "hcfc141b",
            "hcfc142b",
            "hcfc22",
        ),
        "hfc134aeq": (
            "c2f6",
            "c3f8",
            "c4f10",
            "c5f12",
            "c6f14",
            "c7f16",
            "c8f18",
            "cc4f8",
            "cf4",
            "hfc125",
            "hfc134a",
            "hfc143a",
            "hfc152a",
            "hfc227ea",
            "hfc23",
            "hfc236fa",
            "hfc245fa",
            "hfc32",
            "hfc365mfc",
            "hfc4310mee",
            "nf3",
            "sf6",
            "so2f2",
        ),
    },
    # AR6, see RADIATIVE_EFFICIENCIES
    radiative_efficiencies=RADIATIVE_EFFICIENCIES,
)
"""Our definition of the equivalent species"""

CMIP6_DEFINITION = EquivalenceDefinition(
    label="CMIP6",
    # From `create_equivalent_ODS_smallerGHG_conc.m`
    # in the code which produced the CMIP6 concentrations
    # (Meinshausen et al., 2017, https://doi.org/10.5194/gmd-10-2057-2017).
    # The CFC-12 and HFC-134a equivalents are its `allthree` case,
    # the CFC-11 equivalent is its `keythree` case.
    # They are the same groups of gases as ours.
    components={
        "cfc11eq": (
            "hfc134a",
            "hfc23",
            "hfc32",
            "hfc125",
            "hfc143a",
            "hfc152a",
            "hfc227ea",
            "hfc236fa",
            "hfc245fa",
            "hfc365mfc",
            "hfc4310mee",
            "nf3",
            "sf6",
            "so2f2",
            "cf4",
            "c2f6",
            "c3f8",
            "c4f10",
            "c5f12",
            "c6f14",
            "c7f16",
            "c8f18",
            "cc4f8",
            "cfc11",
            "cfc113",
            "cfc114",
            "cfc115",
            "hcfc22",
            "hcfc141b",
            "hcfc142b",
            "ch3ccl3",
            "ccl4",
            "ch3cl",
            "ch2cl2",
            "chcl3",
            "ch3br",
            "halon1211",
            "halon1301",
            "halon2402",
        ),
        "cfc12eq": (
            "cfc12",
            "cfc11",
            "cfc113",
            "cfc114",
            "cfc115",
            "hcfc22",
            "hcfc141b",
            "hcfc142b",
            "ch3ccl3",
            "ccl4",
            "ch3cl",
            "ch2cl2",
            "chcl3",
            "ch3br",
            "halon1211",
            "halon1301",
            "halon2402",
        ),
        "hfc134aeq": (
            "hfc134a",
            "hfc23",
            "hfc32",
            "hfc125",
            "hfc143a",
            "hfc152a",
            "hfc227ea",
            "hfc236fa",
            "hfc245fa",
            "hfc365mfc",
            "hfc4310mee",
            "nf3",
            "sf6",
            "so2f2",
            "cf4",
            "c2f6",
            "c3f8",
            "c4f10",
            "c5f12",
            "c6f14",
            "c7f16",
            "c8f18",
            "cc4f8",
        ),
    },
    # From the `Substances_overview` sheet of `CMIP6_Substances_data_overview.xlsx`,
    # the table `create_equivalent_ODS_smallerGHG_conc.m` reads.
    # That sheet is Appendix 8.A of the IPCC's Fifth Assessment Report
    # (Myhre et al., 2013).
    radiative_efficiencies={
        "cfc11": Q(0.26, "W / m^2 / ppb"),
        "cfc12": Q(0.32, "W / m^2 / ppb"),
        "cfc113": Q(0.3, "W / m^2 / ppb"),
        "cfc114": Q(0.31, "W / m^2 / ppb"),
        "cfc115": Q(0.2, "W / m^2 / ppb"),
        "hcfc22": Q(0.21, "W / m^2 / ppb"),
        "hcfc141b": Q(0.16, "W / m^2 / ppb"),
        "hcfc142b": Q(0.19, "W / m^2 / ppb"),
        "ch3ccl3": Q(0.07, "W / m^2 / ppb"),
        "ccl4": Q(0.17, "W / m^2 / ppb"),
        "ch3cl": Q(0.01, "W / m^2 / ppb"),
        "ch2cl2": Q(0.03, "W / m^2 / ppb"),
        "chcl3": Q(0.08, "W / m^2 / ppb"),
        "ch3br": Q(0.004, "W / m^2 / ppb"),
        "halon1211": Q(0.29, "W / m^2 / ppb"),
        "halon1301": Q(0.3, "W / m^2 / ppb"),
        "halon2402": Q(0.31, "W / m^2 / ppb"),
        "hfc134a": Q(0.16, "W / m^2 / ppb"),
        "hfc23": Q(0.18, "W / m^2 / ppb"),
        "hfc32": Q(0.11, "W / m^2 / ppb"),
        "hfc125": Q(0.23, "W / m^2 / ppb"),
        "hfc143a": Q(0.16, "W / m^2 / ppb"),
        "hfc152a": Q(0.1, "W / m^2 / ppb"),
        "hfc227ea": Q(0.26, "W / m^2 / ppb"),
        "hfc236fa": Q(0.24, "W / m^2 / ppb"),
        "hfc245fa": Q(0.24, "W / m^2 / ppb"),
        "hfc365mfc": Q(0.22, "W / m^2 / ppb"),
        "hfc4310mee": Q(0.42, "W / m^2 / ppb"),
        "nf3": Q(0.2, "W / m^2 / ppb"),
        "sf6": Q(0.57, "W / m^2 / ppb"),
        "so2f2": Q(0.2, "W / m^2 / ppb"),
        "cf4": Q(0.09, "W / m^2 / ppb"),
        "c2f6": Q(0.25, "W / m^2 / ppb"),
        "c3f8": Q(0.28, "W / m^2 / ppb"),
        "c4f10": Q(0.36, "W / m^2 / ppb"),
        "c5f12": Q(0.41, "W / m^2 / ppb"),
        "c6f14": Q(0.44, "W / m^2 / ppb"),
        "c7f16": Q(0.5, "W / m^2 / ppb"),
        "c8f18": Q(0.55, "W / m^2 / ppb"),
        "cc4f8": Q(0.32, "W / m^2 / ppb"),
    },
)
"""The definition of the equivalent species used to produce the CMIP6 dataset"""

IGCC_EXTRA_COLUMNS = {
    "cfc13": "CFC-13",
    "cfc112": "CFC-112",
    "cfc112a": "CFC-112a",
    "cfc113a": "CFC-113a",
    "cfc114a": "CFC-114a",
    "hcfc133a": "HCFC-133a",
    "hcfc31": "HCFC-31",
    "hcfc124": "HCFC-124",
}
"""Column of IGCC's concentrations file of each gas IGCC includes but we don't"""

IGCC_DEFINITION = EquivalenceDefinition(
    label="IGCC",
    # Mirrors `gases_montreal` and `gases_hfcs`
    # in `notebooks/01_trace-gas-global-mean.py`
    # of https://github.com/ClimateIndicator/forcing-timeseries at v6.4.0.
    # IGCC has no CFC-11 equivalent.
    components={
        "cfc12eq": (
            "cfc12",
            "cfc11",
            "cfc113",
            "cfc114",
            "cfc115",
            "cfc13",
            "hcfc22",
            "hcfc141b",
            "hcfc142b",
            "ch3ccl3",
            "ccl4",
            "ch3cl",
            "ch3br",
            "ch2cl2",
            "chcl3",
            "halon1211",
            "halon1301",
            "halon2402",
            "cfc112",
            "cfc112a",
            "cfc113a",
            "cfc114a",
            "hcfc133a",
            "hcfc31",
            "hcfc124",
        ),
        "hfc134aeq": (
            "hfc134a",
            "hfc23",
            "hfc32",
            "hfc125",
            "hfc143a",
            "hfc152a",
            "hfc227ea",
            "hfc236fa",
            "hfc245fa",
            "hfc365mfc",
            "hfc4310mee",
        ),
    },
    # Mirrors `radeff` in the same notebook,
    # which cites Hodnebrog et al. (2020), https://doi.org/10.1029/2019RG000691.
    # Only the gases in our or IGCC's equivalent species are copied.
    # IGCC's C6F14 is n-C6F14, the isomer we match ours to.
    radiative_efficiencies={
        "hfc125": Q(0.23378, "W / m^2 / ppb"),
        "hfc134a": Q(0.16714, "W / m^2 / ppb"),
        "hfc143a": Q(0.168, "W / m^2 / ppb"),
        "hfc152a": Q(0.10174, "W / m^2 / ppb"),
        "hfc227ea": Q(0.27325, "W / m^2 / ppb"),
        "hfc23": Q(0.19111, "W / m^2 / ppb"),
        "hfc236fa": Q(0.25069, "W / m^2 / ppb"),
        "hfc245fa": Q(0.24498, "W / m^2 / ppb"),
        "hfc32": Q(0.11144, "W / m^2 / ppb"),
        "hfc365mfc": Q(0.22813, "W / m^2 / ppb"),
        "hfc4310mee": Q(0.35731, "W / m^2 / ppb"),
        "nf3": Q(0.20448, "W / m^2 / ppb"),
        "c2f6": Q(0.26105, "W / m^2 / ppb"),
        "c3f8": Q(0.26999, "W / m^2 / ppb"),
        "c4f10": Q(0.36874, "W / m^2 / ppb"),
        "c5f12": Q(0.4076, "W / m^2 / ppb"),
        "c6f14": Q(0.44888, "W / m^2 / ppb"),
        "c7f16": Q(0.50312, "W / m^2 / ppb"),
        "c8f18": Q(0.55787, "W / m^2 / ppb"),
        "cf4": Q(0.09859, "W / m^2 / ppb"),
        "cc4f8": Q(0.31392, "W / m^2 / ppb"),
        "sf6": Q(0.56657, "W / m^2 / ppb"),
        "so2f2": Q(0.21074, "W / m^2 / ppb"),
        "ccl4": Q(0.16616, "W / m^2 / ppb"),
        "cfc11": Q(0.25941, "W / m^2 / ppb"),
        "cfc112": Q(0.28192, "W / m^2 / ppb"),
        "cfc112a": Q(0.24564, "W / m^2 / ppb"),
        "cfc113": Q(0.30142, "W / m^2 / ppb"),
        "cfc113a": Q(0.24094, "W / m^2 / ppb"),
        "cfc114": Q(0.31433, "W / m^2 / ppb"),
        "cfc114a": Q(0.29747, "W / m^2 / ppb"),
        "cfc115": Q(0.24625, "W / m^2 / ppb"),
        "cfc12": Q(0.31998, "W / m^2 / ppb"),
        "cfc13": Q(0.27752, "W / m^2 / ppb"),
        "ch2cl2": Q(0.02882, "W / m^2 / ppb"),
        "ch3br": Q(0.00432, "W / m^2 / ppb"),
        "ch3ccl3": Q(0.06454, "W / m^2 / ppb"),
        "ch3cl": Q(0.00466, "W / m^2 / ppb"),
        "chcl3": Q(0.07357, "W / m^2 / ppb"),
        "hcfc124": Q(0.20721, "W / m^2 / ppb"),
        "hcfc133a": Q(0.14995, "W / m^2 / ppb"),
        "hcfc141b": Q(0.16065, "W / m^2 / ppb"),
        "hcfc142b": Q(0.19329, "W / m^2 / ppb"),
        "hcfc22": Q(0.21385, "W / m^2 / ppb"),
        "hcfc31": Q(0.068, "W / m^2 / ppb"),
        "halon1211": Q(0.30014, "W / m^2 / ppb"),
        "halon1301": Q(0.29943, "W / m^2 / ppb"),
        "halon2402": Q(0.31169, "W / m^2 / ppb"),
    },
)
"""The definition of the equivalent species used to produce the IGCC dataset"""


@dataclass(frozen=True)
class EquivalenceDataset:
    """
    A dataset's equivalent species, and everything needed to explain them
    """

    definition: EquivalenceDefinition
    """How the dataset defines its equivalent species"""

    concentrations: pd.DataFrame
    """
    The dataset's global-, annual-mean concentrations of the individual gases

    One column per gas, indexed by (integer) year, in [UNITS][].
    """

    published: pd.DataFrame
    """
    The dataset's published global-, annual-mean equivalent species

    One column per equivalent species, indexed by (integer) year, in [UNITS][].
    """

    @property
    def label(self) -> str:
        """Name of the dataset"""
        return self.definition.label


def to_units(values: Iterable[float], units: str) -> list[float]:
    """
    Convert values to [UNITS][]

    Parameters
    ----------
    values
        Values to convert

    units
        Units `values` are in

    Returns
    -------
    :
        `values`, in [UNITS][]
    """
    return list(Q(list(values), units).to(UNITS).m)


def get_weight(
    gas: str, equivalent_species: str, definition: EquivalenceDefinition
) -> float:
    """
    Get the weight of a gas in an equivalent species

    I.e. the ratio of its radiative efficiency
    to that of the equivalent species' reference gas.

    Parameters
    ----------
    gas
        Gas of interest

    equivalent_species
        Equivalent species of interest

    definition
        Definition whose radiative efficiencies to use

    Returns
    -------
    :
        Weight of `gas` in `equivalent_species`
    """
    res = definition.radiative_efficiencies
    reference = res[definition.get_reference_gas(equivalent_species)]

    return float((res[gas] / reference).to("dimensionless").m)


def calculate_equivalent(
    equivalent_species: str,
    definition: EquivalenceDefinition,
    concentrations: pd.DataFrame,
    components: Iterable[str] | None = None,
) -> pd.Series[float]:
    """
    Calculate an equivalent species

    Parameters
    ----------
    equivalent_species
        Equivalent species to calculate e.g. `cfc12eq`

    definition
        Definition whose radiative efficiencies to use

    concentrations
        Concentrations of the individual gases, one column per gas

    components
        Gases to include

        If `None`, the gases `definition` includes.

    Returns
    -------
    :
        The equivalent species
    """
    if components is None:
        components = definition.components[equivalent_species]

    return sum(  # type: ignore[return-value]
        concentrations[gas] * get_weight(gas, equivalent_species, definition)
        for gas in components
    )


def check_reproduction(
    equivalent_species: str, dataset: EquivalenceDataset
) -> pd.Series[float]:
    """
    Get how far a dataset's definition is from reproducing its published values

    Parameters
    ----------
    equivalent_species
        Equivalent species of interest

    dataset
        Dataset of interest

    Returns
    -------
    :
        Published minus reproduced values, for each year
    """
    reproduced = calculate_equivalent(
        equivalent_species, dataset.definition, dataset.concentrations
    )

    return (dataset.published[equivalent_species] - reproduced).dropna()


DIFFERENCE_PARTS = (
    "difference in included gases",
    "difference in radiative efficiency",
    "other",
)
"""The parts [decompose_difference][] splits a difference into, in order"""


def decompose_difference(
    equivalent_species: str,
    ours: EquivalenceDataset,
    other: EquivalenceDataset,
) -> pd.DataFrame:
    """
    Split the difference between two datasets' equivalent species into its causes

    The difference is ours minus other.
    It is split up by changing other's calculation into ours one step at a time,
    each step's change being one part of the difference:

    1. difference in included gases:
       other's concentrations and radiative efficiencies,
       our group of gases instead of other's
    1. difference in radiative efficiency:
       other's concentrations, our group of gases,
       our radiative efficiencies instead of other's
    1. other: everything that's left,
       i.e. the total difference minus the first two parts,
       so the three parts always add up to the total difference.
       If both definitions reproduce their published values
       (see [check_reproduction][]), this is exactly
       the difference due to the concentrations of the individual gases.

    The parts depend (a little) on the order of the steps,
    because e.g. how much a gas adds depends on the radiative efficiencies
    it is weighted with.
    We include gases first so that other's radiative efficiencies are used
    to weight the gases only it includes: we don't have radiative efficiencies
    for gases we don't include.

    Parameters
    ----------
    equivalent_species
        Equivalent species of interest

    ours
        Our dataset

    other
        Dataset to compare against

    Returns
    -------
    :
        For each year both datasets cover, each dataset's published values,
        the total difference and each of [DIFFERENCE_PARTS][]
    """
    other_reproduced = calculate_equivalent(
        equivalent_species, other.definition, other.concentrations
    )
    our_gases_other_res = calculate_equivalent(
        equivalent_species,
        other.definition,
        other.concentrations,
        components=ours.definition.components[equivalent_species],
    )
    our_gases_our_res = calculate_equivalent(
        equivalent_species, ours.definition, other.concentrations
    )

    res = pd.DataFrame(
        {
            ours.label: ours.published[equivalent_species],
            other.label: other.published[equivalent_species],
        }
    ).dropna()
    res["total difference"] = res[ours.label] - res[other.label]
    res[DIFFERENCE_PARTS[0]] = our_gases_other_res - other_reproduced
    res[DIFFERENCE_PARTS[1]] = our_gases_our_res - our_gases_other_res
    res[DIFFERENCE_PARTS[2]] = (
        res["total difference"] - res[DIFFERENCE_PARTS[0]] - res[DIFFERENCE_PARTS[1]]
    )

    return res


def get_contributions_by_gas(
    equivalent_species: str,
    ours: EquivalenceDataset,
    other: EquivalenceDataset,
    year: int,
) -> pd.DataFrame:
    """
    Split each part of [decompose_difference][] into each gas' contribution

    Parameters
    ----------
    equivalent_species
        Equivalent species of interest

    ours
        Our dataset

    other
        Dataset to compare against

    year
        Year of interest

    Returns
    -------
    :
        Each gas' contribution to each of [DIFFERENCE_PARTS][] in `year`,
        sorted by the size of its total contribution.
        The "other" part is the contribution due to concentrations,
        so it only sums to the "other" part if both datasets' definitions
        reproduce their published values.
    """
    our_gases = ours.definition.components[equivalent_species]
    other_gases = other.definition.components[equivalent_species]

    rows = {}
    for gas in sorted(set(our_gases) | set(other_gases)):
        other_conc = other.concentrations.loc[year, gas]
        # Weight of the gas in the other dataset's calculation
        other_weight = get_weight(gas, equivalent_species, other.definition)

        if gas in our_gases and gas not in other_gases:
            included = other_conc * other_weight
        elif gas in other_gases and gas not in our_gases:
            included = -other_conc * other_weight
        else:
            included = 0.0

        if gas in our_gases:
            our_weight = get_weight(gas, equivalent_species, ours.definition)
            radiative_efficiency = other_conc * (our_weight - other_weight)
            concentration = (
                ours.concentrations.loc[year, gas] - other_conc
            ) * our_weight
        else:
            radiative_efficiency = 0.0
            concentration = 0.0

        rows[gas] = {
            DIFFERENCE_PARTS[0]: included,
            DIFFERENCE_PARTS[1]: radiative_efficiency,
            DIFFERENCE_PARTS[2]: concentration,
        }

    res = pd.DataFrame(rows).T
    res["total"] = res.sum(axis="columns")

    return res.sort_values("total", key=lambda s: s.abs(), ascending=False)


def get_gases(definitions: Iterable[EquivalenceDefinition]) -> list[str]:
    """
    Get every gas any of the definitions includes in any equivalent species

    Parameters
    ----------
    definitions
        Definitions of interest

    Returns
    -------
    :
        Every gas included, sorted
    """
    return sorted(
        {
            gas
            for definition in definitions
            for components in definition.components.values()
            for gas in components
        }
    )


def load_cmip7(bundle_dir: Path = DEFAULT_BUNDLE_DIR) -> EquivalenceDataset:
    """
    Load our dataset

    Parameters
    ----------
    bundle_dir
        Directory which holds the original run's bundle

    Returns
    -------
    :
        Our dataset
    """

    def load(gas: str) -> pd.Series[float]:
        da = get_only_data_variable(
            xr.load_dataset(
                bundle_dir
                / "data"
                / "interim"
                / gas
                / f"{gas}_global-mean_annual-mean.nc"
            )
        )
        return pd.Series(
            to_units(da.values, da.attrs["units"]),
            index=da["year"].values.astype(int),
            name=gas,
        )

    return EquivalenceDataset(
        definition=CMIP7_DEFINITION,
        concentrations=pd.concat(
            [load(gas) for gas in get_gases([CMIP7_DEFINITION])], axis="columns"
        ),
        published=pd.concat([load(gas) for gas in EQUIVALENT_SPECIES], axis="columns"),
    )


def load_cmip6() -> EquivalenceDataset:
    """
    Load the CMIP6 dataset

    Returns
    -------
    :
        The CMIP6 dataset
    """

    def load(gas: str) -> pd.Series[float]:
        da = get_cmip6_spatial_means(gas, "yr").sel(sector="Global")
        return pd.Series(
            to_units(da.values, da.attrs["units"]),
            # Times are the middle of each year
            index=(da["time"].values - 0.5).round().astype(int),
            name=gas,
        )

    return EquivalenceDataset(
        definition=CMIP6_DEFINITION,
        concentrations=pd.concat(
            [load(gas) for gas in get_gases([CMIP6_DEFINITION])], axis="columns"
        ),
        published=pd.concat([load(gas) for gas in EQUIVALENT_SPECIES], axis="columns"),
    )


IGCC_EQUIVALENT_COLUMNS = {
    "cfc12eq": "CFC[CFC-12-eq]",
    "hfc134aeq": "HFC[HFC-134a-eq]",
}
"""Column of IGCC's concentrations file which holds each of its equivalent species"""


def load_igcc() -> EquivalenceDataset:
    """
    Load the IGCC dataset

    Like in the results figures, only 1850 onwards:
    the 1750 value is IGCC's pre-industrial reference, not part of the record.

    Returns
    -------
    :
        The IGCC dataset
    """
    raw = pd.read_csv(
        ensure_file_downloaded(IGCC_CONCENTRATIONS_URL, IGCC_CONCENTRATIONS_FILE),
        index_col="YYYY",
    )
    first_year = 1850
    raw = raw.loc[raw.index >= first_year]
    raw.index = raw.index.astype(int)

    columns = {**IGCC_COLUMNS, **IGCC_EXTRA_COLUMNS}
    concentrations = pd.DataFrame(
        {
            gas: to_units(
                raw[columns[gas]],
                get_igcc_units(gas) if gas in IGCC_COLUMNS else "ppt",
            )
            # Gases in our definition too, so we can swap our groups of gases
            # into IGCC's calculation
            for gas in get_gases([IGCC_DEFINITION, CMIP7_DEFINITION])
        },
        index=raw.index,
    )
    published = pd.DataFrame(
        {
            gas: to_units(raw[column], "ppt")
            for gas, column in IGCC_EQUIVALENT_COLUMNS.items()
        },
        index=raw.index,
    )

    return EquivalenceDataset(
        definition=IGCC_DEFINITION,
        concentrations=concentrations,
        published=published,
    )


def get_radiative_effect(equivalent_species: str, concentration: float) -> float:
    """
    Get the approximate radiative effect of a concentration of an equivalent species

    Uses our radiative efficiencies, like the results figures.

    Parameters
    ----------
    equivalent_species
        Equivalent species of interest

    concentration
        Concentration, in [UNITS][]

    Returns
    -------
    :
        Approximate radiative effect, in W / m^2
    """
    return float(
        (Q(concentration, UNITS) * RADIATIVE_EFFICIENCIES[equivalent_species])
        .to("W / m^2")
        .m
    )


DIFFERENCE_PART_COLOURS = {
    "total difference": OKABE_ITO["black"],
    DIFFERENCE_PARTS[0]: OKABE_ITO["orange"],
    DIFFERENCE_PARTS[1]: OKABE_ITO["blue"],
    DIFFERENCE_PARTS[2]: OKABE_ITO["bluish green"],
}
"""Colour to draw the total difference and each of its parts in"""


def plot_decomposition(
    equivalent_species: str,
    decomposition: pd.DataFrame,
    ours: EquivalenceDataset,
    other: EquivalenceDataset,
) -> matplotlib.figure.Figure:
    """
    Plot two datasets' equivalent species and their difference, split into its parts

    Parameters
    ----------
    equivalent_species
        Equivalent species of interest

    decomposition
        Output of [decompose_difference][]

    ours
        Our dataset

    other
        Dataset compared against

    Returns
    -------
    :
        The figure
    """
    fig, (ax_values, ax_diff) = plt.subplots(
        nrows=2, sharex=True, figsize=(8, 7), layout="constrained"
    )

    for label, linestyle in ((ours.label, "-"), (other.label, "--")):
        ax_values.plot(
            decomposition.index,
            decomposition[label],
            label=label,
            linestyle=linestyle,
        )

    ax_values.set_ylabel(label_name(f"{equivalent_species} [{UNITS}]"))
    ax_values.legend()

    for column, colour in DIFFERENCE_PART_COLOURS.items():
        ax_diff.plot(
            decomposition.index,
            decomposition[column],
            label=column,
            color=colour,
            linewidth=2.0 if column == "total difference" else 1.0,
        )

    ax_diff.axhline(0.0, color="0.6", linewidth=0.5, zorder=0)
    ax_diff.set_ylabel(f"{ours.label} - {other.label} [{UNITS}]")
    ax_diff.set_xlabel("year")
    ax_diff.legend()

    return fig
