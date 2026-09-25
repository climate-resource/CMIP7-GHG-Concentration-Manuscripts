"""
Generation of the tex inputs which summarise the gases processed like CFC-12

These gases all follow the same method,
but each takes its own inputs and choices
(which networks observe it, where its global-mean comes from,
its pre-industrial value etc.).
Rather than write thirty-four sets of these into the text,
we pull them out of the original run's bundle
and write them into tex which the manuscript's text refers to.

The gas names are written the way the manuscript's source text writes them
(e.g. `C_2F_6`, `CFC12`), so the manuscript's replacements
format them exactly as they are formatted in the text around them.
"""

from __future__ import annotations

from pathlib import Path

import xarray as xr
from loguru import logger

from local.historical_ghg_forcing_for_cmip7.cfc12_like_methods_figure import (
    CFC12_LIKE_GASES,
    SOURCE_BIBKEYS,
    get_cfc12_like_all_data_with_bins,
    get_global_mean_supplement,
    get_step_config,
    interim_dir,
    supplement_replaces_obs_network,
)

MANUSCRIPT_GAS_NAMES = {
    "c2f6": "C_2F_6",
    "c3f8": "C_3F_8",
    "ccl4": "CCl_4",
    "cf4": "CF_4",
    "cfc11": "CFC-11",
    "cfc113": "CFC-113",
    "cfc114": "CFC-114",
    "cfc115": "CFC-115",
    "cfc12": "CFC12",
    "ch2cl2": "CH_2Cl_2",
    "ch3br": "CH_3Br",
    "ch3ccl3": "CH_3CCl_3",
    "ch3cl": "CH_3Cl",
    "chcl3": "CHCl_3",
    "halon1211": "Halon-1211",
    "halon1301": "Halon-1301",
    "halon2402": "Halon-2402",
    "hcfc141b": "HCFC-141b",
    "hcfc142b": "HCFC-142b",
    "hcfc22": "HCFC-22",
    "hfc125": "HFC-125",
    "hfc134a": "HFC-134a",
    "hfc143a": "HFC-143a",
    "hfc152a": "HFC-152a",
    "hfc227ea": "HFC-227ea",
    "hfc23": "HFC-23",
    "hfc236fa": "HFC-236fa",
    "hfc245fa": "HFC-245fa",
    "hfc32": "HFC-32",
    "hfc365mfc": "HFC-365mfc",
    "hfc4310mee": "HFC-4310mee",
    "nf3": "NF_3",
    "sf6": "SF_6",
    "so2f2": "SO_2F_2",
}
"""How the manuscript's source text writes each gas' name

These are keys of the manuscript's `replacements.yaml`,
which is where they are turned into tex.
"""

NETWORK_LABELS = {
    "AGAGE": "AGAGE",
    "NOAA": "NOAA HATS",
}
"""How each observational network is written in the manuscript

Keyed by the network's name in the `network` column
of the original run's binned observational network data.
"""

PRE_INDUSTRIAL_SOURCE_CITATIONS = {
    "M17": r"\citet{meinshausen_historical_2017}",
    "Velders et al., 2022": r"\citet{velders_2022}",
    "Velders et al., 2022 (with adjustments to support interpolation)": (
        r"\citet{velders_2022}\textsuperscript{a}"
    ),
}
"""Citation for each pre-industrial source, as the original run's config names it

The adjusted source is marked with a footnote, which the table explains.
The adjustment is HFC-134a's pre-industrial year, which the original run
moved from 1980 to 1988 because, with 1980, every gap-filling method
in `1304_sf6-like_create-global-annual-mean` fails
(the fits overshoot or go negative).
"""

PRE_INDUSTRIAL_UNIT = "ppt"
"""Unit the pre-industrial values are given in

Every gas in this group uses it, so it goes in the column heading.
"""


def get_networks(gas: str, bundle_dir: Path) -> tuple[str, ...]:
    """
    Get the observational networks which observe a gas

    Parameters
    ----------
    gas
        Gas of interest

    bundle_dir
        Directory which holds the original run's bundle

    Returns
    -------
    :
        The networks' names, as the manuscript writes them, sorted

    Raises
    ------
    AssertionError
        The NOAA data is not all from HATS,
        in which case the label we use for NOAA would be wrong
    """
    all_data_with_bins = get_cfc12_like_all_data_with_bins(gas, bundle_dir)

    noaa_sources = set(
        all_data_with_bins.loc[all_data_with_bins["network"] == "NOAA", "source"]
    )
    if noaa_sources - {"hats"}:
        msg = f"Expected only NOAA HATS data for {gas=}, found {noaa_sources}"
        raise AssertionError(msg)

    return tuple(
        sorted(NETWORK_LABELS[v] for v in all_data_with_bins["network"].unique())
    )


def generate_cfc12_like_obs_network_sources_list(
    outfile: Path,
    bundle_dir: Path,
    force_rerun: bool = False,
) -> Path:
    """
    Generate the list of which networks observe which gases processed like CFC-12

    Parameters
    ----------
    outfile
        File in which to write the list

    bundle_dir
        Directory which holds the original run's bundle

    force_rerun
        Re-generate the list, even if the output file already exists

    Returns
    -------
    :
        `outfile`
    """
    if outfile.exists() and not force_rerun:
        logger.info(f"Using existing {outfile}")
        return outfile

    gases_by_networks: dict[tuple[str, ...], list[str]] = {}
    for gas in CFC12_LIKE_GASES:
        gases_by_networks.setdefault(get_networks(gas, bundle_dir), []).append(gas)

    items = []
    # Most networks first, which also puts the group most gases are in first
    for networks, gases in sorted(
        gases_by_networks.items(), key=lambda kv: (-len(kv[0]), kv[0])
    ):
        if len(networks) == 1:
            networks_label = f"{networks[0]} only"
        else:
            networks_label = " and ".join(networks)

        gas_names = ", ".join(MANUSCRIPT_GAS_NAMES[gas] for gas in gases)
        items.append(f"    \\item {networks_label}: {gas_names}")

    # Gas names can't be broken at their hyphens,
    # so a justified line has too few places to stretch
    res = "\n".join(["\\begin{itemize}", "\\raggedright", *items, "\\end{itemize}", ""])

    outfile.parent.mkdir(exist_ok=True, parents=True)
    logger.info(f"Writing {outfile}")
    outfile.write_text(res)

    return outfile


def get_global_mean_source_cell(gas: str, bundle_dir: Path) -> str:
    """
    Get the table cell which says where a gas' global-mean comes from

    Parameters
    ----------
    gas
        Gas of interest

    bundle_dir
        Directory which holds the original run's bundle

    Returns
    -------
    :
        Table cell
    """
    global_mean_supplement = get_global_mean_supplement(gas, bundle_dir)
    if global_mean_supplement is None:
        return "Obs. network"

    label, supplement = global_mean_supplement
    citation = f"\\citet{{{SOURCE_BIBKEYS[label]}}}"

    global_mean_extended = xr.load_dataset(
        interim_dir(gas, bundle_dir) / f"{gas}_global-annual-mean_allyears.nc"
    )
    if supplement_replaces_obs_network(
        supplement, int(global_mean_extended["year"].max())
    ):
        return citation

    return f"{citation}\\textsuperscript{{b}}"


def get_per_gas_row(gas: str, bundle_dir: Path) -> str:
    """
    Get a gas' row in the per-gas table

    Parameters
    ----------
    gas
        Gas of interest

    bundle_dir
        Directory which holds the original run's bundle

    Returns
    -------
    :
        The row, as tex

    Raises
    ------
    AssertionError
        The gas' pre-industrial value is not in [`PRE_INDUSTRIAL_UNIT`][]
    """
    global_mean_from_obs_network = xr.load_dataset(
        interim_dir(gas, bundle_dir)
        / f"{gas}_observational-network_global-annual-mean.nc"
    )
    obs_network_years = (
        f"{int(global_mean_from_obs_network['year'].min())}--"
        f"{int(global_mean_from_obs_network['year'].max())}"
    )

    step_config = get_step_config(gas, bundle_dir)
    long_poleward_extension = (
        "Yes" if step_config["allow_long_poleward_extension"] else "No"
    )

    pre_industrial = step_config["pre_industrial"]
    pre_industrial_value, pre_industrial_unit = pre_industrial["value"]
    if pre_industrial_unit != PRE_INDUSTRIAL_UNIT:
        msg = (
            f"Expected {gas=}'s pre-industrial value in {PRE_INDUSTRIAL_UNIT}, "
            f"got {pre_industrial_unit}"
        )
        raise AssertionError(msg)

    cells = (
        MANUSCRIPT_GAS_NAMES[gas],
        obs_network_years,
        long_poleward_extension,
        get_global_mean_source_cell(gas, bundle_dir),
        str(pre_industrial["year"]),
        f"{pre_industrial_value:g}",
        PRE_INDUSTRIAL_SOURCE_CITATIONS[pre_industrial["source"]],
    )

    return " & ".join(cells) + r" \\"


def generate_cfc12_like_per_gas_table(
    outfile: Path,
    bundle_dir: Path,
    force_rerun: bool = False,
) -> Path:
    """
    Generate the table of per-gas inputs and choices for the gases processed like CFC-12

    Parameters
    ----------
    outfile
        File in which to write the table

    bundle_dir
        Directory which holds the original run's bundle

    force_rerun
        Re-generate the table, even if the output file already exists

    Returns
    -------
    :
        `outfile`
    """
    if outfile.exists() and not force_rerun:
        logger.info(f"Using existing {outfile}")
        return outfile

    rows = [f"    {get_per_gas_row(gas, bundle_dir)}" for gas in CFC12_LIKE_GASES]

    # The three pre-industrial columns share a heading, which saves the width
    # that writing "Pre-industrial" over each of them would take up.
    header = " & ".join(
        (
            "Gas",
            "Obs. network",
            "Long poleward",
            "Global-mean",
            r"\multicolumn{3}{c}{Pre-industrial}",
        )
    )
    sub_header = " & ".join(
        (
            "",
            "years",
            "extension",
            "source",
            "year",
            f"value ({PRE_INDUSTRIAL_UNIT})",
            "source",
        )
    )

    res = "\n".join(
        [
            r"\begin{table*}[t]",
            r"\caption{",
            "    Per-gas inputs and choices for the gases processed like CFC12",
            r"    (Section \ref{ssec:methods-cfc12-like}).",
            "    Obs. network years: years covered by the global-, annual-mean",
            "    derived from the observational network.",
            "    Global-mean source: replaces the global-, annual-mean",
            "    derived from the observational network, unless noted.",
            "    Pre-industrial: the year from which, and value at which,",
            "    the global-, annual-mean is held constant, and their source.",
            "}",
            r"\label{tab:cfc12-like-per-gas}",
            # Thirty-four rows only fit on a page at this size
            r"\footnotesize",
            r"\begin{tabular}{lcclccl}",
            r"\tophline",
            f"{header} \\\\",
            f"{sub_header} \\\\",
            r"\middlehline",
            *rows,
            r"\bottomhline",
            r"\end{tabular}",
            r"\belowtable{",
            "    \\textsuperscript{a} 1988 rather than the 1980 used for",
            r"    the other \citet{velders_2022} pre-industrial years,",
            "    as from 1980, no gap-filling method (Step 4)",
            "    can follow the rapid rise at the start of their data.",
            "    \\textsuperscript{b} Harmonised to the observational network",
            "    in the first year of the observational network",
            "    (with a 100-year transition),",
            "    then the observational network is used from that year on.",
            "}",
            r"\end{table*}",
            "",
        ]
    )

    outfile.parent.mkdir(exist_ok=True, parents=True)
    logger.info(f"Writing {outfile}")
    outfile.write_text(res)

    return outfile
