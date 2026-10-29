"""
Summary of the original run's MAGICC output

The original run's MAGICC output database holds every ensemble member
(600 per scenario and MAGICC version), which is GBs.
We only need a few percentiles of a few variables,
so we summarise it once and use the summary from then on.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from loguru import logger

MAGICC_VERSION = "MAGICCv7.6.0a3"
"""Version of MAGICC the concentrations are based on"""

MAGICC_AR6_VERSION = "MAGICCv7.5.3"
"""Version of MAGICC used in AR6 (also run, but not used for the concentrations)"""

MAGICC_SCENARIO_NAMES = {
    "SSP1 - Very Low Emissions": "vl",
    "SSP2 - Low Overshoot_a": "ln",
    "SSP2 - Low Emissions": "l",
    "SSP2 - Medium-Low Emissions": "ml",
    "SSP2 - Medium Emissions": "m",
    "SSP5 - Medium-Low Emissions_a": "hl",
    "SSP3 - High Emissions": "h",
    "VL-cf": "vl-cf",
}
"""CMIP scenario name of each scenario in the MAGICC output"""

MAGICC_GAS_NAMES = {
    "c2f6": "C2F6",
    "c3f8": "C3F8",
    "c4f10": "C4F10",
    "c5f12": "C5F12",
    "c6f14": "C6F14",
    "c7f16": "C7F16",
    "c8f18": "C8F18",
    "cc4f8": "cC4F8",
    "ccl4": "CCl4",
    "cf4": "CF4",
    "cfc11": "CFC11",
    "cfc113": "CFC113",
    "cfc114": "CFC114",
    "cfc115": "CFC115",
    "cfc12": "CFC12",
    "ch2cl2": "CH2Cl2",
    "ch3br": "CH3Br",
    "ch3ccl3": "CH3CCl3",
    "ch3cl": "CH3Cl",
    "ch4": "CH4",
    "chcl3": "CHCl3",
    "co2": "CO2",
    "halon1211": "Halon1211",
    "halon1301": "Halon1301",
    "halon2402": "Halon2402",
    "hcfc141b": "HCFC141b",
    "hcfc142b": "HCFC142b",
    "hcfc22": "HCFC22",
    "hfc125": "HFC125",
    "hfc134a": "HFC134a",
    "hfc143a": "HFC143a",
    "hfc152a": "HFC152a",
    "hfc227ea": "HFC227ea",
    "hfc23": "HFC23",
    "hfc236fa": "HFC236fa",
    "hfc245fa": "HFC245fa",
    "hfc32": "HFC32",
    "hfc365mfc": "HFC365mfc",
    "hfc4310mee": "HFC4310mee",
    "n2o": "N2O",
    "nf3": "NF3",
    "sf6": "SF6",
    "so2f2": "SO2F2",
}
"""Name of each gas in MAGICC's output (and in the emissions scenarios)"""

NPP_VARIABLE = "CO2_CURRENT_NPP"
"""MAGICC's net primary productivity (NPP) variable"""

QUANTILES = (0.05, 0.5, 0.95)
"""Quantiles of the ensemble we keep"""


def get_concentration_variable(gas: str) -> str:
    """Get the MAGICC variable which holds a gas' concentrations"""
    return f"Atmospheric Concentrations|{MAGICC_GAS_NAMES[gas]}"


def summarise_magicc_output(db_dir: Path, out_file: Path) -> Path:
    """
    Summarise the MAGICC output database

    Parameters
    ----------
    db_dir
        Directory of the MAGICC output database
        (feather backend, as written by `pandas_openscm`)

    out_file
        File in which to write the summary

    Returns
    -------
    :
        `out_file`
    """
    variables = {
        NPP_VARIABLE,
        *(get_concentration_variable(g) for g in MAGICC_GAS_NAMES),
    }

    index = pd.read_feather(db_dir / "index.feather").reset_index()
    index = index[index["variable"].isin(variables)]
    missing = variables - set(index["variable"])
    if missing:
        msg = f"Not in the MAGICC output: {sorted(missing)}"
        raise AssertionError(msg)

    group_levels = ["climate_model", "scenario", "variable", "unit"]
    res_l = []
    file_ids = index["file_id"].unique()
    logger.info(f"Summarising {len(file_ids)} files from {db_dir}")
    for file_id in file_ids:
        raw = pd.read_feather(db_dir / f"{file_id}.feather")
        raw = raw[raw.index.get_level_values("variable").isin(variables)]
        quantiles = raw.groupby(level=group_levels).quantile(list(QUANTILES))
        quantiles.index = quantiles.index.set_names("quantile", level=-1)
        res_l.append(quantiles)

    res = pd.concat(res_l).sort_index()
    res.columns = res.columns.astype(int)
    res.columns.name = "year"
    # Switch to CMIP scenario names, which is how we refer to scenarios everywhere else
    res = res.rename(index=MAGICC_SCENARIO_NAMES, level="scenario")
    if res.index.duplicated().any():
        raise AssertionError

    out_file.parent.mkdir(parents=True, exist_ok=True)
    res.to_feather(out_file)
    logger.info(f"Wrote {out_file}")

    return out_file


def get_magicc_quantile(
    summary: pd.DataFrame,
    variable: str,
    scenario: str,
    quantile: float,
    climate_model: str = MAGICC_VERSION,
) -> pd.Series:
    """
    Get one quantile of one variable in one scenario from the summary

    Parameters
    ----------
    summary
        Summary, as written by [summarise_magicc_output][]

    variable
        MAGICC variable

    scenario
        CMIP scenario name

    quantile
        Quantile, one of [QUANTILES][]

    climate_model
        MAGICC version

    Returns
    -------
    :
        Values, indexed by year, with the units as the series' name
    """
    selected = summary.xs(
        (climate_model, scenario, variable, quantile),
        level=("climate_model", "scenario", "variable", "quantile"),
    )
    if selected.shape[0] != 1:
        raise AssertionError(selected)

    units = selected.index.get_level_values("unit")[0]

    return pd.Series(selected.iloc[0].to_numpy(), index=selected.columns, name=units)
