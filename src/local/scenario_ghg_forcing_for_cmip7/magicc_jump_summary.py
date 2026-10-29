"""
Comparing the MAGICC runs with and without the 2014/2015 jump in MAGICC's inputs

The runs are described (and run) in
[local.scenario_ghg_forcing_for_cmip7.magicc_jump][].
Here we check that our `published` runs reproduce the published output
and compare the two set-ups' ensemble medians.
The result is summarised in [MAGICC_JUMP_SUMMARY][],
so the value checks don't have to load the runs.
"""

from __future__ import annotations

import json
from functools import cache
from pathlib import Path

import numpy as np
import pandas as pd

from local.scenario_ghg_forcing_for_cmip7.bundle import BUNDLE_DIR
from local.scenario_ghg_forcing_for_cmip7.magicc_jump import (
    EXPERIMENTS,
    MAGICC_JUMP_DIR,
    SETUPS,
    get_out_file,
)
from local.scenario_ghg_forcing_for_cmip7.magicc_output import (
    MAGICC_SCENARIO_NAMES,
    MAGICC_VERSION,
)
from local.scenario_ghg_forcing_for_cmip7.rerun_projections import MAGICC_OUTPUT_DB_DIR

REPRODUCTION_VARIABLES = (
    "Surface Air Temperature Change",
    "Effective Radiative Forcing",
    "Atmospheric Concentrations|CO2",
    "Atmospheric Concentrations|CH4",
    "Atmospheric Concentrations|N2O",
)
"""Variables which must match the published output for a run to count as reproduced

If these match, we can safely assume the rest does too.
"""

JUMP_VARIABLES = {
    "Surface Air Temperature Change": "gsat",
    "Effective Radiative Forcing": "erf",
    "Effective Radiative Forcing|Greenhouse Gases": "erf-ghg",
    "Effective Radiative Forcing|Aerosols": "erf-aerosols",
    "Atmospheric Concentrations|CO2": "co2",
    "Atmospheric Concentrations|CH4": "ch4",
    "Atmospheric Concentrations|N2O": "n2o",
}
"""Variables whose difference between the set-ups we report, with their short names

The GHG and aerosol ERF show which of the two drives the difference in total ERF.
The concentrations are MAGICC's (i.e. before we harmonise them to the history).
"""

JUMP_PERIODS = {
    "history": (1750, 2014),
    "to-2100": (2015, 2100),
    "to-2500": (2015, 2500),
}
"""Periods over which we report the difference between the set-ups, with their names

Over the history, the set-ups differ in their inputs
(CMIP6 vs. CMIP7 concentrations and emissions).
From 2015 (the first year of the scenarios), their inputs are the same,
so the difference is the effect of the jump.
"""

MAGICC_JUMP_SUMMARY = MAGICC_JUMP_DIR / "summary.json"
"""Summary of the runs, relative to the bundle's root (see [summarise_magicc_jump][])"""


@cache
def load_run(name: str, setup: str, bundle_dir: Path = BUNDLE_DIR) -> pd.DataFrame:
    """
    Load one experiment's run in one set-up

    Parameters
    ----------
    name
        CMIP scenario name of the experiment

    setup
        Set-up, one of [local.scenario_ghg_forcing_for_cmip7.magicc_jump.SETUPS][]

    bundle_dir
        Directory which holds the bundle

    Returns
    -------
    :
        Every ensemble member, one row per member and variable, one column per year
    """
    if setup not in SETUPS:
        raise ValueError(setup)

    raw = pd.read_feather(bundle_dir / get_out_file(name, setup))
    year_columns = [c for c in raw.columns if c.isdigit()]
    res = raw.set_index([c for c in raw.columns if c not in year_columns])
    res.columns = res.columns.astype(int)

    return res


def load_published(
    name: str, variables: tuple[str, ...], bundle_dir: Path = BUNDLE_DIR
) -> pd.DataFrame:
    """
    Load the published MAGICC output for one experiment

    Parameters
    ----------
    name
        CMIP scenario name of the experiment

    variables
        Variables to load

    bundle_dir
        Directory which holds the bundle

    Returns
    -------
    :
        Every ensemble member, one row per member and variable, one column per year
    """
    db_dir = bundle_dir / MAGICC_OUTPUT_DB_DIR
    index = pd.read_feather(db_dir / "index.feather").reset_index()
    index = index[
        (index["climate_model"] == MAGICC_VERSION)
        & (index["scenario"].map(MAGICC_SCENARIO_NAMES) == name)
        & (index["variable"].isin(variables))
    ]
    res = pd.concat(
        pd.read_feather(db_dir / f"{file_id}.feather")
        for file_id in index["file_id"].unique()
    )
    res.columns = res.columns.astype(int)

    return res


def get_reproduction_differences(
    name: str, bundle_dir: Path = BUNDLE_DIR
) -> pd.Series[float]:
    """
    Get how far our `published` run is from the published output

    Compared member by member, for [REPRODUCTION_VARIABLES][].

    Parameters
    ----------
    name
        CMIP scenario name of the experiment

    bundle_dir
        Directory which holds the bundle

    Returns
    -------
    :
        Largest absolute difference, by variable (in the variable's units)
    """
    ours = load_run(name, "published", bundle_dir)
    ours = ours[ours.index.get_level_values("variable").isin(REPRODUCTION_VARIABLES)]
    published = load_published(name, REPRODUCTION_VARIABLES, bundle_dir)

    def key(df: pd.DataFrame) -> pd.DataFrame:
        return df.reset_index().set_index(["variable", "run_id"])[list(df.columns)]

    ours_k = key(ours)
    published_k = key(published).loc[ours_k.index, ours_k.columns]
    diff = (ours_k - published_k).abs().max(axis="columns")

    return diff.groupby(level="variable").max()


def get_median(name: str, setup: str, variable: str, bundle_dir: Path) -> pd.Series:
    """
    Get the ensemble median of one variable for one experiment in one set-up
    """
    run = load_run(name, setup, bundle_dir)
    run = run[run.index.get_level_values("variable") == variable]

    return run.median(axis="index")


def get_median_difference(name: str, variable: str, bundle_dir: Path) -> pd.Series:
    """
    Get the difference between the set-ups' ensemble medians

    `cmip7-history` minus `published`.
    """
    return get_median(name, "cmip7-history", variable, bundle_dir) - get_median(
        name, "published", variable, bundle_dir
    )


def get_max_median_difference(
    variable: str, period: tuple[int, int], bundle_dir: Path = BUNDLE_DIR
) -> dict[str, float | str | int]:
    """
    Get the largest difference between the two set-ups' ensemble medians

    Over every experiment and the years in `period`.

    Parameters
    ----------
    variable
        Variable of interest

    period
        First and last year (inclusive) to consider

    bundle_dir
        Directory which holds the bundle

    Returns
    -------
    :
        Largest absolute difference (`value`, in the variable's units),
        the difference with its sign (`signed`),
        and the experiment (`experiment`) and year (`year`) in which it occurs
    """
    largest: dict[str, float | str | int] = {"value": -np.inf}
    for _, _, name in EXPERIMENTS:
        diff = get_median_difference(name, variable, bundle_dir).loc[
            period[0] : period[1]
        ]
        year = int(diff.abs().idxmax())
        if abs(diff.loc[year]) > largest["value"]:
            largest = {
                "value": float(abs(diff.loc[year])),
                "signed": float(diff.loc[year]),
                "experiment": name,
                "year": year,
            }

    return largest


def summarise_magicc_jump(bundle_dir: Path = BUNDLE_DIR) -> Path:
    """
    Summarise the runs, so the value checks don't have to load them

    The summary holds

    - the reproduction check for every experiment
      (see [get_reproduction_differences][])
    - the largest difference between the set-ups' medians for [JUMP_VARIABLES][]
      over each of [JUMP_PERIODS][] (see [get_max_median_difference][])
    - the difference in the GHG and aerosol ERF where the difference in total ERF
      is largest, for each of [JUMP_PERIODS][]

    Parameters
    ----------
    bundle_dir
        Directory which holds the bundle

    Returns
    -------
    :
        The summary file
    """
    summary = {
        "n_runs": int(
            load_run(EXPERIMENTS[0][-1], "published", bundle_dir)
            .index.get_level_values("run_id")
            .nunique()
        ),
        "reproduction": {
            name: {
                k: float(v)
                for k, v in get_reproduction_differences(name, bundle_dir).items()
            }
            for _, _, name in EXPERIMENTS
        },
        "max_median_difference": {
            short_name: {
                "variable": variable,
                **{
                    period_name: get_max_median_difference(variable, period, bundle_dir)
                    for period_name, period in JUMP_PERIODS.items()
                },
            }
            for variable, short_name in JUMP_VARIABLES.items()
        },
        "erf_breakdown_at_max_erf_difference": {},
    }
    for period_name in JUMP_PERIODS:
        largest = summary["max_median_difference"]["erf"][period_name]
        summary["erf_breakdown_at_max_erf_difference"][period_name] = {
            "experiment": largest["experiment"],
            "year": largest["year"],
            **{
                short_name: float(
                    get_median_difference(
                        largest["experiment"], variable, bundle_dir
                    ).loc[largest["year"]]
                )
                for variable, short_name in JUMP_VARIABLES.items()
                if short_name.startswith("erf")
            },
        }

    out_file = bundle_dir / MAGICC_JUMP_SUMMARY
    out_file.parent.mkdir(exist_ok=True, parents=True)
    out_file.write_text(json.dumps(summary, indent=2))

    return out_file
