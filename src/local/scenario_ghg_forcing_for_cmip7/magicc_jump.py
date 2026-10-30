"""
Running MAGICC to quantify the effect of the 2014/2015 jump in MAGICC's inputs

The published MAGICC runs use the CMIP6 historical concentrations
(and MAGICC's own CMIP6-based historical emissions) up to 2014,
then each experiment's CMIP7 emissions from 2015
(see the methods, Step 1 for CO2).
To quantify the effect of this jump, we run MAGICC for every experiment in two set-ups
(see the notebook, [MAGICC_JUMP_NOTEBOOK][]):

- `published`: as the original run, which also checks that we reproduce
  the published output
- `cmip7-history`: with the CMIP7 historical concentrations and emissions instead,
  so there is no jump

The difference between the two is the largest effect we could expect from the jump
(had we recalibrated MAGICC on the CMIP7 history, it would be smaller).

Everything is written below [MAGICC_JUMP_DIR][] in the bundle.
The runs are compared in [local.scenario_ghg_forcing_for_cmip7.magicc_jump_summary][]
(kept separate so that changing the comparison doesn't re-run MAGICC).
"""

from __future__ import annotations

from pathlib import Path
from tempfile import TemporaryDirectory

import jupytext
import yaml

from local.cmip_ghg_generation import ensure_bundle_available, run_in_bundle
from local.scenario_ghg_forcing_for_cmip7.bundle import (
    BUNDLE_DIR,
    MODIFIED_NOTEBOOKS_DIR,
    RERUN_NOTEBOOKS_DIR,
    TARBALL_STRIP_COMPONENTS,
    ZENODO_RECORD_ID,
)
from local.scenario_ghg_forcing_for_cmip7.rerun_projections import (
    COMPLETE_EMISSIONS_DIR,
    HISTORICAL_CONCENTRATIONS_ROOT,
    MARKERS,
)
from local.scenario_ghg_forcing_for_cmip7.rerun_vl_cf import VL_CF_MARKER

EXPERIMENTS = (*MARKERS, VL_CF_MARKER)
"""Experiments we run, as (model, scenario, CMIP scenario name)"""

SETUPS = ("published", "cmip7-history")
"""The two set-ups we run each experiment in, see the module's docstring"""

N_RUNS = 0
"""Number of ensemble members to run (0 means all 600)

Set to a small number to try things out quickly.
"""

N_MAGICC_WORKERS = 6
"""Number of MAGICC processes to run in parallel"""

MAGICC_TARBALL = "magicc--magicc-v7.6.0a3.tar.gz"
"""The bundle's tarball of the MAGICC version the concentrations are based on"""

MAGICC_DIR = Path("magicc") / "magicc-v7.6.0a3"
"""Where the MAGICC tarball extracts to, relative to the bundle's root"""

MAGICC_EXE = MAGICC_DIR / "bin" / "magicc-darwin-arm64"
"""MAGICC executable, relative to the bundle's root (as in the original run)"""

MAGICC_PROB_DISTRIBUTION = (
    MAGICC_DIR / "configs" / "magicc-ar7-fast-track-drawnset-v0-3-0.json"
)
"""MAGICC's probabilistic configuration, relative to the bundle's root"""

MAGICC_JUMP_NOTEBOOK = MODIFIED_NOTEBOOKS_DIR / "0021_run-magicc_cmip7-history-check.py"
"""The notebook which runs one experiment in one set-up

Written by us (from the original run's `0021` and `1200` notebooks), tracked in git.
"""

MAGICC_JUMP_DIR = Path("manuscript-outputs") / "magicc-jump"
"""Where the runs are written, relative to the bundle's root"""


def get_complete_emissions_file(model: str, scenario: str) -> Path:
    """
    Get an experiment's complete emissions file, relative to the bundle's root

    Parameters
    ----------
    model
        Model which produced the experiment's emissions

    scenario
        Scenario name used by the model

    Returns
    -------
    :
        The file (1750 to 2500)
    """
    stem = f"{scenario}_{model}".replace(" ", "_").replace(".", "-")

    return COMPLETE_EMISSIONS_DIR / f"{stem}.feather"


def get_out_file(name: str, setup: str) -> Path:
    """
    Get the file which holds one experiment's run in one set-up

    Parameters
    ----------
    name
        CMIP scenario name of the experiment

    setup
        Set-up, one of [SETUPS][]

    Returns
    -------
    :
        The file, relative to the bundle's root
    """
    return MAGICC_JUMP_DIR / f"{name}_{setup}.feather"


def ensure_magicc(bundle_dir: Path = BUNDLE_DIR) -> Path:
    """
    Ensure the bundle's copy of MAGICC is available

    Parameters
    ----------
    bundle_dir
        Directory which holds the bundle

    Returns
    -------
    :
        The MAGICC executable
    """
    ensure_bundle_available(
        files_to_get=(),
        files_to_get_tarred=(MAGICC_TARBALL,),
        bundle_dir=bundle_dir,
        record_id=ZENODO_RECORD_ID,
        tar_strip_components=TARBALL_STRIP_COMPONENTS,
    )

    return bundle_dir / MAGICC_EXE


def run_magicc_jump(  # noqa: PLR0913
    model: str,
    scenario: str,
    name: str,
    setup: str,
    n_runs: int = N_RUNS,
    bundle_dir: Path = BUNDLE_DIR,
) -> Path:
    """
    Run MAGICC for one experiment in one set-up

    Runs [MAGICC_JUMP_NOTEBOOK][] in the bundle's own (pixi) environment.

    Parameters
    ----------
    model
        Model which produced the experiment's emissions

    scenario
        Scenario name used by the model

    name
        CMIP scenario name of the experiment

    setup
        Set-up, one of [SETUPS][]

    n_runs
        Number of ensemble members to run (0 means all)

    bundle_dir
        Directory which holds the bundle

    Returns
    -------
    :
        The file which holds the run's output
    """
    if setup not in SETUPS:
        raise ValueError(setup)

    ensure_magicc(bundle_dir)

    run_name = f"0021_run-magicc_cmip7-history-check_{name}_{setup}"
    to_run = bundle_dir / RERUN_NOTEBOOKS_DIR / f"{run_name}.ipynb"
    executed = bundle_dir / RERUN_NOTEBOOKS_DIR / f"{run_name}_executed.ipynb"
    to_run.parent.mkdir(exist_ok=True, parents=True)
    jupytext.write(jupytext.read(MAGICC_JUMP_NOTEBOOK), to_run, fmt="ipynb")

    out_file = get_out_file(name, setup)
    parameters = {
        "model": model,
        "scenario": scenario,
        "complete_file": get_complete_emissions_file(model, scenario).as_posix(),
        "setup": setup,
        "n_runs": n_runs,
        "magicc_exe": MAGICC_EXE.as_posix(),
        "magicc_prob_distribution": MAGICC_PROB_DISTRIBUTION.as_posix(),
        "n_magicc_workers": N_MAGICC_WORKERS,
        "historical_concentrations_root": (
            HISTORICAL_CONCENTRATIONS_ROOT.resolve().as_posix()
        ),
        "out_file": out_file.as_posix(),
    }
    # MAGICC's Fortran filename fields are limited to 255 characters. A deep
    # checkout plus a long model/scenario name can silently truncate SCEN7 paths.
    # A short-lived alias keeps exactly the same executable and input files.
    with TemporaryDirectory(prefix="cmip7-magicc-", dir="/tmp") as temporary:
        alias = Path(temporary) / "model"
        alias.symlink_to((bundle_dir / MAGICC_DIR).resolve(), target_is_directory=True)
        parameters["magicc_exe"] = (alias / "bin" / MAGICC_EXE.name).as_posix()
        run_in_bundle(
            [
                "pixi",
                "run",
                "--frozen",
                "papermill",
                "--cwd",
                str(bundle_dir.resolve()),
                str(to_run.resolve()),
                str(executed.resolve()),
                "-y",
                yaml.safe_dump(parameters),
            ],
            bundle_dir=bundle_dir,
        )

    return bundle_dir / out_file
