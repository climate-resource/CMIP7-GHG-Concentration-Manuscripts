"""
Re-running notebooks from the original CMIP7 ScenarioMIP GHG concentration run

The run is archived at [https://zenodo.org/records/21501391]().
As for the historical run (see [local.cmip_ghg_generation][]),
the archive doesn't contain every intermediate object the notebooks built,
nor the notebooks as executed,
so we keep a copy of the executed notebooks we need in this repository
and re-run them with an extra save cell bolted on.

Unlike the historical run's notebooks, which take a config file,
the scenario run's notebooks take their inputs and outputs directly as parameters.
"""

from __future__ import annotations

from pathlib import Path

import openscm_units
import pandas as pd
import pint

from local.cmip_ghg_generation import (
    ensure_bundle_available,
    ensure_bundle_environment,
    run_notebook_from_bundle_dir,
    write_modified_notebook,
)
from local.paths import DATA_RAW_DIR

ZENODO_RECORD_ID = "21501391"
"""ID of the Zenodo record which holds the original run (v1.1.0)"""

TARBALL_STRIP_COMPONENTS = 0
"""Leading path components to strip when extracting the tarballs

Unlike the historical run's, the scenario run's tarballs have no prefix.
"""

BUNDLE_DIR = DATA_RAW_DIR / "cmip7-scenariomip-ghg-concentrations" / "v1.1.0"
"""Where we put our copy of the original run's bundle

Not tracked by git, it is all re-downloadable.
"""

BUNDLE_FILES = ("pyproject.toml", "pixi.lock", "README.md")
"""Files we need from the bundle

The README is needed because `pyproject.toml` names it as the package's readme,
so pixi can't install the bundle's package without it.
"""

BUNDLE_TARBALLS = ("src.tar.gz", "data--interim--monthly-means.tar.gz")
"""Tarballs we need from the bundle

We deliberately avoid the MAGICC tarballs:
they are big (GBs) and none of the notebooks we re-run need them.
"""

EXECUTED_NOTEBOOKS_DIR = (
    DATA_RAW_DIR / "scenario-ghg-forcing-for-cmip7" / "notebooks-executed"
)
"""Where we keep our copy of the original run's notebooks, as executed

Tracked by git, because these notebooks aren't on Zenodo.
"""

MODIFIED_NOTEBOOKS_DIR = (
    DATA_RAW_DIR / "scenario-ghg-forcing-for-cmip7" / "notebooks-modified"
)
"""Where we keep the modified notebooks, in jupytext's `py:percent` format

Tracked by git, so we can see how they change over time.
"""

RERUN_NOTEBOOKS_DIR = Path("notebooks-rerun")
"""Where we write the notebooks we re-run, before and after execution

Relative to the bundle's root directory.
"""

INVERSE_EMISSIONS_GHG = "ccl4"
"""Gas for which we re-run the inverse emissions notebook"""

INVERSE_EMISSIONS_NOTEBOOK = EXECUTED_NOTEBOOKS_DIR / (
    "1020_calculate-inverse-emissions_single-concentration-projection"
    f"_{INVERSE_EMISSIONS_GHG}_inverse-emissions.ipynb"
)
"""The inverse emissions notebook, as executed in the original run"""

MONTHLY_MEAN_FILE = (
    Path("data")
    / "interim"
    / "monthly-means"
    / f"single-concentration-projection_{INVERSE_EMISSIONS_GHG}_monthly-mean.nc"
)
"""The inverse emissions notebook's input

Relative to the bundle's root directory,
because that is the notebook's working directory.
"""

INVERSE_EMISSIONS_FILE = (
    Path("manuscript-outputs") / f"{INVERSE_EMISSIONS_GHG}_inverse-emissions.feather"
)
"""Where the re-run inverse emissions notebook writes its usual output

Relative to the bundle's root directory.
Not where the original run wrote it, so we never overwrite the original run's data.
"""

INVERSE_EMISSIONS_INTERMEDIATE_FILE = (
    Path("manuscript-outputs")
    / f"{INVERSE_EMISSIONS_GHG}_inverse-emissions-intermediate.csv"
)
"""Where the re-run inverse emissions notebook saves the intermediates we want

Relative to the bundle's root directory.
"""

ONE_BOX_RERUN_UNIT = "ppt"
"""Unit in which the one-box model's re-run of the concentrations is saved"""

INVERSE_EMISSIONS_SAVE_CELL = f"""
# Added for the CMIP7 ScenarioMIP GHG manuscript.
# The original run never saved the one-box model's re-run of the concentrations
# (it is only plotted), but we need it for the manuscript.
manuscript_out_file = Path("{INVERSE_EMISSIONS_INTERMEDIATE_FILE.as_posix()}")
manuscript_out_file.parent.mkdir(exist_ok=True, parents=True)
pd.DataFrame(
    {{"year": years, "run_res": run_res.to("{ONE_BOX_RERUN_UNIT}").m}}
).to_csv(manuscript_out_file, index=False)
manuscript_out_file
"""
"""Cell appended to the inverse emissions notebook to save the intermediates we want"""


def ensure_bundle(bundle_dir: Path = BUNDLE_DIR) -> Path:
    """
    Ensure the parts of the bundle we need are available locally

    Parameters
    ----------
    bundle_dir
        Directory in which to put the bundle

    Returns
    -------
    :
        `bundle_dir`, for convenience
    """
    return ensure_bundle_available(
        files_to_get=BUNDLE_FILES,
        files_to_get_tarred=BUNDLE_TARBALLS,
        bundle_dir=bundle_dir,
        record_id=ZENODO_RECORD_ID,
        tar_strip_components=TARBALL_STRIP_COMPONENTS,
    )


def rerun_inverse_emissions_notebook(
    executed_notebook: Path = INVERSE_EMISSIONS_NOTEBOOK,
    modified_notebook_py: Path | None = None,
    bundle_dir: Path = BUNDLE_DIR,
) -> Path:
    """
    Re-run the inverse emissions notebook, saving the intermediates we want

    The bundle must already be available (see [ensure_bundle][]).
    Its environment is installed if needed,
    which can be slow the first time (see
    [ensure_bundle_environment][local.cmip_ghg_generation.ensure_bundle_environment]).

    Parameters
    ----------
    executed_notebook
        The notebook, as executed in the original run

    modified_notebook_py
        Where to write the modified notebook in `py:percent` format

        If `None`, it goes in [MODIFIED_NOTEBOOKS_DIR][].

    bundle_dir
        Directory which holds the bundle

    Returns
    -------
    :
        Path to the executed, modified notebook
    """
    if modified_notebook_py is None:
        modified_notebook_py = MODIFIED_NOTEBOOKS_DIR / f"{executed_notebook.stem}.py"

    ensure_bundle_environment(bundle_dir)

    ipynb_to_run = bundle_dir / RERUN_NOTEBOOKS_DIR / f"{executed_notebook.stem}.ipynb"
    to_run = write_modified_notebook(
        start_from=executed_notebook,
        out_py=modified_notebook_py,
        out_ipynb=ipynb_to_run,
        extra_cells=[INVERSE_EMISSIONS_SAVE_CELL],
        injected_parameters=[
            f'ghg = "{INVERSE_EMISSIONS_GHG}"',
            f'monthly_mean_file = "{MONTHLY_MEAN_FILE.as_posix()}"',
            f'out_file = "{INVERSE_EMISSIONS_FILE.as_posix()}"',
        ],
    )

    return run_notebook_from_bundle_dir(to_run, ipynb_to_run, bundle_dir=bundle_dir)


def load_inverse_emissions(inverse_emissions_file: Path, unit: str) -> pd.Series:
    """
    Load the inverse emissions written by the inverse emissions notebook

    The notebook writes the emissions in two units,
    a mixing ratio per year (e.g. ppt / yr) and a mass per year (e.g. t CCl4 / yr).

    Parameters
    ----------
    inverse_emissions_file
        File written by the notebook

    unit
        Unit to return the emissions in

        Must be compatible with exactly one of the units the notebook writes.

    Returns
    -------
    :
        The inverse emissions in `unit`, indexed by year

    Raises
    ------
    AssertionError
        Not exactly one of the units the notebook writes is compatible with `unit`
    """
    inverse_emissions = pd.read_feather(inverse_emissions_file)

    res = []
    for row_index, values in inverse_emissions.iterrows():
        row_unit = dict(zip(inverse_emissions.index.names, row_index))["unit"]
        try:
            converted = openscm_units.unit_registry.Quantity(
                values.to_numpy(), row_unit
            ).to(unit)
        except pint.DimensionalityError:
            continue

        res.append(pd.Series(converted.m, index=values.index, name=unit))

    if len(res) != 1:
        msg = (
            f"Expected exactly one row of {inverse_emissions_file} "
            f"to be convertible to {unit!r}, found {len(res)}. "
            f"Units in the file: {inverse_emissions.index.get_level_values('unit')}"
        )
        raise AssertionError(msg)

    res[0].index.name = "year"
    return res[0]


def load_one_box_rerun(intermediate_file: Path) -> pd.Series:
    """
    Load the one-box model's re-run of the concentrations

    Parameters
    ----------
    intermediate_file
        File written by [INVERSE_EMISSIONS_SAVE_CELL][]

    Returns
    -------
    :
        The re-run concentrations in [ONE_BOX_RERUN_UNIT][], indexed by year
    """
    return pd.read_csv(intermediate_file).set_index("year")["run_res"]
