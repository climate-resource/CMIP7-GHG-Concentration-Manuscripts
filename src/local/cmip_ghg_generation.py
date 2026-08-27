"""
Re-running notebooks from the original CMIP7 GHG concentration generation run

The original run is archived at [https://zenodo.org/records/14892947]().
That archive contains everything needed to reproduce the run,
but it does not contain every intermediate object the notebooks built.
When we need one of those intermediate objects for a figure,
the cheapest reliable way to get it
is to re-run the relevant notebook with an extra save call bolted on,
rather than to re-implement the logic that built it.

The one thing the archive does *not* contain is the notebooks *as executed*.
Those live only on the machine which did the original run,
so we keep a copy of the ones we need in this repository
(see [`ensure_executed_notebook_available`][local.cmip_ghg_generation.ensure_executed_notebook_available]).
"""  # noqa: E501

from __future__ import annotations

import os
import shutil
import subprocess
from collections.abc import Iterable
from pathlib import Path

import jupytext
import nbformat
from loguru import logger
from nbformat.notebooknode import NotebookNode

from local.paths import DATA_RAW_DIR
from local.zenodo import download_files, extract_tar_gz

ZENODO_RECORD_ID = "14892947"
"""ID of the Zenodo record which holds the original run"""

TARBALL_STRIP_COMPONENTS = 2
"""Leading path components to strip when extracting the tarballs

Members are prefixed with `output-bundles/v1.0.0/`, which we don't want.
"""

BUNDLE_CONFIG_FILE = "v1.0.0-config-raw.yaml"
"""Config file to use when re-running notebooks

We use the raw config rather than `v1.0.0-config.yaml`
because the latter has the original run's absolute paths baked into it.
The raw config's paths are relative to the bundle's root directory,
which is why notebooks must be re-run with that as their working directory.
"""

DEFAULT_BUNDLE_DIR = DATA_RAW_DIR / "cmip-ghg-concentration-generation" / "v1.0.0"
"""Where we put our copy of the original run's bundle

This is deliberately not tracked by git: it is ~50MB and it is all re-downloadable.
"""

EXECUTED_NOTEBOOKS_DIR = (
    DATA_RAW_DIR / "historical-ghg-forcing-for-cmip7" / "notebooks-executed"
)
"""Where we keep our copy of the original run's notebooks, as executed

This *is* tracked by git, because these notebooks aren't on Zenodo.
"""

MODIFIED_NOTEBOOKS_DIR = (
    DATA_RAW_DIR / "historical-ghg-forcing-for-cmip7" / "notebooks-modified"
)
"""Where we keep the modified notebooks, in jupytext's `py:percent` format

These are generated, so they can go stale.
We track them anyway because being able to see how they change over time
is worth more than the risk of them drifting.
"""

DEFAULT_ORIGINAL_RUN_NOTEBOOKS_DIR = Path(
    "/Users/znicholls/Documents/repos/CMIP-GHG-Concentration-Generation"
    "/output-bundles/v1.0.0/notebooks-executed"
)
"""Where the original run's executed notebooks live on the machine that ran it

Only consulted if the notebook we want isn't already in
[`EXECUTED_NOTEBOOKS_DIR`][local.cmip_ghg_generation.EXECUTED_NOTEBOOKS_DIR].
"""

ENV_VARS_TO_DROP_WHEN_RUNNING_NOTEBOOKS = (
    "VIRTUAL_ENV",
    "PYTHONPATH",
    "PYTHONHOME",
    "CONDA_PREFIX",
    "UV_PROJECT_ENVIRONMENT",
)
"""Environment variables to drop before handing over to pixi

This matters more than it looks.
This repository and the original run both ship a package called `local`,
and this code generally runs inside this repository's virtual environment.
If that environment leaked into the pixi subprocess,
the notebook could silently import the wrong `local`.
"""


def ensure_executed_notebook_available(
    notebook: str | Path,
    original_run_notebooks_dir: Path = DEFAULT_ORIGINAL_RUN_NOTEBOOKS_DIR,
) -> Path:
    """
    Ensure that a notebook from the original run, as executed, is in this repository

    Parameters
    ----------
    notebook
        Notebook to make available,
        as a path relative to the original run's `notebooks-executed` directory
        e.g. `"calculate_n2o_monthly_fifteen_degree_pieces/only/1000_n2o_bin.ipynb"`

    original_run_notebooks_dir
        The original run's `notebooks-executed` directory

        Only used if we don't already have a copy of `notebook`.

    Returns
    -------
        Path to our copy of the notebook

    Raises
    ------
    FileNotFoundError
        We don't have a copy of the notebook
        and it isn't in `original_run_notebooks_dir` either.
    """
    our_copy = EXECUTED_NOTEBOOKS_DIR / notebook
    if our_copy.exists():
        return our_copy

    original = original_run_notebooks_dir / notebook
    if not original.exists():
        msg = (
            f"Could not find {notebook} as executed in the original run. "
            f"It is not in this repository ({our_copy}) "
            f"and it is not in {original_run_notebooks_dir} either. "
            "The executed notebooks are not archived on Zenodo, "
            "so this copy has to be committed to this repository "
            "by someone who has the original run on their machine "
            "(see local.cmip_ghg_generation.DEFAULT_ORIGINAL_RUN_NOTEBOOKS_DIR)."
        )
        raise FileNotFoundError(msg)

    logger.info(f"Copying {original} to {our_copy}")
    our_copy.parent.mkdir(exist_ok=True, parents=True)
    shutil.copy2(original, our_copy)

    return our_copy


def ensure_bundle_available(
    files_to_get: tuple[str, ...],
    files_to_get_tarred: tuple[str, ...],
    bundle_dir: Path = DEFAULT_BUNDLE_DIR,
    record_id: str = ZENODO_RECORD_ID,
    tar_strip_components: int = TARBALL_STRIP_COMPONENTS,
) -> Path:
    """
    Ensure that the original run's bundle is available locally

    Downloads whatever is missing from Zenodo and extracts the tarballs.
    Already downloaded files are not downloaded again.

    Parameters
    ----------
    files_to_get
        Files to get from the bundle

    files_to_get_tarred
        Tar files to get and extract from the bundle

    bundle_dir
        Directory in which to put the bundle

    record_id
        ID of the record from which to get the bundle

    tar_strip_components
        Path components to strip when extracting tarballs

    Returns
    -------
        `bundle_dir`, for convenience
    """
    download_files(
        record_id=record_id,
        filenames=files_to_get,
        dest_dir=bundle_dir,
    )

    tarball_dir = bundle_dir / "zenodo-tarballs"
    tarballs = download_files(
        record_id=record_id,
        filenames=files_to_get_tarred,
        dest_dir=tarball_dir,
    )
    for tarball in tarballs:
        marker = tarball.with_suffix(tarball.suffix + ".extracted")
        if marker.exists():
            logger.debug(f"{tarball} has already been extracted")
            continue

        extract_tar_gz(
            tarball,
            bundle_dir,
            strip_components=tar_strip_components,
        )
        marker.touch()

    return bundle_dir


def ensure_bundle_environment(bundle_dir: Path = DEFAULT_BUNDLE_DIR) -> None:
    """
    Ensure that the original run's environment is installed

    We use pixi with the bundle's own lock file,
    which reproduces the original run's environment exactly.

    Be warned that the first notebook run after this installs the environment
    can take the better part of an hour, almost all of it in the first import.
    That is the operating system checking a few thousand freshly downloaded
    binaries, not anything the notebook is doing; later runs take seconds.

    Parameters
    ----------
    bundle_dir
        Directory which holds the bundle
    """
    logger.info(f"Ensuring the pixi environment in {bundle_dir} is installed")
    # `--frozen` so pixi uses the archived lock file as-is,
    # rather than trying to re-solve it against today's package index.
    run_in_bundle(
        ["pixi", "install", "--frozen"],
        bundle_dir=bundle_dir,
    )


def run_in_bundle(
    command: list[str],
    bundle_dir: Path,
    env_vars_to_drop: tuple[str, ...] = ENV_VARS_TO_DROP_WHEN_RUNNING_NOTEBOOKS,
) -> None:
    """
    Run a command with the bundle as its working directory

    The environment this repository runs in is stripped out first,
    see [`ENV_VARS_TO_DROP_WHEN_RUNNING_NOTEBOOKS`][local.cmip_ghg_generation.ENV_VARS_TO_DROP_WHEN_RUNNING_NOTEBOOKS].

    Parameters
    ----------
    command
        Command to run

    bundle_dir
        Directory which holds the bundle

    env_vars_to_drop
        Environment variables to drop when running the command
    """  # noqa: E501
    env = {
        k: v
        for k, v in os.environ.items()
        # Make this a function argument rather than a global
        if k not in ENV_VARS_TO_DROP_WHEN_RUNNING_NOTEBOOKS
    }

    logger.debug(f"Running {command} in {bundle_dir}")
    # The command is built by this module, it is not user input, hence noqa S603
    subprocess.run(command, cwd=bundle_dir, env=env, check=True)  # noqa: S603


def write_modified_notebook(  # noqa: PLR0913
    start_from: Path,
    out_py: Path,
    out_ipynb: Path,
    extra_cells: Iterable[str],
    config_file: str = BUNDLE_CONFIG_FILE,
    step_config_id: str | None = None,
) -> Path:
    """
    Write a modified copy of a notebook from the original run

    The modifications are:

    1. the injected parameters are rewritten
       so the notebook uses our copy of the bundle rather than the original run's
    1. all outputs are cleared
    1. `extra_cells` are appended

    Two copies are written: a `py:percent` copy, which is what we track in git
    because it is the copy a human can read a diff of,
    and an `ipynb` copy, which is what we actually intend to be executed.

    Parameters
    ----------
    start_from
        Notebook to modify

    out_py
        Path in which to write the modified notebook in `py:percent` format

    out_ipynb
        Path in which to write the modified notebook in `ipynb` format

    extra_cells
        Source of each code cell to append to the notebook

    config_file
        Config file for the notebook to use

        This is interpreted relative to the notebook's working directory,
        which is the bundle's root directory.

    step_config_id
        Step config ID for the notebook to use

        If `None`, whatever the original run used is kept.

    Returns
    -------
        `out_ipynb`, i.e. the notebook to execute
    """
    notebook = jupytext.read(start_from)

    injected_parameters = [f'config_file = "{config_file}"']
    if step_config_id is not None:
        injected_parameters.append(f'step_config_id = "{step_config_id}"')

    _set_injected_parameters(notebook, injected_parameters)
    _clear_run_artefacts(notebook)

    for extra_cell in extra_cells:
        notebook.cells.append(nbformat.v4.new_code_cell(extra_cell.strip()))

    for out_file, fmt in ((out_py, "py:percent"), (out_ipynb, "ipynb")):
        out_file.parent.mkdir(exist_ok=True, parents=True)
        logger.info(f"Writing {out_file}")
        jupytext.write(notebook, out_file, fmt=fmt)

    return out_ipynb


def _set_injected_parameters(
    notebook: NotebookNode, parameter_lines: list[str]
) -> None:
    """
    Overwrite a notebook's injected parameters cell

    Parameters
    ----------
    notebook
        Notebook to modify, in place

    parameter_lines
        Lines to use as the injected parameters

    Raises
    ------
    ValueError
        The notebook has no cell tagged `injected-parameters`.
    """
    for cell in notebook.cells:
        if "injected-parameters" in cell.get("metadata", {}).get("tags", []):
            cell["source"] = "\n".join(["# Parameters", *parameter_lines])
            return

    msg = (
        "Could not find a cell tagged 'injected-parameters'. "
        "This function expects a notebook as executed by papermill."
    )
    raise ValueError(msg)


def _clear_run_artefacts(notebook: NotebookNode) -> None:
    """
    Clear the traces of the original run from a notebook

    Parameters
    ----------
    notebook
        Notebook to modify, in place
    """
    for cell in notebook.cells:
        # Timings from the original run mean nothing for our re-run,
        # and they make the tracked `py:percent` copy unreadable.
        cell.get("metadata", {}).pop("papermill", None)
        cell.get("metadata", {}).pop("execution", None)

        if cell["cell_type"] != "code":
            continue

        cell["outputs"] = []
        cell["execution_count"] = None


def run_notebook_from_bundle_dir(
    notebook: Path,
    out_notebook: Path,
    bundle_dir: Path = DEFAULT_BUNDLE_DIR,
) -> Path:
    """
    Run a notebook in the original run's environment

    Parameters
    ----------
    notebook
        Notebook to run

    out_notebook
        Path in which to write the executed notebook

    bundle_dir
        Directory which holds the bundle

    Returns
    -------
        `out_notebook`
    """
    out_notebook.parent.mkdir(exist_ok=True, parents=True)

    logger.info(f"Running {notebook}")
    run_in_bundle(
        [
            "pixi",
            "run",
            "--frozen",
            "papermill",
            # Show the progress by default.
            # Can add a flag later if we need to be able to toggle this.
            # "--no-progress-bar",
            # The config file's paths are relative to the bundle's root directory
            "--cwd",
            str(bundle_dir.resolve()),
            str(notebook.resolve()),
            str(out_notebook.resolve()),
        ],
        bundle_dir=bundle_dir,
    )

    return out_notebook
