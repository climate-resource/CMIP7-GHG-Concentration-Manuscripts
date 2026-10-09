"""
Getting the variance each EOF explains out of the original run

The original run only ever saved the EOFs it kept.
The ones it dropped are exactly the ones which say
that dropping them costs almost nothing,
so they are what a reader needs to see to judge the choice,
and they only ever existed inside the notebook which calculated them.

Every gas is in the same position and the way out is the same for all of them:
re-run that notebook with a cell appended which saves the full decomposition.
The recipe therefore lives here rather than once per gas.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd
from loguru import logger

from local.cmip_ghg_generation import (
    DEFAULT_BUNDLE_DIR,
    DEFAULT_ORIGINAL_RUN_NOTEBOOKS_DIR,
    MODIFIED_NOTEBOOKS_DIR,
    ensure_bundle_available,
    ensure_bundle_environment,
    ensure_executed_notebook_available,
    run_notebook_from_bundle_dir,
    write_modified_notebook,
)

EOF_COLUMN = "eof"
"""Column which holds the EOF in the variance explained table"""

VARIANCE_EXPLAINED_COLUMN = "variance_explained_fraction"
"""Column which holds the fraction of the variance explained"""


@dataclass(frozen=True)
class DecompositionToSave:
    """
    One EOF decomposition to get out of a re-run notebook
    """

    eofs_pcs_variable: str
    """Name of the notebook variable which holds the full EOFs and PCs

    This is the decomposition before the original run selects
    the EOFs it keeps, i.e. the one with every EOF still in it.
    """

    variance_explained_file: Path
    """Where to save the variance each EOF explains

    Relative to the bundle's root directory,
    because that is the notebook's working directory.
    """

    full_eofs_pcs_file: Path
    """Where to save the full decomposition

    We keep the whole thing, not just the variance explained we derive from
    it, so the decomposition can be looked at again without paying for
    another re-run.

    Relative to the bundle's root directory,
    because that is the notebook's working directory.
    """


def build_save_cell(decomposition: DecompositionToSave) -> str:
    """
    Build the cell which saves a decomposition out of a re-run notebook

    Parameters
    ----------
    decomposition
        Decomposition the cell should save

    Returns
    -------
    :
        Source of the cell, to append to the notebook
    """
    return f"""
# Added for the CMIP7 GHG manuscript.
# The original run only ever saved the EOFs it keeps,
# but we want to show how much of the variance every EOF explains,
# so that keeping only the first few can be justified.
from pathlib import Path

import numpy as np
import pandas as pd

# The EOFs are the right singular vectors of the residuals, so they are
# orthonormal, which makes the principal components uncorrelated:
# the principal components' cross-product is `diag(D) ** 2`,
# i.e. the singular values squared with nothing off the diagonal.
# Bound to a name of our own, so nothing here can shadow the notebook's.
manuscript_eofs_pcs = {decomposition.eofs_pcs_variable}
pcs = manuscript_eofs_pcs["principal-components"].transpose("year", "eof").data.m
singular_values_squared = pcs.T @ pcs

off_diagonal = singular_values_squared - np.diag(np.diag(singular_values_squared))
if not np.allclose(off_diagonal, 0.0, atol=1e-10 * np.trace(singular_values_squared)):
    msg = "The principal components are not uncorrelated, so this is not an SVD"
    raise AssertionError(msg)

variance_explained = np.diag(singular_values_squared) / np.trace(
    singular_values_squared
)

full_eofs_pcs_out_file = Path("{decomposition.full_eofs_pcs_file.as_posix()}")
full_eofs_pcs_out_file.parent.mkdir(exist_ok=True, parents=True)
to_save = manuscript_eofs_pcs.pint.dequantify()
# A decomposition can be over a stacked dimension -- CO2's seasonality change
# is over latitude and month together -- and netCDF cannot hold a MultiIndex.
# Unstacking it into its levels keeps everything, in a form netCDF can write.
multi_indexes = [
    name
    for name, index in to_save.indexes.items()
    if isinstance(index, pd.MultiIndex)
]
if multi_indexes:
    to_save = to_save.reset_index(multi_indexes)

to_save.to_netcdf(full_eofs_pcs_out_file)

manuscript_out_file = Path("{decomposition.variance_explained_file.as_posix()}")
manuscript_out_file.parent.mkdir(exist_ok=True, parents=True)
pd.DataFrame(
    {{
        "{EOF_COLUMN}": manuscript_eofs_pcs["eof"].values,
        "{VARIANCE_EXPLAINED_COLUMN}": variance_explained,
    }}
).to_csv(manuscript_out_file, index=False)
manuscript_out_file
"""


def get_variance_explained(  # noqa: PLR0913
    base_notebook: Path,
    decompositions: tuple[DecompositionToSave, ...],
    step_config_id: str = "only",
    bundle_dir: Path = DEFAULT_BUNDLE_DIR,
    original_run_notebooks_dir: Path = DEFAULT_ORIGINAL_RUN_NOTEBOOKS_DIR,
    force_rerun: bool = False,
) -> tuple[pd.DataFrame, ...]:
    """
    Get the variance explained by each EOF of one or more decompositions

    All of `decompositions` must come out of the same notebook.
    They are fetched together because the notebook is re-run as a whole:
    asking for them one at a time would run it once for each of them.

    Parameters
    ----------
    base_notebook
        Notebook to re-run,
        as a path relative to the original run's `notebooks-executed` directory

    decompositions
        Decompositions to get out of `base_notebook`

    step_config_id
        Step config ID for the notebook to use

        This is the gas itself for the steps which run once per gas.

    bundle_dir
        Directory in which to keep the original run's bundle

    original_run_notebooks_dir
        The original run's `notebooks-executed` directory

        Only used if we don't already have a copy of the notebook we need.

    force_rerun
        Re-run the notebook even if its output is already there

    Returns
    -------
    :
        The variance explained by each EOF, one table per decomposition,
        in the order the decompositions were given
    """
    out_files = [
        bundle_dir / decomposition.variance_explained_file
        for decomposition in decompositions
    ]
    needed = [
        path
        for decomposition in decompositions
        for path in (
            bundle_dir / decomposition.variance_explained_file,
            bundle_dir / decomposition.full_eofs_pcs_file,
        )
    ]
    if all(path.exists() for path in needed) and not force_rerun:
        logger.info(f"Using existing {', '.join(str(path) for path in out_files)}")
        return tuple(pd.read_csv(path) for path in out_files)

    start_from = ensure_executed_notebook_available(
        base_notebook,
        original_run_notebooks_dir=original_run_notebooks_dir,
    )
    ensure_bundle_available(
        files_to_get=(
            "pyproject.toml",
            "pixi.lock",
            "v1.0.0-config-raw.yaml",
        ),
        files_to_get_tarred=(
            "src.tar.gz",
            "data--interim.tar.gz",
        ),
        bundle_dir=bundle_dir,
    )
    ensure_bundle_environment(bundle_dir)

    notebook_name = base_notebook.stem
    ipynb_to_run = bundle_dir / "notebooks-rerun" / f"{notebook_name}.ipynb"
    # The steps which run once per gas share one notebook between all of them,
    # so the tracked copy is of whichever gas was asked for last. What we track
    # it for is the cell we append, which is the same for every gas, so one
    # copy says everything a reader needs.
    to_run = write_modified_notebook(
        start_from=start_from,
        out_py=MODIFIED_NOTEBOOKS_DIR / f"{notebook_name}.py",
        out_ipynb=ipynb_to_run,
        extra_cells=[
            build_save_cell(decomposition) for decomposition in decompositions
        ],
        step_config_id=step_config_id,
    )
    run_notebook_from_bundle_dir(
        to_run,
        ipynb_to_run,
        bundle_dir=bundle_dir,
    )

    return tuple(pd.read_csv(path) for path in out_files)
