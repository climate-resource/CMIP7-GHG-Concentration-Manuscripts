"""
Expensive tasks: getting the original run's bundle and re-running its notebooks

These are kept in their own module on purpose.
pytask considers a task changed if the module it is defined in changes,
so if these sat with the cheap tasks,
editing a cheap task would re-run the notebooks too.
"""

from pathlib import Path
from typing import Annotated

from pytask import Product

import local.scenario_ghg_forcing_for_cmip7.bundle
from local.scenario_ghg_forcing_for_cmip7.bundle import (
    BUNDLE_DIR,
    ensure_bundle,
    rerun_inverse_emissions_notebook,
)
from local.scenario_ghg_forcing_for_cmip7.config import (
    BUNDLE_INVERSE_EMISSIONS,
    BUNDLE_INVERSE_EMISSIONS_INTERMEDIATE,
    BUNDLE_MONTHLY_MEAN,
    BUNDLE_PIXI_LOCK,
    BUNDLE_PYPROJECT,
    EXECUTED_INVERSE_EMISSIONS_NOTEBOOK,
    MODIFIED_INVERSE_EMISSIONS_NOTEBOOK,
)


def task_download_bundle(
    pyproject: Annotated[Path, Product] = BUNDLE_PYPROJECT,
    pixi_lock: Annotated[Path, Product] = BUNDLE_PIXI_LOCK,
    monthly_mean: Annotated[Path, Product] = BUNDLE_MONTHLY_MEAN,
) -> None:
    """Download the parts of the original run's bundle which we need"""
    ensure_bundle(BUNDLE_DIR)


def task_rerun_inverse_emissions_notebook(  # noqa: PLR0913
    executed_notebook: Path = EXECUTED_INVERSE_EMISSIONS_NOTEBOOK,
    bundle_source: Path = Path(local.scenario_ghg_forcing_for_cmip7.bundle.__file__),
    pixi_lock: Path = BUNDLE_PIXI_LOCK,
    monthly_mean: Path = BUNDLE_MONTHLY_MEAN,
    modified_notebook_py: Annotated[
        Path, Product
    ] = MODIFIED_INVERSE_EMISSIONS_NOTEBOOK,
    intermediate: Annotated[Path, Product] = BUNDLE_INVERSE_EMISSIONS_INTERMEDIATE,
    inverse_emissions: Annotated[Path, Product] = BUNDLE_INVERSE_EMISSIONS,
) -> None:
    """
    Re-run the inverse emissions notebook, saving the intermediates we want

    Runs in the bundle's own (pixi) environment.
    """
    rerun_inverse_emissions_notebook(
        executed_notebook=executed_notebook,
        modified_notebook_py=modified_notebook_py,
        bundle_dir=BUNDLE_DIR,
    )
