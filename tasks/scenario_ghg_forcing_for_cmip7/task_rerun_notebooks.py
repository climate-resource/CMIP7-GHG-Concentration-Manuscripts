"""
Expensive tasks: getting the original run's bundle and re-running its notebooks

These are kept in their own module on purpose.
pytask considers a task changed if the module it is defined in changes,
so if these sat with the cheap tasks,
editing a cheap task would re-run the notebooks too.
"""

from pathlib import Path
from typing import Annotated

from pytask import Product, task

import local.scenario_ghg_forcing_for_cmip7.bundle
import local.scenario_ghg_forcing_for_cmip7.magicc_output
import local.scenario_ghg_forcing_for_cmip7.rerun_projections
from local.cmip_ghg_generation import DEFAULT_BUNDLE_DIR as HISTORICAL_BUNDLE_DIR
from local.cmip_ghg_generation import ensure_bundle_available
from local.scenario_ghg_forcing_for_cmip7.bundle import (
    BUNDLE_DIR,
    ensure_bundle,
    rerun_inverse_emissions_notebook,
)
from local.scenario_ghg_forcing_for_cmip7.config import (
    BUNDLE_EXTRACTED_MARKERS,
    BUNDLE_INVERSE_EMISSIONS,
    BUNDLE_INVERSE_EMISSIONS_INTERMEDIATE,
    BUNDLE_MONTHLY_MEAN,
    BUNDLE_PIXI_LOCK,
    BUNDLE_PYPROJECT,
    EXECUTED_INVERSE_EMISSIONS_NOTEBOOK,
    HISTORICAL_BUNDLE_EXTRACTED_MARKERS,
    HISTORICAL_BUNDLE_TARBALLS,
    MAGICC_SUMMARY,
    MODIFIED_INVERSE_EMISSIONS_NOTEBOOK,
    RERUN_CO2_SEASONALITY_FILES,
    RERUN_EMISSIONS_FILES,
    RERUN_LAT_GRADIENT_GASES,
    get_rerun_lat_gradient_files,
)
from local.scenario_ghg_forcing_for_cmip7.magicc_output import summarise_magicc_output
from local.scenario_ghg_forcing_for_cmip7.rerun_projections import (
    MAGICC_OUTPUT_DB_DIR,
    ONE_EOF_LAT_GRADIENT_NOTEBOOK,
    SEASONALITY_NOTEBOOK,
    TWO_EOF_GASES,
    TWO_EOF_LAT_GRADIENT_NOTEBOOK,
    ensure_rerun_inputs,
    extract_emissions,
    rerun_co2_seasonality,
    rerun_one_eof_lat_gradient,
    rerun_two_eof_lat_gradient,
)

RERUN_SOURCE = Path(local.scenario_ghg_forcing_for_cmip7.rerun_projections.__file__)
"""The module which defines the re-runs (they re-run if it changes)"""


def task_download_bundle(
    pyproject: Annotated[Path, Product] = BUNDLE_PYPROJECT,
    pixi_lock: Annotated[Path, Product] = BUNDLE_PIXI_LOCK,
    monthly_mean: Annotated[Path, Product] = BUNDLE_MONTHLY_MEAN,
    extracted_markers: Annotated[list[Path], Product] = BUNDLE_EXTRACTED_MARKERS,
) -> None:
    """Download the parts of the original run's bundle which we need"""
    ensure_bundle(BUNDLE_DIR)


def task_download_historical_bundle(
    extracted_markers: Annotated[
        list[Path], Product
    ] = HISTORICAL_BUNDLE_EXTRACTED_MARKERS,
) -> None:
    """Download the parts of the historical run's bundle which we need"""
    ensure_bundle_available(
        files_to_get=(),
        files_to_get_tarred=HISTORICAL_BUNDLE_TARBALLS,
        bundle_dir=HISTORICAL_BUNDLE_DIR,
    )


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


def task_summarise_magicc_output(
    extracted_markers: list[Path] = BUNDLE_EXTRACTED_MARKERS,
    magicc_output_source: Path = Path(
        local.scenario_ghg_forcing_for_cmip7.magicc_output.__file__
    ),
    out_file: Annotated[Path, Product] = MAGICC_SUMMARY,
) -> None:
    """Summarise the original run's MAGICC output (GBs, we need a few quantiles)"""
    summarise_magicc_output(BUNDLE_DIR / MAGICC_OUTPUT_DB_DIR, out_file)


def task_extract_rerun_emissions(
    extracted_markers: list[Path] = BUNDLE_EXTRACTED_MARKERS,
    pixi_lock: Path = BUNDLE_PIXI_LOCK,
    rerun_source: Path = RERUN_SOURCE,
    out_files: Annotated[list[Path], Product] = RERUN_EMISSIONS_FILES,
) -> None:
    """Extract the emissions which drive the latitudinal gradients (all markers)"""
    ensure_rerun_inputs(BUNDLE_DIR)
    extract_emissions(BUNDLE_DIR)


def task_rerun_co2_seasonality(
    executed_notebook: Path = SEASONALITY_NOTEBOOK,
    extracted_markers: list[Path] = BUNDLE_EXTRACTED_MARKERS,
    rerun_source: Path = RERUN_SOURCE,
    out_files: Annotated[list[Path], Product] = RERUN_CO2_SEASONALITY_FILES,
) -> None:
    """Re-run the CO2 seasonality projection for all the markers"""
    ensure_rerun_inputs(BUNDLE_DIR)
    rerun_co2_seasonality(BUNDLE_DIR)


for _ghg in RERUN_LAT_GRADIENT_GASES:
    _executed_notebook = (
        TWO_EOF_LAT_GRADIENT_NOTEBOOK
        if _ghg in TWO_EOF_GASES
        else ONE_EOF_LAT_GRADIENT_NOTEBOOK
    )

    @task(id=_ghg)
    def task_rerun_lat_gradient(  # noqa: PLR0913
        ghg: str = _ghg,
        executed_notebook: Path = _executed_notebook,
        emissions_files: list[Path] = RERUN_EMISSIONS_FILES,
        extracted_markers: list[Path] = BUNDLE_EXTRACTED_MARKERS,
        rerun_source: Path = RERUN_SOURCE,
        out_files: Annotated[list[Path], Product] = get_rerun_lat_gradient_files(_ghg),
    ) -> None:
        """Re-run a gas' latitudinal gradient projection for all the markers"""
        ensure_rerun_inputs(BUNDLE_DIR)
        if ghg in TWO_EOF_GASES:
            rerun_two_eof_lat_gradient(ghg, BUNDLE_DIR)
        else:
            rerun_one_eof_lat_gradient(ghg, BUNDLE_DIR)
