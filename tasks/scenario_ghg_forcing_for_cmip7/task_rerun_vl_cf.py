"""
Expensive task: re-running the CO2 seasonality projection for VL-CF

In its own module, for the same reason as `task_rerun_notebooks`.
"""

from pathlib import Path
from typing import Annotated

from pytask import Product

import local.scenario_ghg_forcing_for_cmip7.rerun_vl_cf
from local.scenario_ghg_forcing_for_cmip7.bundle import BUNDLE_DIR
from local.scenario_ghg_forcing_for_cmip7.config import (
    BUNDLE_EXTRACTED_MARKERS,
    RERUN_CO2_SEASONALITY_VL_CF_FILES,
)
from local.scenario_ghg_forcing_for_cmip7.rerun_projections import (
    SEASONALITY_NOTEBOOK,
    ensure_rerun_inputs,
)
from local.scenario_ghg_forcing_for_cmip7.rerun_vl_cf import (
    rerun_co2_seasonality_vl_cf,
)


def task_rerun_co2_seasonality_vl_cf(
    executed_notebook: Path = SEASONALITY_NOTEBOOK,
    extracted_markers: list[Path] = BUNDLE_EXTRACTED_MARKERS,
    rerun_source: Path = Path(
        local.scenario_ghg_forcing_for_cmip7.rerun_vl_cf.__file__
    ),
    out_files: Annotated[list[Path], Product] = RERUN_CO2_SEASONALITY_VL_CF_FILES,
) -> None:
    """Re-run the CO2 seasonality projection for VL-CF"""
    ensure_rerun_inputs(BUNDLE_DIR)
    rerun_co2_seasonality_vl_cf(BUNDLE_DIR)
