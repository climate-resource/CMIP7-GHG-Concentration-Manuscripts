"""
Expensive tasks: running MAGICC with and without the CMIP7 history

Kept in their own module, for the same reason as `task_rerun_notebooks.py`:
pytask considers a task changed if the module it is defined in changes,
and these runs take a while.
"""

from pathlib import Path
from typing import Annotated

from pytask import Product, task

import local.scenario_ghg_forcing_for_cmip7.magicc_jump
import local.scenario_ghg_forcing_for_cmip7.magicc_jump_summary
from local.scenario_ghg_forcing_for_cmip7.bundle import BUNDLE_DIR
from local.scenario_ghg_forcing_for_cmip7.config import (
    BUNDLE_EXTRACTED_MARKERS,
    BUNDLE_PIXI_LOCK,
    HISTORICAL_BUNDLE_EXTRACTED_MARKERS,
)
from local.scenario_ghg_forcing_for_cmip7.magicc_jump import (
    EXPERIMENTS,
    MAGICC_JUMP_NOTEBOOK,
    SETUPS,
    get_out_file,
    run_magicc_jump,
)
from local.scenario_ghg_forcing_for_cmip7.magicc_jump_summary import (
    MAGICC_JUMP_SUMMARY,
    summarise_magicc_jump,
)

MAGICC_JUMP_SOURCE = Path(local.scenario_ghg_forcing_for_cmip7.magicc_jump.__file__)
"""The module which defines the runs (they re-run if it changes)"""

MAGICC_JUMP_SUMMARY_SOURCE = Path(
    local.scenario_ghg_forcing_for_cmip7.magicc_jump_summary.__file__
)
"""The module which compares the runs (only the summary re-runs if it changes)"""

MAGICC_JUMP_OUT_FILES = [
    BUNDLE_DIR / get_out_file(name, setup)
    for _, _, name in EXPERIMENTS
    for setup in SETUPS
]
"""Every run's output"""

for _model, _scenario, _name in EXPERIMENTS:
    for _setup in SETUPS:

        @task(id=f"{_name}-{_setup}")
        def task_run_magicc_jump(  # noqa: PLR0913
            model: str = _model,
            scenario: str = _scenario,
            name: str = _name,
            setup: str = _setup,
            notebook: Path = MAGICC_JUMP_NOTEBOOK,
            magicc_jump_source: Path = MAGICC_JUMP_SOURCE,
            pixi_lock: Path = BUNDLE_PIXI_LOCK,
            extracted_markers: list[Path] = BUNDLE_EXTRACTED_MARKERS,
            historical_extracted_markers: list[
                Path
            ] = HISTORICAL_BUNDLE_EXTRACTED_MARKERS,
            out_file: Annotated[Path, Product] = BUNDLE_DIR
            / get_out_file(_name, _setup),
        ) -> None:
            """Run MAGICC for one experiment in one set-up"""
            run_magicc_jump(model, scenario, name, setup)


def task_summarise_magicc_jump(
    out_files: list[Path] = MAGICC_JUMP_OUT_FILES,
    magicc_jump_summary_source: Path = MAGICC_JUMP_SUMMARY_SOURCE,
    summary: Annotated[Path, Product] = BUNDLE_DIR / MAGICC_JUMP_SUMMARY,
) -> None:
    """Summarise the runs (reproduction check and the jump's effect)"""
    summarise_magicc_jump(BUNDLE_DIR)
