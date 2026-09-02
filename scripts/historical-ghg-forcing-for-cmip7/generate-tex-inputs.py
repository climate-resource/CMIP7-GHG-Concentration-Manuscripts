"""
Compile the tex inputs for the historical GHG manuscript
"""

from pathlib import Path
from typing import Annotated

import typer

from local.cmip_ghg_generation import (
    DEFAULT_BUNDLE_DIR,
    DEFAULT_ORIGINAL_RUN_NOTEBOOKS_DIR,
)
from local.historical_ghg_forcing_for_cmip7 import generate_n2o_methods_figure


def main(
    n2o_methods_figure_file: Annotated[
        Path,
        typer.Option(
            help="Path to in which to write the N2O methods figure. ",
            dir_okay=False,
            file_okay=True,
        ),
    ],
    bundle_dir: Annotated[
        Path,
        typer.Option(
            help=("Directory in which to keep the original run's bundle."),
            dir_okay=True,
            file_okay=False,
        ),
    ] = DEFAULT_BUNDLE_DIR,
    original_run_notebooks_dir: Annotated[
        Path,
        typer.Option(
            help=(
                "Path to the original run's executed notebooks. "
                "Only needed if this repository doesn't already have a copy "
                "of the notebooks we re-run."
            ),
            dir_okay=True,
            file_okay=False,
        ),
    ] = DEFAULT_ORIGINAL_RUN_NOTEBOOKS_DIR,
    force_rerun: Annotated[
        bool,
        typer.Option(
            help=(
                "Re-run the original run's notebooks, "
                "even if we already have the data they produce."
            ),
        ),
    ] = False,
) -> None:
    """
    Create the inputs
    """
    generate_n2o_methods_figure(
        n2o_methods_figure_file,
        bundle_dir=bundle_dir,
        original_run_notebooks_dir=original_run_notebooks_dir,
        force_rerun=force_rerun,
    )


if __name__ == "__main__":
    typer.run(main)
