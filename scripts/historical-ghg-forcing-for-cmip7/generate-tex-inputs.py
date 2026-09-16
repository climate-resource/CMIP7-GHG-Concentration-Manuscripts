"""
Compile the tex inputs for the historical GHG manuscript
"""

from pathlib import Path
from typing import Annotated, Optional

import typer

from local.cmip_ghg_generation import (
    DEFAULT_BUNDLE_DIR,
    DEFAULT_ORIGINAL_RUN_NOTEBOOKS_DIR,
)
from local.historical_ghg_forcing_for_cmip7 import (
    SF6_LIKE_GASES,
    generate_ch4_methods_figure,
    generate_co2_methods_figure,
    generate_n2o_methods_figure,
    generate_sf6_like_methods_figure,
)


def parse_sf6_like_methods_figure_file(value: str) -> tuple[str, Path]:
    """
    Parse a gas and a file to write that gas' methods figure in

    Parameters
    ----------
    value
        Value to parse, as `"<gas>=<file>"`

    Returns
    -------
    :
        The gas and the file to write its figure in

    Raises
    ------
    typer.BadParameter
        `value` is not a gas and a file,
        or the gas is not one which is processed like SF6
    """
    gas, _, figure_file = value.partition("=")
    if not figure_file:
        msg = f"Expected '<gas>=<file>', got {value!r}"
        raise typer.BadParameter(msg)

    if gas not in SF6_LIKE_GASES:
        msg = (
            f"{gas!r} is not processed like SF6. "
            f"Expected one of: {', '.join(SF6_LIKE_GASES)}"
        )
        raise typer.BadParameter(msg)

    return gas, Path(figure_file)


def main(  # noqa: PLR0913
    co2_methods_figure_file: Annotated[
        Path,
        typer.Option(
            help="Path to in which to write the CO2 methods figure. ",
            dir_okay=False,
            file_okay=True,
        ),
    ],
    ch4_methods_figure_file: Annotated[
        Path,
        typer.Option(
            help="Path to in which to write the CH4 methods figure. ",
            dir_okay=False,
            file_okay=True,
        ),
    ],
    n2o_methods_figure_file: Annotated[
        Path,
        typer.Option(
            help="Path to in which to write the N2O methods figure. ",
            dir_okay=False,
            file_okay=True,
        ),
    ],
    sf6_like_methods_figure_file: Annotated[
        Optional[list[str]],
        typer.Option(
            help=(
                "Gas which is processed like SF6 "
                "and the path in which to write its methods figure, "
                "as '<gas>=<file>' e.g. 'sf6=figures/sf6_methods.pdf'. "
                "Repeat the option for each gas you want a figure for."
            ),
        ),
    ] = None,
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
    # Parsed before anything is drawn, so a typo in a gas name
    # is caught now rather than three figures from now.
    sf6_like_figures = [
        parse_sf6_like_methods_figure_file(value)
        for value in sf6_like_methods_figure_file or []
    ]

    generate_co2_methods_figure(
        co2_methods_figure_file,
        bundle_dir=bundle_dir,
        original_run_notebooks_dir=original_run_notebooks_dir,
        force_rerun=force_rerun,
    )

    generate_ch4_methods_figure(
        ch4_methods_figure_file,
        bundle_dir=bundle_dir,
        original_run_notebooks_dir=original_run_notebooks_dir,
        force_rerun=force_rerun,
    )

    generate_n2o_methods_figure(
        n2o_methods_figure_file,
        bundle_dir=bundle_dir,
        original_run_notebooks_dir=original_run_notebooks_dir,
        force_rerun=force_rerun,
    )

    # One figure per gas asked for: the gases processed like SF6
    # all share a figure, but there are thirty-four of them,
    # so which ones we draw is the caller's call.
    for gas, figure_file in sf6_like_figures:
        generate_sf6_like_methods_figure(
            gas,
            figure_file,
            bundle_dir=bundle_dir,
            original_run_notebooks_dir=original_run_notebooks_dir,
            force_rerun=force_rerun,
        )


if __name__ == "__main__":
    typer.run(main)
