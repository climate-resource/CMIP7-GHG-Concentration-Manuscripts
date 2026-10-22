"""
Generate the scenario results figures, for previewing

Not part of the manuscript build (yet).
"""

from typing import Annotated

import typer

from local.scenario_ghg_forcing_for_cmip7.config import FIGURES_OUT_DIR
from local.scenario_ghg_forcing_for_cmip7.results_figure import (
    generate_results_figure,
)

DEFAULT_GASES = ["co2", "ch4", "n2o", "cfc12eq", "hfc134aeq"]
"""Gases to generate figures for if none are given"""


def main(
    gases: Annotated[
        list[str] | None,
        typer.Argument(help="Gases to generate figures for"),
    ] = None,
) -> None:
    """
    Generate the scenario results figures
    """
    for gas in gases or DEFAULT_GASES:
        generate_results_figure(gas, FIGURES_OUT_DIR / f"{gas}_results.pdf")


if __name__ == "__main__":
    typer.run(main)
