"""
Generate the scenario results figures, for previewing

The manuscript build makes them too
(`task_results_figure` in `tasks/scenario_ghg_forcing_for_cmip7/task_tex_inputs.py`),
this is just a quick way to make them without running the build.
"""

from typing import Annotated

import typer

from local.scenario_ghg_forcing_for_cmip7.appendix_figures import (
    RESULTS_FIGURE_GASES,
)
from local.scenario_ghg_forcing_for_cmip7.config import get_results_figure_file
from local.scenario_ghg_forcing_for_cmip7.results_figure import (
    generate_results_figure,
)


def main(
    gases: Annotated[
        list[str] | None,
        typer.Argument(help="Gases to generate figures for"),
    ] = None,
) -> None:
    """
    Generate the scenario results figures
    """
    for gas in gases or RESULTS_FIGURE_GASES:
        generate_results_figure(gas, get_results_figure_file(gas))


if __name__ == "__main__":
    typer.run(main)
