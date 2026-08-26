"""
Compile the tex inputs for the historical GHG manuscript
"""

from pathlib import Path
from typing import Annotated

import typer

from local.historical_ghg_forcing_for_cmip7 import generate_n2o_methods_figure

REPO_ROOT = Path(__file__).parents[1]


def main(
    n2o_methods_figure_file: Annotated[
        Path,
        typer.Option(
            help="Path to in which to write the N2O methods figure. ",
            dir_okay=False,
            file_okay=True,
        ),
    ],
) -> None:
    """
    Create the inputs
    """
    generate_n2o_methods_figure(n2o_methods_figure_file)


if __name__ == "__main__":
    typer.run(main)
