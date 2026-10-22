"""
Tasks which generate the figures, tables and other tex inputs
"""

from pathlib import Path
from typing import Annotated

from pytask import Product

import local.historical_ghg_forcing_for_cmip7.author_contributions
import local.scenario_ghg_forcing_for_cmip7.figures
import local.scenario_ghg_forcing_for_cmip7.tables
from local.historical_ghg_forcing_for_cmip7.author_contributions import (
    generate_author_contributions,
)
from local.scenario_ghg_forcing_for_cmip7.config import (
    AUTHOR_CONTRIBUTIONS,
    BUNDLE_INVERSE_EMISSIONS,
    BUNDLE_INVERSE_EMISSIONS_INTERMEDIATE,
    BUNDLE_MONTHLY_MEAN,
    INVERSE_EMISSIONS_FIGURE,
    METADATA,
    SUMMARY_TABLE,
)
from local.scenario_ghg_forcing_for_cmip7.figures import (
    generate_inverse_emissions_figure,
)
from local.scenario_ghg_forcing_for_cmip7.tables import generate_summary_table


def task_inverse_emissions_figure(
    inverse_emissions: Path = BUNDLE_INVERSE_EMISSIONS,
    intermediate: Path = BUNDLE_INVERSE_EMISSIONS_INTERMEDIATE,
    monthly_mean: Path = BUNDLE_MONTHLY_MEAN,
    figures_source: Path = Path(local.scenario_ghg_forcing_for_cmip7.figures.__file__),
    outfile: Annotated[Path, Product] = INVERSE_EMISSIONS_FIGURE,
) -> None:
    """Generate the CCl4 inverse emissions figure"""
    generate_inverse_emissions_figure(
        outfile=outfile,
        inverse_emissions_file=inverse_emissions,
        intermediate_file=intermediate,
        monthly_mean_file=monthly_mean,
    )


def task_summary_table(
    inverse_emissions: Path = BUNDLE_INVERSE_EMISSIONS,
    intermediate: Path = BUNDLE_INVERSE_EMISSIONS_INTERMEDIATE,
    tables_source: Path = Path(local.scenario_ghg_forcing_for_cmip7.tables.__file__),
    outfile: Annotated[Path, Product] = SUMMARY_TABLE,
) -> None:
    """Generate the CCl4 summary table"""
    generate_summary_table(
        outfile=outfile,
        inverse_emissions_file=inverse_emissions,
        intermediate_file=intermediate,
    )


def task_author_contributions(
    metadata: Path = METADATA,
    author_contributions_source: Path = Path(
        local.historical_ghg_forcing_for_cmip7.author_contributions.__file__
    ),
    outfile: Annotated[Path, Product] = AUTHOR_CONTRIBUTIONS,
) -> None:
    """Generate the author contributions statement from the metadata file"""
    generate_author_contributions(outfile=outfile, metadata_file=metadata)
