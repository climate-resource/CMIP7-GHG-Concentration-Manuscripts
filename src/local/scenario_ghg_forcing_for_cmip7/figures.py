"""
Generation of the scenario GHG manuscript's figures
"""

from __future__ import annotations

from pathlib import Path

import openscm_units
import xarray as xr
from loguru import logger

from local.historical_ghg_forcing_for_cmip7.layout import (
    Panel,
    Row,
    create_figure,
    manuscript_style,
    save_figure,
)
from local.scenario_ghg_forcing_for_cmip7.bundle import (
    ONE_BOX_RERUN_UNIT,
    load_inverse_emissions,
    load_one_box_rerun,
)

INVERSE_EMISSIONS_EMISSIONS_UNIT = "kt CCl4 / yr"
"""Unit in which to plot the inverse emissions"""

INVERSE_EMISSIONS_TITLES = {
    "concentrations": "Concentrations",
    "emissions": "Inverse emissions",
}
"""Title of each panel of the inverse emissions figure"""

INVERSE_EMISSIONS_ROWS = (
    Row(
        panels=(Panel("concentrations"), Panel("emissions")),
        height=1.6,
    ),
)
"""Layout of the inverse emissions figure's panels, top to bottom"""


@manuscript_style
def generate_inverse_emissions_figure(
    outfile: Path,
    inverse_emissions_file: Path,
    intermediate_file: Path,
    monthly_mean_file: Path,
) -> Path:
    """
    Generate the inverse emissions figure

    Parameters
    ----------
    outfile
        File in which to write the figure

    inverse_emissions_file
        Inverse emissions written by the re-run inverse emissions notebook

    intermediate_file
        Intermediates saved by the re-run inverse emissions notebook

    monthly_mean_file
        Monthly-mean concentrations which the notebook inverted

    Returns
    -------
    :
        `outfile`
    """
    monthly_mean = xr.load_dataarray(
        monthly_mean_file, decode_times=xr.coders.CFDatetimeCoder(use_cftime=True)
    ).squeeze("scenario", drop=True)
    one_box_rerun = load_one_box_rerun(intermediate_file)
    inverse_emissions = load_inverse_emissions(
        inverse_emissions_file, unit=INVERSE_EMISSIONS_EMISSIONS_UNIT
    )

    fig, axes = create_figure(INVERSE_EMISSIONS_ROWS)

    ax = axes["concentrations"]
    decimal_year = (
        monthly_mean["time"].dt.year + (monthly_mean["time"].dt.month - 0.5) / 12.0
    )
    monthly_mean_values = openscm_units.unit_registry.Quantity(
        monthly_mean.values, monthly_mean.attrs["units"]
    ).to(ONE_BOX_RERUN_UNIT)
    ax.plot(decimal_year, monthly_mean_values.m, label="Monthly-mean", linewidth=0.8)
    ax.plot(
        one_box_rerun.index + 0.5,
        one_box_rerun.to_numpy(),
        label="One-box model re-run",
        linestyle="--",
        linewidth=0.8,
    )
    ax.set_ylabel(ONE_BOX_RERUN_UNIT)
    ax.legend()

    ax = axes["emissions"]
    ax.plot(inverse_emissions.index + 0.5, inverse_emissions.to_numpy(), linewidth=0.8)
    ax.set_ylabel(INVERSE_EMISSIONS_EMISSIONS_UNIT)

    for ax in axes.values():
        ax.set_xlabel("Year")

    logger.info(f"Writing {outfile}")
    return save_figure(
        fig, axes, INVERSE_EMISSIONS_ROWS, INVERSE_EMISSIONS_TITLES, outfile
    )
