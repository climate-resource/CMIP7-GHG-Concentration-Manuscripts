"""
Generation of the C8F18 methods figure

The pieces this figure shares with the CH4 methods figure
live in [local.historical_ghg_forcing_for_cmip7.plotting][].
What is here is the data this figure loads,
the panels it has and how they are laid out.
"""
# Differences from ch4
# - nothing is derived from measurements here.
#   There have been no new measurements of C8F18 since Ivy et al. (2012),
#   so the original run simply took the CMIP6 values,
#   historical up to 2014 and SSP2-4.5 after it.
#   That leaves nothing to say about an observational network,
#   a decomposition or an extension: the components are what CMIP6 had.
# - seasonality is zero everywhere, so it has no panel.
# - the latitudinal gradient EOF is assumed, as it is for the C4F10-like gases,
#   and the PC that goes with it is the difference between CMIP6's
#   northern and southern hemispheric means.

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from loguru import logger

from local.historical_ghg_forcing_for_cmip7.layout import (
    Panel,
    Row,
    create_figure,
    label_panels,
    lay_out_figure,
)
from local.historical_ghg_forcing_for_cmip7.plotting import (
    BROKEN_SPLIT,
    add_colour_bar,
    get_only_data_variable,
    label_name,
    linear_latitudinal_gradient_eof,
    plot_extension_pieces,
    plot_flying_carpet,
    plot_lat_gradient_eofs,
    plot_monthly_means,
    plot_yearly_means,
)

GAS = "c8f18"
"""Gas this figure is drawn for

The original run's `calculate_c8f18_like_monthly_fifteen_degree_pieces` step
has only this one gas in it, and its notebook refuses to run for any other,
so, unlike the other groups, there is nothing here to parameterise.
"""

CMIP6_HISTORICAL_LAST_YEAR = 2014
"""Last year the CMIP6 historical dataset covers

Everything after this comes from CMIP6's SSP2-4.5 projection.
"""

SOURCE_LABELS = {
    "historical": "CMIP6 historical",
    "projection": "CMIP6 SSP2-4.5",
}
"""How each of the two sources is named to a reader"""

EXTENSION_SPLIT_MARGIN = 30
"""Years to leave before the record starts when breaking an extended panel

Both extended components are zero for every year
before C8F18 was manufactured, and everything these panels are about
happens after that. So the axis is broken just before the first year
with anything in it, rather than at some fixed year:
the flat run-up goes in the left half and the rest gets the right half.
This is how much of the flat run-up is kept on the right,
so the reader can see the line was flat before the record starts.
"""

SH_LAT = -45.0
"""Latitude which stands for the southern hemisphere in the hemispheric means"""

TITLES = {
    "gm-ext": "Extended global-mean",
    "lat-grad-pc-ext": "Extended lat. gradient PC",
    "lat-grad-eof": "Lat. gradient EOF",
    "monthly": "Monthly spatial-means",
    "yearly": "Yearly spatial-means",
    "flying-carpet": "Native resolution",
}
"""Title of each panel"""

ROWS = (
    Row(
        panels=(
            Panel("gm-ext", width=1.5, broken=True, broken_split=BROKEN_SPLIT),
            Panel("lat-grad-pc-ext", width=1.5, broken=True, broken_split=BROKEN_SPLIT),
            Panel("lat-grad-eof"),
        ),
        height=2.3,
    ),
    Row(
        panels=(
            Panel("monthly"),
            Panel("yearly", width=1.5, broken=True, broken_split=BROKEN_SPLIT),
            Panel(
                "flying-carpet",
                aspect=1.0,
                projection="3d",
                colour_bar=True,
                colour_bar_height=0.6,
            ),
        ),
        height=2.8,
    ),
)
"""Layout of the figure's panels, top to bottom

Two rows is all this gas needs:
the components it is built from, and the outputs they make.
There is no observational network row and no extension row,
because the components are taken from CMIP6 whole.
"""


def interim_dir(bundle_dir: Path) -> Path:
    """
    Get the directory which holds C8F18's interim data

    Parameters
    ----------
    bundle_dir
        Directory which holds the original run's bundle

    Returns
    -------
        Directory which holds C8F18's interim data
    """
    return bundle_dir / "data" / "interim" / GAS


def get_source_pieces(
    years: np.typing.NDArray[np.int64],
) -> dict[str, np.typing.NDArray[np.int64]]:
    """
    Get the years which came from each of the two CMIP6 datasets

    Parameters
    ----------
    years
        Years the component covers

    Returns
    -------
        Years which came from each source, labelled by source
    """
    return {
        SOURCE_LABELS["historical"]: years[years <= CMIP6_HISTORICAL_LAST_YEAR],
        SOURCE_LABELS["projection"]: years[years > CMIP6_HISTORICAL_LAST_YEAR],
    }


def get_split_year(da: xr.DataArray, margin: int = EXTENSION_SPLIT_MARGIN) -> int:
    """
    Get the year to break an extended panel's axis at

    Parameters
    ----------
    da
        Component the panel shows, with a `year` dimension

    margin
        Years of the flat run-up to keep on the right of the break

    Returns
    -------
        Year to break the axis at
    """
    non_zero = np.nonzero(da.values)[0]

    return int(da["year"].values[non_zero[0]]) - margin


def get_latitudinal_gradient_pc(bundle_dir: Path) -> xr.DataArray:
    """
    Get C8F18's latitudinal gradient principal component

    The original run never wrote this out for C8F18,
    but it does not have to be re-derived either:
    the EOF is normalised so that its hemispheric means are plus and minus
    a half, which makes the PC exactly the difference between
    the northern and southern hemispheric means, which is how the original
    run built it and what the hemispheric mean output holds.

    Parameters
    ----------
    bundle_dir
        Directory which holds the original run's bundle

    Returns
    -------
        C8F18's latitudinal gradient principal component, by year
    """
    hemispheric_mean = get_only_data_variable(
        xr.load_dataset(
            interim_dir(bundle_dir) / f"{GAS}_hemispheric-mean_annual-mean.nc"
        )
    )

    northern = hemispheric_mean.sel(lat=-SH_LAT, drop=True)
    southern = hemispheric_mean.sel(lat=SH_LAT, drop=True)

    res = northern - southern
    res.attrs = hemispheric_mean.attrs

    return res


def generate_c8f18_methods_figure(
    outfile: Path,
    bundle_dir: Path,
    force_rerun: bool = False,
) -> Path:
    """
    Generate the C8F18 methods figure

    Parameters
    ----------
    outfile
        File in which to write the figure

    bundle_dir
        Directory which holds the original run's bundle

    force_rerun
        Re-generate the figure, even if the output file already exists

    Returns
    -------
    :
        `outfile`
    """
    if outfile.exists() and not force_rerun:
        logger.info(f"Using existing {outfile}")
        return outfile

    gas_dir = interim_dir(bundle_dir)

    fig, axes = create_figure(ROWS)

    global_mean = get_only_data_variable(
        xr.load_dataset(gas_dir / f"{GAS}_global-mean_annual-mean.nc")
    )
    split_year = get_split_year(global_mean)

    plot_extension_pieces(
        global_mean,
        axes["gm-ext-l"],
        axes["gm-ext-r"],
        pieces=get_source_pieces(global_mean["year"].values),
        split_year=split_year,
    )

    lat_gradient_pc = get_latitudinal_gradient_pc(bundle_dir)
    plot_extension_pieces(
        lat_gradient_pc,
        axes["lat-grad-pc-ext-l"],
        axes["lat-grad-pc-ext-r"],
        pieces=get_source_pieces(lat_gradient_pc["year"].values),
        split_year=split_year,
    )

    # The PC carries the units, so the EOF it multiplies is dimensionless
    lat_gradient_eof = linear_latitudinal_gradient_eof("dimensionless")
    plot_lat_gradient_eofs(
        lat_gradient_eof.assign_coords(eof=0).expand_dims({"eof": [0]}),
        axes["lat-grad-eof"],
    )

    native_resolution = xr.load_dataset(gas_dir / f"{GAS}_fifteen-degree_monthly.nc")
    max_year = int(native_resolution["year"].max())
    flying_carpet_mesh = plot_flying_carpet(
        native_resolution.sel(year=range(max_year - 9, max_year + 1)),
        axes["flying-carpet"],
    )

    gm_monthly = xr.load_dataset(gas_dir / f"{GAS}_global-mean_monthly.nc")
    gm_monthly = gm_monthly.assign_coords(lat=["Global"])
    hm_monthly = xr.load_dataset(gas_dir / f"{GAS}_hemispheric-mean_monthly.nc")
    hm_monthly = hm_monthly.assign_coords(
        lat=[
            "Southern hemisphere" if v == SH_LAT else "Northern hemisphere"
            for v in hm_monthly["lat"]
        ]
    )
    pda = xr.concat([gm_monthly, hm_monthly], "lat")
    pda = pda.rename({"lat": "Region"})
    max_year = int(pda["year"].max())
    plot_monthly_means(
        pda.sel(year=range(max_year - 4, max_year + 1)), ax=axes["monthly"]
    )

    gm_yearly = xr.load_dataset(gas_dir / f"{GAS}_global-mean_annual-mean.nc")
    gm_yearly = gm_yearly.assign_coords(lat=["Global"])
    hm_yearly = xr.load_dataset(gas_dir / f"{GAS}_hemispheric-mean_annual-mean.nc")
    hm_yearly = hm_yearly.assign_coords(
        lat=[
            "Southern hemisphere" if v == SH_LAT else "Northern hemisphere"
            for v in hm_yearly["lat"]
        ]
    )
    pda = xr.concat([gm_yearly, hm_yearly], "lat")
    pda = pda.rename({"lat": "Region"})
    plot_yearly_means(pda, ax_left=axes["yearly-l"], ax_right=axes["yearly-r"])

    add_colour_bar(
        fig,
        flying_carpet_mesh,
        cax=axes["flying-carpet-colour-bar"],
        # The panel's own vertical axis carries no label, so this says
        # both what is plotted and what its units are.
        label=label_name(
            f"{GAS} [{get_only_data_variable(native_resolution).attrs['units']}]"
        ),
        label_on_top=True,
    )

    label_panels(ROWS, axes, TITLES)
    # Last, because it needs to know how much room everything takes up
    lay_out_figure(fig, axes, ROWS)

    outfile.parent.mkdir(exist_ok=True, parents=True)
    logger.info(f"Writing {outfile}")
    fig.savefig(outfile)
    plt.close(fig)

    return outfile
