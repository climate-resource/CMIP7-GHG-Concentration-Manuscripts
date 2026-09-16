"""
Generation of the methods figure for the gases processed like C4F10

The pieces this figure shares with the CH4 methods figure
live in [local.historical_ghg_forcing_for_cmip7.plotting][].
What is here is the data this figure loads,
the panels it has and how they are laid out.

One function serves every gas in this group,
because the method is the same for all of them
(see the C4F10-like section of the manuscript's methods).
"""
# Differences from ch4
# - there is no observational network, only two sites' worth of data
#   from Droste et al. (2020), so the panel which says what went in
#   is those two timeseries and nothing else, and the binning,
#   the maps and the interpolation panels all go with it
# - the latitudinal gradient EOF is assumed rather than derived,
#   and the PC is whatever makes the gradient match the two sites
# - the global-mean is then whatever makes the two sites' values come back
#   once the latitudinal gradient is added to it, so it is derived
#   after the latitudinal gradient rather than before it
# - seasonality is assumed zero, so there is no seasonality panel
# - both extensions are the same two assumptions
#   (zero before Droste et al. starts, linear extrapolation after it ends),
#   so there is no regression panel either

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
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
    plot_global_mean_from_obs_network,
    plot_input_timeseries,
    plot_lat_gradient_pieces_from_obs_network,
    plot_monthly_means,
    plot_yearly_means,
)

C4F10_LIKE_GASES = (
    "c4f10",
    "c5f12",
    "c6f14",
    "c7f16",
    "cc4f8",
)
"""Gases which are processed the way C4F10 is

These are the step config IDs of the original run's
`calculate_c4f10_like_monthly_fifteen_degree_pieces` step,
i.e. exactly the gases this figure can be drawn for.
"""

DROSTE_FILE = Path("droste-et-al-2020") / "droste_et_al_2020.csv"
"""Where the input data lives, relative to the bundle's `data/interim` directory"""

DROSTE_LABEL = "Droste et al. (2020)"
"""How the input data's source is named to a reader"""

EXTENSION_SPLIT_MARGIN = 30
"""Years to leave before the input data starts when breaking an extended panel

Both extended components are zero for every year
before the input data's first year, and everything these panels are about
happens after it. So the axis is broken just before that year
rather than at some fixed year which suits no gas in particular:
the flat run-up goes in the left half and the rest gets the right half.
This is how much of the flat run-up is kept on the right,
so the reader can see the line was flat before the data starts.
"""

TITLES = {
    "inputs": f"Inputs: {DROSTE_LABEL}",
    "lat-grad-eof": "Assumed lat. gradient",
    "lat-grad-pc": "Derived lat. gradient PC",
    "gm": "Derived global-mean",
    "gm-ext": "Extended global-mean",
    "lat-grad-pc-ext": "Extended lat. gradient PC",
    "monthly": "Monthly spatial-means",
    "yearly": "Yearly spatial-means",
    "flying-carpet": "Native resolution",
}
"""Title of each panel"""

ROWS = (
    Row(
        panels=(Panel("inputs"),),
        height=2.5,
    ),
    Row(
        panels=(
            Panel("lat-grad-eof"),
            Panel("lat-grad-pc"),
            Panel("gm"),
        ),
        height=2.3,
    ),
    Row(
        panels=(
            Panel("gm-ext", broken=True, broken_split=BROKEN_SPLIT),
            Panel("lat-grad-pc-ext", broken=True, broken_split=BROKEN_SPLIT),
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

The rows follow the steps of the method,
so the panels are labelled in reading order.

- The inputs: what was measured, where and when.
- The components the inputs are decomposed into,
  in the order they are derived:
  the latitudinal gradient is assumed and then scaled to the inputs,
  and the global-mean is whatever is left over.
  There is no seasonality panel: seasonality is assumed zero for these gases.
- Extending each of those back and forward in time.
- The outputs, including the flying carpet,
  which is square and so sets its row's height.

Each row is laid out independently of the others,
see [local.historical_ghg_forcing_for_cmip7.layout][],
so rows need not have the same number of panels.
"""


def interim_dir(gas: str, bundle_dir: Path) -> Path:
    """
    Get the directory which holds a gas' interim data

    Parameters
    ----------
    gas
        Gas of interest

    bundle_dir
        Directory which holds the original run's bundle

    Returns
    -------
        Directory which holds `gas`' interim data
    """
    return bundle_dir / "data" / "interim" / gas


def get_droste_data(gas: str, bundle_dir: Path) -> pd.DataFrame:
    """
    Get a gas' input data

    Unlike the gases with an observational network,
    everything these gases are built from is a single source
    which the original run wrote into the bundle,
    so there is no notebook to re-run here.

    Parameters
    ----------
    gas
        Gas of interest

    bundle_dir
        Directory which holds the original run's bundle

    Returns
    -------
        The input data for `gas`

    Raises
    ------
    AssertionError
        The source has no data for `gas`
    """
    droste_file = bundle_dir / "data" / "interim" / DROSTE_FILE
    droste = pd.read_csv(droste_file)

    res = droste[droste["gas"] == gas]
    if res.empty:
        msg = f"No {DROSTE_LABEL} data for {gas=}, check {droste_file}"
        raise AssertionError(msg)

    return res


def get_extension_pieces(
    years: np.typing.NDArray[np.int64], input_years: np.typing.NDArray[np.int64]
) -> dict[str, np.typing.NDArray[np.int64]]:
    """
    Get the years which came from each source in an extended component

    Both extended components are extended the same way,
    because both are derived from the same source over the same years:
    zero before it starts (which is what the source itself says,
    and the only value a concentration and hence a gradient can take
    before the gas is manufactured),
    and a linear extrapolation after it ends.

    Parameters
    ----------
    years
        Years the extended component covers

    input_years
        Years the input data covers

    Returns
    -------
        Years which came from each source, labelled by source
    """
    return {
        "Assumed zero": years[years < input_years.min()],
        f"{DROSTE_LABEL} based": years[
            (years >= input_years.min()) & (years <= input_years.max())
        ],
        "Linear extrapolation": years[years > input_years.max()],
    }


def generate_c4f10_like_methods_figure(  # noqa: PLR0915
    gas: str,
    outfile: Path,
    bundle_dir: Path,
) -> Path:
    """
    Generate the methods figure for a gas which is processed like C4F10

    Parameters
    ----------
    gas
        Gas to draw the figure for

        Must be one of [`C4F10_LIKE_GASES`][].

    outfile
        File in which to write the figure

    bundle_dir
        Directory which holds the original run's bundle

    Returns
    -------
        `outfile`

    Raises
    ------
    AssertionError
        `gas` is not processed like C4F10,
        or its latitudinal gradient has more than one EOF,
        or its latitudinal gradient EOF is not the one we assume
    """
    if gas not in C4F10_LIKE_GASES:
        msg = f"{gas=} is not processed like C4F10, expected one of {C4F10_LIKE_GASES}"
        raise AssertionError(msg)

    gas_dir = interim_dir(gas, bundle_dir)

    droste = get_droste_data(gas, bundle_dir)
    input_years = np.sort(droste["year"].unique())

    fig, axes = create_figure(ROWS)

    plot_input_timeseries(droste, axes["inputs"])

    lat_gradient = xr.load_dataset(gas_dir / f"{gas}_allyears-lat-gradient-eofs-pcs.nc")
    expected_n_eofs = 1
    if len(lat_gradient["eof"]) != expected_n_eofs:
        msg = (
            f"{gas} has {len(lat_gradient['eof'])} EOFs, "
            f"these gases should only ever have {expected_n_eofs}"
        )
        raise AssertionError(msg)

    # The EOF is assumed rather than derived, and the panel says so,
    # so it is worth checking that the assumption in the file
    # is the one we describe.
    eofs = lat_gradient["eofs"]
    assumed_eof = linear_latitudinal_gradient_eof(eofs.attrs["units"])
    if not np.allclose(eofs.squeeze("eof").values, assumed_eof.values):
        msg = (
            f"{gas}'s stored latitudinal gradient EOF is not the linear EOF "
            "we say is assumed for these gases"
        )
        raise AssertionError(msg)

    # The PC is only derived for the years the input data covers;
    # every other year in this file came from the extension,
    # which has its own panel further down.
    plot_lat_gradient_pieces_from_obs_network(
        lat_gradient.sel(year=input_years),
        {
            "pcs": axes["lat-grad-pc"],
            "eofs": axes["lat-grad-eof"],
        },
    )

    global_mean_extended = xr.load_dataset(
        gas_dir / f"{gas}_global-annual-mean_allyears.nc"
    )
    plot_global_mean_from_obs_network(
        global_mean_extended.sel(year=input_years), axes["gm"]
    )

    split_year = int(input_years.min()) - EXTENSION_SPLIT_MARGIN

    global_mean_extended_da = get_only_data_variable(global_mean_extended)
    plot_extension_pieces(
        global_mean_extended_da,
        axes["gm-ext-l"],
        axes["gm-ext-r"],
        pieces=get_extension_pieces(
            global_mean_extended_da["year"].values, input_years
        ),
        split_year=split_year,
    )

    pcs_extended = lat_gradient["principal-components"].squeeze("eof", drop=True)
    plot_extension_pieces(
        pcs_extended,
        axes["lat-grad-pc-ext-l"],
        axes["lat-grad-pc-ext-r"],
        pieces=get_extension_pieces(pcs_extended["year"].values, input_years),
        split_year=split_year,
    )

    native_resolution = xr.load_dataset(gas_dir / f"{gas}_fifteen-degree_monthly.nc")
    max_year = int(native_resolution["year"].max())
    flying_carpet_mesh = plot_flying_carpet(
        native_resolution.sel(year=range(max_year - 9, max_year + 1)),
        axes["flying-carpet"],
    )

    gm_monthly = xr.load_dataset(gas_dir / f"{gas}_global-mean_monthly.nc")
    gm_monthly = gm_monthly.assign_coords(lat=["Global"])
    hm_monthly = xr.load_dataset(gas_dir / f"{gas}_hemispheric-mean_monthly.nc")
    sh_lat = -45.0
    hm_monthly = hm_monthly.assign_coords(
        lat=[
            "Southern hemisphere" if v == sh_lat else "Northern hemisphere"
            for v in hm_monthly["lat"]
        ]
    )
    pda = xr.concat([gm_monthly, hm_monthly], "lat")
    pda = pda.rename({"lat": "Region"})
    max_year = int(pda["year"].max())
    plot_monthly_means(
        pda.sel(year=range(max_year - 4, max_year + 1)), ax=axes["monthly"]
    )

    gm_yearly = xr.load_dataset(gas_dir / f"{gas}_global-mean_annual-mean.nc")
    gm_yearly = gm_yearly.assign_coords(lat=["Global"])
    hm_yearly = xr.load_dataset(gas_dir / f"{gas}_hemispheric-mean_annual-mean.nc")
    hm_yearly = hm_yearly.assign_coords(
        lat=[
            "Southern hemisphere" if v == sh_lat else "Northern hemisphere"
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
            f"{gas} [{get_only_data_variable(native_resolution).attrs['units']}]"
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
