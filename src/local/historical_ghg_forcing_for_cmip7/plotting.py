"""
Shared pieces of the historical GHG manuscript's methods figures

The CH4 and N2O methods figures tell the same story about different gases,
so they are the same figure with different data in it.
Everything which does not depend on which gas we are plotting lives here;
what is left in each gas' module is the data it loads,
the panels it has and how they are laid out.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import TYPE_CHECKING

import cartopy.crs as ccrs
import matplotlib.axes
import matplotlib.colorbar
import matplotlib.colors
import matplotlib.figure
import matplotlib.lines
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np
import openscm_units
import pandas as pd
import seaborn as sns
import xarray as xr

from local.xarray_time import convert_time_to_year_month, convert_year_month_to_time

if TYPE_CHECKING:
    import matplotlib.cm
    import matplotlib.collections

LAT_BIN_BOUNDS = np.arange(-90, 91, 15)
"""Bounds of the latitudinal bins used by the original run

This mirrors `local.binning.LAT_BIN_BOUNDS` in the original run.
"""

LAT_BIN_CENTRES = (LAT_BIN_BOUNDS[:-1] + LAT_BIN_BOUNDS[1:]) / 2.0
"""Centres of the latitudinal bins used by the original run
"""

LON_BIN_BOUNDS = np.arange(-180, 181, 60)
"""Bounds of the longitudinal bins used by the original run

This mirrors `local.binning.LON_BIN_BOUNDS` in the original run.
"""

LON_BIN_CENTRES = (LON_BIN_BOUNDS[:-1] + LON_BIN_BOUNDS[1:]) / 2.0
"""Centres of the longitudinal bins used by the original run
"""

LAT_AXIS_LIMITS = (-91, 91)
"""Limits of the latitude axis

Used by both the map and the counts panel, so the two can be read against
each other. A little room past the poles, so the polar stations have
somewhere to sit.
"""

MAP_ASPECT = (LAT_AXIS_LIMITS[1] - LAT_AXIS_LIMITS[0]) / 360.0
"""Height of a map's data box divided by its width

The maps are plate carrée, so one degree is the same length in either direction
and the map is as tall, relative to its width,
as its latitude range is relative to the 360 degrees of longitude.
"""

NETWORK_GROUPS = {
    "NOAA": "NOAA",
    "AGAGE": "AGAGE",
    "GAGE": "AGAGE",
    "ALE": "AGAGE",
}
"""How the networks are grouped for plotting

ALE, GAGE and AGAGE are successive instruments at what are largely
the same physical sites, so we show them as one network.
The group goes by the name of the current instrument,
because the full name of all three is wider
than the panel whose legend has to carry it.
"""

MOVING_SUFFIX = " (moving)"
"""Suffix which marks a network group whose stations move

Some networks report from ships as well as from fixed sites.
A ship is a different kind of observation from a fixed site
(it has no one location to put on the map),
so it is shown as its own group rather than folded in with the fixed sites.
"""

NETWORK_GROUP_COLOURS = {
    "NOAA": "#0072b2",
    "NOAA (moving)": "#3a3b3a",
    "AGAGE": "#d55e00",
}
"""Colour to use for each group of observational networks

These are from the Okabe-Ito palette, i.e. they are colour-blind safe.
The exception is the moving group, which is drawn in grey:
it is a qualifier on NOAA rather than a network in its own right,
and there is no fourth Okabe-Ito colour left which reads as such.

This covers every group either gas has.
A gas whose data has no shipboard observations
simply never asks for the moving group.
"""

NETWORK_GROUP_MARKERS = {
    "NOAA": "o",
    "NOAA (moving)": ".",
    "AGAGE": "^",
}
"""Marker to use for each group of observational networks

Groups are distinguished by marker as well as colour
so the figure still works in greyscale.
"""

NETWORK_GROUP_MARKER_SIZES = {
    "NOAA": 50.0,
    "NOAA (moving)": 10.0,
    "AGAGE": 25.0,
}
"""Marker size to use for each group of observational networks on the map

Several sites host both groups, so the markers are drawn
as concentric outlines, largest first, to keep both visible.
"""

NETWORK_GROUP_COLUMN = "network_group"
"""Column we add to hold the network group"""

LATITUDE_NORMALISATION = matplotlib.colors.Normalize(vmin=-90.0, vmax=90.0)
"""Range of latitudes the latitude colour map is stretched over

Fixed to the whole globe, rather than to the latitudes a gas happens
to have been measured at, so that a colour means the same latitude
in every panel and in both gases' figures.
"""

LATITUDE_COLOUR_MAP = "coolwarm"
"""Colour map used to show latitude

Diverging about the equator, and distinguishable with the common forms
of colour vision deficiency. Deliberately different from the colour map
used for the observation counts, so the two scales don't get confused.
"""

INTERPOLATION_COLOUR_MAP = "autumn_r"

LEGEND_MARKER_COLOUR = "0.35"
"""Colour to draw legend markers in where colour means something else

In the timeseries panel colour shows latitude, so the legend
can only speak about the marker shape.
"""


def ghg(pdf: pd.DataFrame, ghg_col: str = "gas") -> str:
    """
    Get GHG of a pandas.DataFrame
    """
    res_l = pdf[ghg_col].unique()
    if len(res_l) != 1:
        raise AssertionError(res_l)

    return res_l[0]


def unit(pdf: pd.DataFrame, unit_col: str = "unit") -> str:
    """
    Get unit of a pandas.DataFrame
    """
    res_l = pdf[unit_col].unique()
    if len(res_l) != 1:
        raise AssertionError(res_l)

    return res_l[0]


GHG_LABELS = {
    "co2": "CO$_2$",
    "ch4": "CH$_4$",
    "n2o": "N$_2$O",
}
"""How each gas' name is written when it is shown to a reader

Keyed by the name the gas goes by in the data,
which is what the figures have to hand.
"""


def label_name(ghg: str) -> str:
    """
    Get the label name for a given string

    Parameters
    ----------
    ghg
        String in which to write out any gas names

        This is a string rather than a bare gas name
        because the gas name is usually already part of a label
        by the time we get here.

    Returns
    -------
        `ghg`, with any gas names written out for a reader
    """
    res = ghg
    for old, new in GHG_LABELS.items():
        res = res.replace(old, new)

    return res


def get_only_data_variable(ds: xr.Dataset) -> xr.DataArray:
    """
    Get only data variable in a dataset
    """
    ds_vs = ds.data_vars
    if len(ds_vs) != 1:
        raise AssertionError(ds_vs)

    res = ds_vs[next(iter(ds.data_vars))]

    return res


def get_interpolated_input_coverage_info(
    all_data_with_bins: pd.DataFrame,
    interpolated_obs: xr.Dataset,
) -> dict[str, tuple[int, int]]:
    """
    Get interpolated input coverage information

    Returns the year and month for the time point
    with the most input points and the least input points.
    """
    # Reproduce what was done in the workflow.
    # Break this out if we need it elsewhere.
    interpolated_obs_nan_free = convert_year_month_to_time(
        convert_time_to_year_month(interpolated_obs).dropna("year")
    )
    nan_free_min_year = int(interpolated_obs_nan_free["time"].dt.year.min())
    nan_free_max_year = int(interpolated_obs_nan_free["time"].dt.year.max())
    all_data_with_bins_relevant = all_data_with_bins[
        (all_data_with_bins["year"] >= nan_free_min_year)
        & (all_data_with_bins["year"] <= nan_free_max_year)
    ].copy()
    all_data_with_bins_relevant["lat-lon-bin"] = (
        all_data_with_bins["lat_bin"].astype(str)
        + "__"
        + all_data_with_bins["lon_bin"].astype(str)
    )
    input_point_counts = all_data_with_bins_relevant.groupby(
        [
            "year",
            "month",
        ]
    )["lat-lon-bin"].count()

    res = {
        "most": input_point_counts.idxmax(),
        "least": input_point_counts.idxmin(),
    }

    return res


def get_decimal_year(indf: pd.DataFrame) -> pd.Series[float]:
    """
    Get the decimal year of each observation

    The data has separate year and month columns rather than a time axis.

    Parameters
    ----------
    indf
        Data for which to calculate the decimal year

    Returns
    -------
        Decimal year of each row of `indf`
    """
    return indf["year"] + (indf["month"] - 0.5) / 12.0


def add_network_group(
    indf: pd.DataFrame, surf_or_ship_col: str = "surf_or_ship"
) -> pd.DataFrame:
    """
    Add the network group of each observation

    Parameters
    ----------
    indf
        Data to which to add the network group

    surf_or_ship_col
        Column which says whether each observation was taken
        at the surface or from a ship

        Not every gas' data has this column.
        Where it is missing, we take every station to be a fixed site.

    Returns
    -------
        `indf`, with the network group added
    """
    out = indf.copy()
    out[NETWORK_GROUP_COLUMN] = out["network"].map(NETWORK_GROUPS)

    unmapped = out[NETWORK_GROUP_COLUMN].isna()
    if unmapped.any():
        msg = f"No network group for {sorted(set(out.loc[unmapped, 'network']))}"
        raise AssertionError(msg)

    if surf_or_ship_col not in out:
        # No shipboard observations for this gas, so every station is fixed
        # and there is nothing to split out.
        return out

    out[NETWORK_GROUP_COLUMN] = out[NETWORK_GROUP_COLUMN] + out[surf_or_ship_col].apply(
        get_moving_suffix
    )

    return pd.concat(
        [
            unify_station_metadata(station_df)
            for _, station_df in out.groupby(["network", "station"])
        ]
    )


def get_moving_suffix(surf_or_ship: str | float) -> str:
    """
    Get the suffix which marks whether a station moves

    Parameters
    ----------
    surf_or_ship
        Whether the observation was taken at the surface or from a ship

        Null means we were not told, in which case we assume a fixed site:
        a moving platform is the thing the data goes out of its way to flag.

    Returns
    -------
        Suffix to add to the network group's name
    """
    if pd.isnull(surf_or_ship):
        return ""

    if surf_or_ship == "surface":
        return ""

    if surf_or_ship == "shipboard":
        return MOVING_SUFFIX

    raise NotImplementedError(surf_or_ship)


def unify_station_metadata(
    station_df: pd.DataFrame, location_tolerance: float = 0.65
) -> pd.DataFrame:
    """
    Give one station one network group and, where we can, one location

    The data reports a network group and a location per observation
    rather than per station, and those can disagree within a station
    (a station which reports from both a fixed site and a ship,
    or a fixed site whose reported position wobbles by a fraction of a degree
    between reports).
    The map draws one marker per station, so it needs one answer for each.

    Parameters
    ----------
    station_df
        One station's observations

    location_tolerance
        How far a station's reported positions may spread, in degrees,
        and still be treated as one position

        Bigger than the wobble in a fixed site's reported position,
        much smaller than the distance a ship covers.

    Returns
    -------
        `station_df`, with its network group, and where possible its location,
        made consistent
    """
    out = station_df.copy()

    group_counts = out[NETWORK_GROUP_COLUMN].value_counts()
    if group_counts.shape[0] > 1:
        # Majority rules
        out[NETWORK_GROUP_COLUMN] = group_counts.index.values[0]

    location_cols = ["latitude", "longitude"]
    locations = out[location_cols]
    if (locations.std() < location_tolerance).all(axis=None):
        for col in location_cols:
            out.loc[:, col] = out[col].median()

    elif (
        locations.drop_duplicates().shape[0] > 1
        and not out[NETWORK_GROUP_COLUMN].str.endswith(MOVING_SUFFIX).all()
    ):
        msg = (
            "A station which does not move reports from more than one place: "
            f"{out['network'].iloc[0]} {out['station'].iloc[0]}"
        )
        raise NotImplementedError(msg)

    return out


def get_network_groups_largest_first(indf: pd.DataFrame) -> list[str]:
    """
    Get the network groups, ordered so the biggest markers are drawn first

    Parameters
    ----------
    indf
        Data from which to get the network groups

    Returns
    -------
        Network groups, largest marker first
    """
    return sorted(
        indf[NETWORK_GROUP_COLUMN].unique(),
        key=lambda group: NETWORK_GROUP_MARKER_SIZES[group],
        reverse=True,
    )


def add_compact_legend(
    ax: matplotlib.axes.Axes, fontsize: str = "x-small", **kwargs: object
) -> None:
    """
    Add a legend which takes as little of its panel as it can

    The panels are small, so a legend at matplotlib's defaults
    is a serious fraction of the panel it sits in,
    and a legend which is taller than its panel spills out over
    the panel's title and its neighbours.

    Parameters
    ----------
    ax
        Axes to add the legend to

    fontsize
        Size to draw the legend's text at

    **kwargs
        Passed on to `ax.legend`
    """
    ax.legend(
        **{
            "fontsize": fontsize,
            "framealpha": 0.9,
            "handletextpad": 0.2,
            "columnspacing": 1.0,
            "labelspacing": 0.3,
            "borderpad": 0.3,
            "borderaxespad": 0.3,
            # Anything the caller asks for wins over what we ask for here
            **kwargs,
        }
    )


def compact_existing_legend(ax: matplotlib.axes.Axes, **kwargs: object) -> None:
    """
    Re-draw a legend which seaborn has already put on a panel, compactly

    Seaborn draws its own legend, at matplotlib's default size,
    so the only way to get a compact one is to ask for it again
    with the handles seaborn worked out.

    Parameters
    ----------
    ax
        Axes whose legend to re-draw

        If it has no legend, nothing happens.

    **kwargs
        Passed on to [add_compact_legend][]
    """
    legend = ax.get_legend()
    if legend is None:
        return

    handles = legend.legend_handles
    labels = [text.get_text() for text in legend.get_texts()]
    title = legend.get_title().get_text()

    add_compact_legend(ax, handles=handles, labels=labels, title=title, **kwargs)
    ax.get_legend().get_title().set_fontsize("x-small")


def latitude_colour(latitude: float) -> tuple[float, float, float, float]:
    """
    Get the colour which stands for a latitude

    Parameters
    ----------
    latitude
        Latitude to get the colour of

    Returns
    -------
        Colour to draw `latitude` in
    """
    return plt.get_cmap(LATITUDE_COLOUR_MAP)(LATITUDE_NORMALISATION(latitude))


def add_latitude_legend(
    ax: matplotlib.axes.Axes,
    latitudes: np.typing.ArrayLike,
    ncols: int = 2,
    every: int = 2,
    **kwargs: object,
) -> None:
    """
    Add a legend which says which latitude each colour stands for

    Parameters
    ----------
    ax
        Axes to add the legend to

    latitudes
        Latitudes which appear on `ax`

        Listed north first, so the legend runs the same way up as a map does.

    ncols
        Number of columns to lay the legend out in

    every
        Show every nth latitude rather than all of them

        Colour stands for latitude here rather than for a category,
        so the legend is a scale, and a scale only has to be sampled
        finely enough to be read off:
        an entry for every latitudinal bin is more entries
        than a panel this size has room for,
        and says no more than every second one does.

    **kwargs
        Passed on to [add_compact_legend][]
    """
    add_compact_legend(
        ax,
        fontsize="xx-small",
        ncols=ncols,
        handles=[
            matplotlib.lines.Line2D(
                [],
                [],
                linestyle="none",
                marker="o",
                markersize=3,
                color=latitude_colour(latitude),
                label=f"{latitude:.0f}",
            )
            for latitude in sorted(np.asarray(latitudes), reverse=True)[::every]
        ],
        title=r"lat [$^{\circ}$N]",
        **kwargs,
    )
    ax.get_legend().get_title().set_fontsize("xx-small")


def plot_station_timeseries(  # noqa: PLR0913
    indf: pd.DataFrame,
    ax: matplotlib.axes.Axes,
    # In axis co-ords
    inset_x0: float = 0.65,
    inset_y0: float = 0.15,
    inset_width: float = 0.3,
    inset_height: float = 0.3,
    # consider making inset x and y injectable too
) -> matplotlib.collections.PathCollection:
    """
    Plot the monthly mean measured at each station

    This is a scatter plot rather than a line plot,
    so that months in which more than one station reported are visible.
    Colour shows the station's latitude and the marker shows its network,
    so the latitudinal spread of the network can be read directly.

    Parameters
    ----------
    indf
        Observational network data

    ax
        Axes on which to plot

    Returns
    -------
        The last scatter drawn, so a colour bar can be added for it
    """
    inset_xlim = (indf["year"].max() - 3, indf["year"].max())
    inset_values = indf[
        (indf["year"] >= inset_xlim[0]) & (indf["year"] <= inset_xlim[1])
    ]
    inset_values_range = inset_values["value"].max() - inset_values["value"].min()
    inset_ylim = (
        np.floor(max(inset_values["value"].min() - 0.05 * inset_values_range, 0.0)),
        np.ceil(inset_values["value"].max() + 0.05 * inset_values_range),
    )
    ax_inset = ax.inset_axes(
        [inset_x0, inset_y0, inset_width, inset_height],
        xlim=inset_xlim,
        ylim=inset_ylim,
        # xticklabels=[],
        # yticklabels=[],
    )
    ax_inset.tick_params(labelsize="small")
    ax_inset.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))

    ax.indicate_inset_zoom(ax_inset, edgecolor="black")

    scatter = None
    for group, group_df in indf.groupby(NETWORK_GROUP_COLUMN):
        for axh in [ax, ax_inset]:
            scatter = axh.scatter(
                get_decimal_year(group_df),
                group_df["value"],
                c=group_df["latitude"],
                cmap=LATITUDE_COLOUR_MAP,
                norm=LATITUDE_NORMALISATION,
                marker=NETWORK_GROUP_MARKERS[group],
                s=10.0,
                alpha=0.4,
                linewidths=0.0,
            )

    if scatter is None:
        msg = "No data to plot"
        raise AssertionError(msg)

    units = unit(indf)
    # gas = ghg(indf)

    ax.set_ylabel(label_name(f"[{units}]"), fontsize="small")
    ax.set_xlabel("year", fontsize="small")
    ax.tick_params(labelsize="small")

    # Colour is doing latitude here, so the legend can only speak about shape
    add_compact_legend(
        ax,
        loc="upper left",
        handles=[
            matplotlib.lines.Line2D(
                [],
                [],
                linestyle="none",
                marker=NETWORK_GROUP_MARKERS[group],
                color=LEGEND_MARKER_COLOUR,
                markersize=5,
                label=group,
            )
            for group in sorted(indf[NETWORK_GROUP_COLUMN].unique())
        ],
    )

    return scatter


def plot_station_locations(
    indf: pd.DataFrame,
    ax: matplotlib.axes.Axes,
) -> None:
    """
    Plot the location of each station in the observational network

    Parameters
    ----------
    indf
        Observational network data

    ax
        Axes on which to plot

        This must have been created with a cartopy projection.
    """
    ax.coastlines(linewidth=0.4, color="0.55")
    ax.set_global()

    # The bins the observations are binned into,
    # so this panel can be read against the counts panel
    for lat_bound in LAT_BIN_BOUNDS:
        ax.axhline(lat_bound, linewidth=0.4, color="0.85", zorder=0)

    for lon_bound in LON_BIN_BOUNDS:
        ax.axvline(lon_bound, linewidth=0.4, color="0.85", zorder=0)

    stations = indf[
        [NETWORK_GROUP_COLUMN, "station", "latitude", "longitude"]
    ].drop_duplicates()

    # Biggest markers first, so co-located networks nest rather than hide each other
    for group in get_network_groups_largest_first(stations):
        group_stations = stations[stations[NETWORK_GROUP_COLUMN] == group]
        ax.scatter(
            group_stations["longitude"],
            group_stations["latitude"],
            transform=ccrs.PlateCarree(),
            facecolors="none",
            edgecolors=NETWORK_GROUP_COLOURS[group],
            marker=NETWORK_GROUP_MARKERS[group],
            s=NETWORK_GROUP_MARKER_SIZES[group],
            linewidths=1.1,
            label=group,
            zorder=3,
            # Cartopy clips to the projection boundary,
            # which would cut the South Pole station in half
            clip_on=False,
        )

    ax.set_yticks(LAT_BIN_BOUNDS[::2], crs=ccrs.PlateCarree())
    ax.set_ylabel(r"latitude [$^{\circ}$N]", fontsize="small")
    ax.set_xticks(LON_BIN_BOUNDS[::2], crs=ccrs.PlateCarree())
    ax.set_xlabel(r"longitude [$^{\circ}$E]", fontsize="small")
    # The map's aspect is fixed, so it can't fill its cell.
    # Anchoring it to the top keeps it up against its panel label.
    ax.set_anchor("N")
    ax.tick_params(labelsize="small")
    ax.set_ylim(LAT_AXIS_LIMITS)
    # On the panel's title line, which has room to spare and is the only
    # place around this panel that does: the map cannot be drawn over
    # (that would cover stations), its row is only as tall as the map itself
    # so there is no band under it, and the gap beside it is a gap between
    # two columns rather than space belonging to this panel.
    add_compact_legend(
        ax,
        fontsize="xx-small",
        loc="lower right",
        bbox_to_anchor=(1.0, 1.0),
        # All on one row, so the legend is no taller than the title line it shares.
        # The maps are wide enough that it still stays clear of the title.
        ncols=len(get_network_groups_largest_first(stations)),
        handlelength=1.0,
        handletextpad=0.3,
    )
    # The gap beside the map is there whether or not we put the legend in it,
    # so the legend is placed by hand and kept out of the layout engine's sums.
    # Left in them, it would ask for its width from the next panel's column,
    # which every panel in that column would then pay for.
    ax.get_legend().set_in_layout(False)


def plot_observation_counts(
    indf: pd.DataFrame,
    ax: matplotlib.axes.Axes,
) -> matplotlib.collections.QuadMesh:
    """
    Plot the number of input data points in each month-latitudinal bin

    The networks report monthly averages, so each row of the data
    is one station's average for one month, and that is what is counted here.
    Many more measurements sit underneath each of those averages,
    but the networks do not report them, so we cannot count them.
    Counts are summed over all longitudes.

    Parameters
    ----------
    indf
        Observational network data

    ax
        Axes on which to plot

    Returns
    -------
        The mesh which was drawn, so a colour bar can be added for it
    """
    # The bins the observations are binned into
    for lat_bound in LAT_BIN_BOUNDS:
        ax.axhline(lat_bound, linewidth=0.4, color="0.85", zorder=3)

    counts = (
        indf.groupby(["year", "month", "lat_bin"]).size().rename("count").reset_index()
    )
    counts["decimal_year"] = get_decimal_year(counts)

    grid = counts.pivot_table(
        index="lat_bin", columns="decimal_year", values="count", fill_value=0
    )
    # Make sure every latitudinal bin has a row, even the empty ones
    lat_bin_centres = (LAT_BIN_BOUNDS[:-1] + LAT_BIN_BOUNDS[1:]) / 2
    grid = grid.reindex(lat_bin_centres, fill_value=0)

    # Every month in the record, so gaps show up as gaps rather than being closed up
    month_starts = np.arange(
        np.floor(counts["decimal_year"].min()),
        np.ceil(counts["decimal_year"].max()) - 1 / 24,
        1 / 12,
    )
    grid = grid.reindex(
        columns=grid.columns.union(month_starts + 1 / 24), fill_value=0
    ).sort_index(axis="columns")

    x_bounds = np.append(grid.columns.to_numpy() - 1 / 24, grid.columns[-1] + 1 / 24)

    # Counts are small integers, so use a discrete scale rather than a continuous one
    max_count = int(grid.to_numpy().max())
    mesh = ax.pcolormesh(
        x_bounds,
        LAT_BIN_BOUNDS,
        np.ma.masked_equal(grid.to_numpy(), 0),
        cmap=plt.get_cmap("YlOrRd", max_count),
        norm=matplotlib.colors.BoundaryNorm(np.arange(0.5, max_count + 1.0), max_count),
        shading="flat",
    )
    ax.set_yticks(LAT_BIN_BOUNDS[::2])
    ax.set_ylim(LAT_AXIS_LIMITS)
    ax.set_ylabel(r"latitude [$^{\circ}$N]", fontsize="small")
    ax.tick_params(labelsize="small")

    return mesh


def plot_coverage_and_interpolated(
    input_data: pd.DataFrame,
    interpolated: xr.Dataset,
    year_month: tuple[int, int],
    ax: matplotlib.axes.Axes,
    # return type hint is wrong
) -> None:
    """
    Plot the interpolated values and the coverage of input data
    """
    ax.coastlines(linewidth=0.6, color="0.3", zorder=2.0)
    ax.set_global()

    for lat_bound in LAT_BIN_BOUNDS:
        ax.axhline(lat_bound, linewidth=0.4, color="0.85", zorder=0)

    for lon_bound in LON_BIN_BOUNDS:
        ax.axvline(lon_bound, linewidth=0.4, color="0.85", zorder=0)

    lon_grid, lat_grid = np.meshgrid(
        LON_BIN_CENTRES,
        LAT_BIN_CENTRES,
    )

    interpolated_ym = convert_time_to_year_month(interpolated).sel(
        year=year_month[0], month=year_month[1]
    )
    interpolated_ym_da = get_only_data_variable(interpolated_ym)
    # interpolated_ym_vs = interpolated_ym.data_vars
    # if len(interpolated_ym_vs) != 1:
    #     raise AssertionError
    #
    # interpolated_ym_da = interpolated_ym_vs[next(iter(interpolated_ym.data_vars))]

    mesh = ax.pcolormesh(
        lon_grid,
        lat_grid,
        interpolated_ym_da.T,
        shading="auto",
        cmap=INTERPOLATION_COLOUR_MAP,
    )

    input_data_ym = input_data[
        (input_data["year"] == year_month[0]) & (input_data["month"] == year_month[1])
    ]
    ax.scatter(
        input_data_ym["lon_bin"],
        input_data_ym["lat_bin"],
        transform=ccrs.PlateCarree(),
        c="k",
        marker="o",
        s=10.0,
        label="Input point",
        zorder=3,
    )

    ax.set_yticks(LAT_BIN_BOUNDS[::2], crs=ccrs.PlateCarree())
    ax.set_ylabel(r"latitude [$^{\circ}$N]", fontsize="small")
    ax.set_xticks(LON_BIN_BOUNDS[::2], crs=ccrs.PlateCarree())
    ax.set_xlabel(r"longitude [$^{\circ}$E]", fontsize="small")
    # The map's aspect is fixed, so it can't fill its cell.
    # Anchoring it to the top keeps it up against its panel label.
    ax.set_anchor("N")
    ax.tick_params(labelsize="small")
    # Same as the locations map, so all the maps are the same shape
    ax.set_ylim(LAT_AXIS_LIMITS)

    return mesh


def plot_global_mean_from_obs_network(gm: xr.Dataset, ax: matplotlib.axes.Axes) -> None:
    """
    Plot global-mean derived from the observational network
    """
    gm_da = get_only_data_variable(gm)
    ax.scatter(
        gm_da["year"].values.squeeze(),
        gm_da.values.squeeze(),
        s=15,
    )
    ax.set_ylabel(f"[{gm_da.attrs['units']}]", fontsize="small")
    ax.tick_params(labelsize="small")


def plot_seasonality_from_obs_network(
    seasonality: xr.Dataset,
    ax: matplotlib.axes.Axes,
    assumed_units: str,
    # assumed_ghg: str,
) -> matplotlib.axes.Axes:
    """
    Plot seasonality derived from the observational network

    There is one series here for each of the twelve latitudinal bins.
    Latitude is shown with the same colour map, over the same range,
    as the timeseries panel uses, so a colour means the same latitude
    everywhere in the figure.
    """
    seasonality_da = get_only_data_variable(seasonality)
    pdf = seasonality_da.to_pandas().stack().rename("value").to_frame().reset_index()
    ax.scatter(
        pdf["month"],
        pdf["value"],
        c=pdf["lat"],
        cmap=LATITUDE_COLOUR_MAP,
        norm=LATITUDE_NORMALISATION,
        s=30.0,
        linewidths=0.0,
    )
    ax.set_ylabel(f"[{seasonality_da.attrs['units']}]", fontsize="small")
    ax.set_xlabel("month", fontsize="small")
    ax.set_xticks(np.arange(1, 12 + 1, 3))
    ax.tick_params(labelsize="small")
    add_latitude_legend(ax, pdf["lat"].unique(), loc="best")

    return ax


def plot_lat_gradient_pieces_from_obs_network(
    lat_gradient_info: xr.Dataset,
    axes: Mapping[str, matplotlib.axes.Axes],
    pcs_name: str = "principal-components",
    eofs_name: str = "eofs",
) -> dict[str, matplotlib.axes.Axes]:
    """
    Plot seasonality derived from the observational network
    """
    da_pcs = lat_gradient_info[pcs_name]
    pdf_pcs = da_pcs.to_pandas().stack().rename("value").to_frame().reset_index()
    sns.scatterplot(
        pdf_pcs,
        x="year",
        y="value",
        hue="eof",
        ax=axes["pcs"],
    )
    axes["pcs"].set_ylabel(f"[{da_pcs.attrs['units']}]", fontsize="small")
    axes["pcs"].set_xlabel("year", fontsize="small")
    axes["pcs"].tick_params(labelsize="small")
    compact_existing_legend(axes["pcs"], loc="best")

    da_eofs = lat_gradient_info[eofs_name]
    pdf_eofs = da_eofs.to_pandas().stack().rename("value").to_frame().reset_index()
    sns.scatterplot(
        pdf_eofs,
        x="value",
        y="lat",
        hue="eof",
        ax=axes["eofs"],
    )
    axes["eofs"].set_yticks(LAT_BIN_BOUNDS[::2])
    axes["eofs"].set_ylabel(r"latitude [$^{\circ}$N]", fontsize="small")
    axes["eofs"].set_xlabel(f"[{da_eofs.attrs['units']}]", fontsize="small")
    axes["eofs"].tick_params(labelsize="small")
    compact_existing_legend(axes["eofs"], loc="best")

    return axes


def add_break_lines_and_setup(  # noqa: PLR0913
    ax_left: matplotlib.axes.Axes,
    ax_right: matplotlib.axes.Axes,
    min_year: int,
    split_year: int,
    max_year: int,
    units: str,
    legend_loc: str = "best",
) -> None:
    """
    Add break lines to axes and do other general setup for broken axes

    Parameters
    ----------
    ax_left
        Left half of the broken axis

    ax_right
        Right half of the broken axis

    min_year
        First year the left half covers

    split_year
        Year the axis is broken at

    max_year
        Last year the right half covers

    units
        Units of the values plotted

    legend_loc
        Where to put the pair's legend

        Passed on to `ax.legend`, so `"best"` leaves it to matplotlib.
        Worth naming where matplotlib's answer is not a good one:
        it places the legend where the data is thinnest,
        which on a panel whose data hugs one edge
        is on top of the axis' own tick labels.
    """
    # One legend for the pair, on the left half:
    # the two halves are one panel as far as a reader is concerned.
    if ax_right.get_legend() is not None:
        ax_right.get_legend().remove()

    if ax_left.get_legend() is None:
        add_compact_legend(ax_left, loc=legend_loc)
    else:
        compact_existing_legend(ax_left, loc=legend_loc)

    ax_left.set_ylabel(f"[{units}]", fontsize="small")
    ax_right.set_ylabel("")
    for ax in (ax_left, ax_right):
        ax.set_xlabel("year", fontsize="small")
        ax.tick_params(labelsize="small")
    # ax.set_xlabel("year")
    ax_left.set_xlim(xmin=min_year, xmax=split_year)
    ax_right.set_xlim(xmin=split_year, xmax=max_year)

    # Hide the spines
    ax_left.spines["right"].set_visible(False)
    ax_right.spines["left"].set_visible(False)

    # The two halves are one panel, so they share one vertical scale,
    # and the right half's copy of it would say nothing the left half
    # has not already said. Hiding it is not only tidier:
    # a panel's tick labels are the most expensive thing about it
    # as far as the layout engine is concerned,
    # and the halves are the only panels which can give that cost up.
    ax_right.set_ylim(ax_left.get_ylim())
    ax_right.yaxis.set_visible(False)

    # Potential option for cutout lines.
    # The important
    # thing to know here is that in axes coordinates, which are always
    # between 0-1, spine endpoints are at these locations (0, 0), (0, 1),
    # (1, 0), and (1, 1).  Thus, we just need to put the diagonals in the
    # appropriate corners of each of our axes, and so long as we use the
    # right transform and disable clipping.

    clear_ticks_near_break(ax_left, at="right")
    clear_ticks_near_break(ax_right, at="left")

    d = 0.015  # how big to make the diagonal lines in axes coordinates
    # arguments to pass plot, just so we don't keep repeating them
    kwargs = dict(transform=ax_left.transAxes, color="k", clip_on=False)
    ax_left.plot((1 - d, 1 + d), (-d, +d), **kwargs)
    ax_left.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)

    kwargs.update(transform=ax_right.transAxes)  # switch to the bottom axes
    ax_right.plot((-d, +d), (1 - d, 1 + d), **kwargs)
    ax_right.plot((-d, +d), (-d, +d), **kwargs)


def clear_ticks_near_break(
    ax: matplotlib.axes.Axes, at: str, keep_clear: float = 0.06
) -> None:
    """
    Drop the tick labels which sit right against a broken axis' break

    The two halves of a broken axis are drawn hard up against each other,
    so a tick label at the very end of one half
    lands on top of a tick label at the very start of the other.

    Parameters
    ----------
    ax
        Half of a broken axis

    at
        Which end of `ax` the break is at, `"left"` or `"right"`

    keep_clear
        How much of the axis to keep clear of tick labels,
        as a fraction of the axis' range

    Raises
    ------
    ValueError
        `at` is neither `"left"` nor `"right"`
    """
    if at not in ("left", "right"):
        msg = f"at must be 'left' or 'right', got {at!r}"
        raise ValueError(msg)

    x_min, x_max = ax.get_xlim()
    keep_clear_from = keep_clear * (x_max - x_min)

    if at == "right":
        keep = [tick for tick in ax.get_xticks() if tick < x_max - keep_clear_from]
    else:
        keep = [tick for tick in ax.get_xticks() if tick > x_min + keep_clear_from]

    ax.set_xticks(keep)
    # set_xticks resets the limits it was given the ticks for
    ax.set_xlim(x_min, x_max)


def plot_global_mean_extension(
    gm: xr.Dataset,
    ax_left: matplotlib.axes.Axes,
    ax_right: matplotlib.axes.Axes,
    input_sources=Mapping[str, pd.DataFrame],
    split_year: int = 1950,
) -> None:
    """
    Plot global-mean derived from the observational network
    """
    gm_da = get_only_data_variable(gm)

    for i, ax in enumerate((ax_left, ax_right)):
        ax.plot(
            gm_da["year"].values.squeeze(),
            gm_da.values.squeeze(),
            label="Extended global-mean" if i < 1 else None,
            linewidth=2,
            # s=30,
        )

        for label, pdf in input_sources.items():
            ax.plot(
                pdf["year"],
                pdf["value"],
                label=label if i < 1 else None,
                linewidth=2,
                # s=30,
            )

    add_break_lines_and_setup(
        ax_left,
        ax_right,
        gm_da["year"].min(),
        split_year,
        gm_da["year"].max(),
        gm_da.attrs["units"],
    )


def plot_pc_timeseries_regression(  # noqa: PLR0913
    lat_grad_info: xr.Dataset,
    timeseries_data: xr.Dataset,
    timeseries_name: str,
    regression_info: dict[str, tuple[float, str]],
    ax: matplotlib.axes.Axes,
    x_unit: str,
    pcs_key: str = "principal-components",
    eof: int = 0,
    ur=openscm_units.unit_registry,
) -> None:
    """
    Plot the regression between a PC and a timeseries
    """
    pc_da = lat_grad_info[pcs_key].sel(eof=eof)

    common_years = np.intersect1d(pc_da["year"], timeseries_data["year"])

    timeseries_da = get_only_data_variable(timeseries_data.sel(year=common_years))
    timeseries_da.to_pandas()
    timeseries_da_units = timeseries_da.attrs["units"]
    conversion_factor = ur(timeseries_da_units).to(x_unit).m
    timeseries_values = timeseries_da.values * conversion_factor

    pc_units = pc_da.attrs["units"]
    pc_values = pc_da.sel(year=common_years).values

    ax.scatter(
        x=timeseries_values,
        y=pc_values,
        label="raw data",
        marker="x",
        s=30,
        color="tab:blue",
        alpha=0.7,
    )
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()

    regression_gradient = (
        ur.Quantity(regression_info["m"][0], regression_info["m"][1])
        .to(f"{pc_units} / ({x_unit})")
        .m
    )
    regression_y_int = (
        ur.Quantity(regression_info["c"][0], regression_info["c"][1]).to(pc_units).m
    )
    ax.axline(
        xy1=(0, regression_y_int),
        slope=regression_gradient,
        label="regression",
        linestyle="-",
        color="tab:orange",
        alpha=0.9,
    )
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    ax.set_xlabel(f"{timeseries_name} [{x_unit}]", fontsize="small")
    ax.set_ylabel(f"PC{eof} [{pc_units}]", fontsize="small")
    ax.tick_params(labelsize="small")

    add_compact_legend(ax, loc="best")


def plot_pcs_extended(  # noqa: PLR0913
    eof_pieces: xr.Dataset,
    ax_left: matplotlib.axes.Axes,
    ax_right: matplotlib.axes.Axes,
    pcs_key: str = "principal-components",
    split_year: int = 1950,
    # Should be numpy array of int, anyway
    pieces: dict[int, dict[str, list[int]]] | None = None,
) -> None:
    """
    Plot extended latitudinal gradient PCs
    """
    pcs_da = eof_pieces[pcs_key]
    pcs_df = pcs_da.to_pandas().stack().rename("value").to_frame().reset_index()

    for eof, info in pieces.items():
        for source, years in info.items():
            pcs_df.loc[
                (pcs_df["eof"] == eof) & (pcs_df["year"].isin(years)), "source"
            ] = source

    if pcs_df["source"].isnull().any(axis=None):
        raise AssertionError

    for i, ax in enumerate((ax_left, ax_right)):
        sns.scatterplot(
            pcs_df,
            x="year",
            y="value",
            # hue="eof",
            # style="source",
            style="eof",
            hue="source",
            ax=ax,
            s=25,
            edgecolor=None,
            alpha=0.7,
        )

    add_break_lines_and_setup(
        ax_left,
        ax_right,
        pcs_da["year"].min(),
        split_year,
        pcs_da["year"].max(),
        pcs_da.attrs["units"],
    )

    legend = ax_left.get_legend()
    if legend is None:
        raise AssertionError

    handles = legend.legend_handles
    labels = [text.get_text() for text in legend.get_texts()]
    # Add fake handles and legends for other EOFs so 2 col splits as we want.
    for _ in range(len(pcs_df["source"].unique()) - len(pcs_df["eof"].unique())):
        handles.append(matplotlib.lines.Line2D([], [], color="none", label=""))
        labels.append("")

    title = legend.get_title().get_text()

    add_compact_legend(ax_left, handles=handles, labels=labels, title=title, ncols=2)
    ax_left.get_legend().get_title().set_fontsize("x-small")


def plot_flying_carpet(
    native_resolution: xr.Dataset, ax: matplotlib.axes.Axes, tick_nearest: int = 5
) -> None:
    """
    Plot flying carpet
    """
    tmp = convert_year_month_to_time(get_only_data_variable(native_resolution))
    tmp = tmp.assign_coords(
        time=tmp["time"].dt.year + tmp["time"].dt.month / 12 - 1 / 24
    )
    mesh = tmp.plot.surface(
        x="time",
        y="lat",
        ax=ax,
        cmap=INTERPOLATION_COLOUR_MAP,
        levels=30,
        add_colorbar=False,
        # alpha=0.7,
    )
    ax.view_init(15, -135, 0)  # type: ignore
    # No label on the vertical axis. Matplotlib works out which way up to
    # write a 3D axis' label from the shape of the axes' box, and on this
    # panel it gets the vertical one upside down as often as not.
    # The colour bar carries what the label would have said.
    ax.set_zlabel("")  # type: ignore
    ax.set_ylabel(r"latitude [$^{\circ}$N]", fontsize="small")
    ax.set_yticks([-45, 45])
    ax.set_xlabel("time", fontsize="small")
    max_year = int(native_resolution["year"].max())
    other_tick = tick_nearest * np.floor(max_year / tick_nearest) - 5
    ax.set_xticks([other_tick, max_year])
    ax.tick_params(labelsize="small", pad=0.0)

    return mesh


def plot_monthly_means(
    pda: xr.Dataset,
    ax: matplotlib.axes.Axes,
) -> None:
    """
    Plot monthly means
    """
    tmp = convert_year_month_to_time(get_only_data_variable(pda))
    tmp = tmp.assign_coords(
        time=tmp["time"].dt.year + tmp["time"].dt.month / 12 - 1 / 24
    )
    pdf_l = []
    for region, rda in tmp.groupby("Region"):
        tmp_df = (
            rda.sel(Region=region).to_pandas().rename("value").to_frame().reset_index()
        )
        tmp_df["Region"] = region
        pdf_l.append(tmp_df)

    pdf = pd.concat(pdf_l)

    sns.scatterplot(
        pdf,
        x="time",
        y="value",
        hue="Region",
        ax=ax,
        s=15,
        alpha=0.7,
        edgecolor=None,
    )

    ax.set_ylabel(f"[{tmp.attrs['units']}]", fontsize="small")
    ax.set_xlabel("time", fontsize="small")
    ax.tick_params(labelsize="small")
    ax.set_xlim(pda["year"].min(), pda["year"].max())
    ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
    compact_existing_legend(ax, loc="best")


def plot_yearly_means(
    pda: xr.Dataset,
    ax_left: matplotlib.axes.Axes,
    ax_right: matplotlib.axes.Axes,
    split_year: int = 1850,
) -> None:
    """
    Plot global-mean derived from the observational network
    """
    tmp = get_only_data_variable(pda)
    pdf_l = []
    for region, rda in tmp.groupby("Region"):
        tmp_df = (
            rda.sel(Region=region).to_pandas().rename("value").to_frame().reset_index()
        )
        tmp_df["Region"] = region
        pdf_l.append(tmp_df)

    pdf = pd.concat(pdf_l)

    for i, ax in enumerate((ax_left, ax_right)):
        sns.scatterplot(
            pdf,
            x="year",
            y="value",
            hue="Region",
            ax=ax,
            s=15,
            alpha=0.7,
            edgecolor=None,
            # legend=False,
        )

    add_break_lines_and_setup(
        ax_left,
        ax_right,
        pda["year"].min(),
        split_year,
        pda["year"].max(),
        tmp.attrs["units"],
        # The record is flat and low until its last century, so the top of
        # the left half is the only part of this panel with room in it,
        # and centring keeps the legend off the axis' tick labels.
        legend_loc="upper center",
    )


def thin_ticks(ticks: np.typing.ArrayLike, max_ticks: int) -> np.typing.NDArray:
    """
    Take every nth tick, so that no more than `max_ticks` are left

    A colour bar which stands for a count wants a tick per count,
    right up until there are more counts than the colour bar has room for,
    at which point the labels run into each other
    and the colour bar says less than it would with a handful of ticks.

    Parameters
    ----------
    ticks
        Ticks we would show if there were room for all of them

    max_ticks
        The most ticks to leave

    Returns
    -------
        Every nth tick of `ticks`, for the smallest n which leaves
        no more than `max_ticks` of them
    """
    ticks = np.asarray(ticks)
    if ticks.size <= max_ticks:
        return ticks

    return ticks[:: int(np.ceil(ticks.size / max_ticks))]


def add_colour_bar(  # noqa: PLR0913
    fig: matplotlib.figure.Figure,
    mappable: matplotlib.cm.ScalarMappable,
    cax: matplotlib.axes.Axes,
    label: str,
    ticks: np.typing.ArrayLike | None = None,
    max_ticks: int = 6,
    label_on_top: bool = False,
    **kwargs: object,
) -> matplotlib.colorbar.Colorbar:
    """
    Add a colour bar to the figure

    Parameters
    ----------
    fig
        Figure to add the colour bar to

    mappable
        Mappable to draw the colour bar for

    cax
        Axes in which to draw the colour bar

        Where this sits is up to the figure's layout,
        see [local.historical_ghg_forcing_for_cmip7.layout][].

    label
        Label for the colour bar

    ticks
        Ticks to show on the colour bar

        If not supplied, the ticks are left to matplotlib.

    max_ticks
        The most ticks to actually show out of `ticks`

        See [thin_ticks][].

    label_on_top
        Put the label above the colour bar rather than beside it

        Useful where there is no width to spare beside the colour bar.

    **kwargs
        Passed on to `fig.colorbar`

    Returns
    -------
        The colour bar which was added
    """
    colour_bar = fig.colorbar(mappable, cax=cax, **kwargs)
    if ticks is not None:
        colour_bar.set_ticks(thin_ticks(ticks, max_ticks=max_ticks))

    if label_on_top:
        colour_bar.ax.set_title(label, fontsize="small", loc="left")
    else:
        colour_bar.set_label(label, fontsize="small")

    colour_bar.ax.tick_params(labelsize="small")

    return colour_bar
