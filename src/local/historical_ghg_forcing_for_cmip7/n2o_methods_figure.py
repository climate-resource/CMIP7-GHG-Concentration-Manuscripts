"""
Generation of the N2O methods figure

The figure shows the observational network which underpins
the N2O concentrations: where the stations are, what they measured over time,
and how many observations fall in each month-latitude bin.

The data behind it was never saved by the original run
(only the binned averages were), so we re-run the notebook which built it
with an extra save call bolted on.
See [local.cmip_ghg_generation][] for how that works.
"""

from __future__ import annotations

import string
from pathlib import Path

import cartopy.crs as ccrs
import matplotlib.axes
import matplotlib.colors
import matplotlib.figure
import matplotlib.lines
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from loguru import logger

from local.cmip_ghg_generation import (
    DEFAULT_BUNDLE_DIR,
    DEFAULT_ORIGINAL_RUN_NOTEBOOKS_DIR,
    MODIFIED_NOTEBOOKS_DIR,
    ensure_bundle_available,
    ensure_bundle_environment,
    ensure_executed_notebook_available,
    run_notebook_from_bundle_dir,
    write_modified_notebook,
)
from local.xarray_time import convert_time_to_year_month, convert_year_month_to_time

ALL_DATA_WITH_BINS_FILE = Path("manuscript-outputs") / "n2o_all-data-with-bins.csv"
"""Where the re-run notebook saves the data we want

Relative to the bundle's root directory,
because that is the notebook's working directory.
"""

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

NETWORK_GROUPS = {
    "NOAA": "NOAA",
    "AGAGE": "AGAGE/GAGE/ALE",
    "GAGE": "AGAGE/GAGE/ALE",
    "ALE": "AGAGE/GAGE/ALE",
}
"""How the networks are grouped for plotting

ALE, GAGE and AGAGE are successive instruments at what are largely
the same physical sites, so we show them as one network.
"""

NETWORK_GROUP_COLOURS = {
    "NOAA": "#0072b2",
    "AGAGE/GAGE/ALE": "#d55e00",
}
"""Colour to use for each group of observational networks

These are from the Okabe-Ito palette, i.e. they are colour-blind safe.
"""

NETWORK_GROUP_MARKERS = {
    "NOAA": "o",
    "AGAGE/GAGE/ALE": "^",
}
"""Marker to use for each group of observational networks

Groups are distinguished by marker as well as colour
so the figure still works in greyscale.
"""

NETWORK_GROUP_MARKER_SIZES = {
    "NOAA": 50.0,
    "AGAGE/GAGE/ALE": 25.0,
}
"""Marker size to use for each group of observational networks on the map

Several sites host both groups, so the markers are drawn
as concentric outlines, largest first, to keep both visible.
"""

NETWORK_GROUP_COLUMN = "network_group"
"""Column we add to hold the network group"""

LATITUDE_COLOUR_MAP = "coolwarm"
"""Colour map used to show latitude

Diverging about the equator, and distinguishable with the common forms
of colour vision deficiency. Deliberately different from the colour map
used for the observation counts, so the two scales don't get confused.
"""

INTERPOLATION_COLOUR_MAP = "autumn"

LEGEND_MARKER_COLOUR = "0.35"
"""Colour to draw legend markers in where colour means something else

In the timeseries panel colour shows latitude, so the legend
can only speak about the marker shape.
"""

PANELS = (
    ("timeseries", "Observational network values"),
    ("locations", "Obs. locations"),
    ("counts", "Obs. counts"),
    ("interpolated-most", "Interpolation: most inputs"),
    ("interpolated-least", "Interpolation: fewest inputs"),
    ("gm", "Obs. global-mean"),
    # "seasonality",
    # "lat-grad-eof",
    # "lat-grad-pc",
    # "gm-ext",
    # "lat-grad-pc-emms",
    # "lat-grad-pc-ext",
    # "flying-carpet",
    # "yearly",
    # "monthly",
)
"""The figure's panels, in the order in which they are labelled

The labels follow the order of the steps in the method,
which is not the order in which the panels are laid out.
"""

MOSAIC = [
    ["timeseries", "timeseries", "interpolated-most", "gm"],
    ["timeseries", "timeseries", "interpolated-least", "seasonality"],
    ["locations", "counts", "lat-grad-eof", "lat-grad-pc"],
    ["gm-ext", "gm-ext", "lat-grad-pc-emms", "lat-grad-pc-ext"],
    ["flying-carpet", "yearly", "yearly", "monthly"],
]
"""Layout of the figure's panels

Four columns, because the columns have to be wide enough
to carry a panel's labels and its colour bar.
With more (hence narrower) columns, the layout engine runs out of room
and gives up, which is what leaves the panels sitting
in matplotlib's default positions rather than in the positions we asked for.
If a panel needs a width that four columns cannot express,
it is cheaper to move panels between rows than to split the columns further.

There is no point setting height ratios here:
the rows which hold a map are shrunk onto the height their map needs
while the figure is being laid out,
see [fit_rows_to_fixed_aspect_panels][].
"""

MAP_PANELS = ("locations", "interpolated-most", "interpolated-least")
"""The panels which show a map

These are the panels with a fixed aspect ratio,
which is what makes the figure's spacing fiddly.
There will not be any more of them.
"""

FIGURE_SIZE = (15.0, 11.2)
"""Size of the figure in inches

The width is set by the page, the height by what the panels need.
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


def label_name(ghg: str) -> str:
    """Get the label name for a given string"""
    replacements = {
        "n2o": "N2O",
    }

    res = ghg
    for old, new in replacements.items():
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


def get_n2o_all_data_with_bins(
    bundle_dir: Path = DEFAULT_BUNDLE_DIR,
    original_run_notebooks_dir: Path = DEFAULT_ORIGINAL_RUN_NOTEBOOKS_DIR,
    force_rerun: bool = False,
) -> pd.DataFrame:
    """
    Get the N2O observational network data, as it went into the binning

    This re-runs the original run's binning notebook if it needs to.
    That is slow the first time (the original run's environment has to be
    downloaded and installed), so the result is re-used if it is already there.

    Parameters
    ----------
    bundle_dir
        Directory in which to keep the original run's bundle

    original_run_notebooks_dir
        The original run's `notebooks-executed` directory

        Only used if we don't already have a copy of the notebook we need.

    force_rerun
        Re-run the notebook even if its output is already there

    Returns
    -------
        The observational network data, with the latitudinal and longitudinal
        bin of each observation added
    """
    out_file = bundle_dir / ALL_DATA_WITH_BINS_FILE
    if out_file.exists() and not force_rerun:
        logger.info(f"Using existing {out_file}")
        return pd.read_csv(out_file)

    base_notebook = (
        Path("calculate_n2o_monthly_fifteen_degree_pieces")
        / "only"
        / "1000_n2o_bin-observational-network.ipynb"
    )

    start_from = ensure_executed_notebook_available(
        base_notebook,
        original_run_notebooks_dir=original_run_notebooks_dir,
    )
    ensure_bundle_available(
        files_to_get=(
            "pyproject.toml",
            "pixi.lock",
            "v1.0.0-config-raw.yaml",
        ),
        files_to_get_tarred=(
            "src.tar.gz",
            "data--interim.tar.gz",
        ),
        bundle_dir=bundle_dir,
    )
    ensure_bundle_environment(bundle_dir)

    save_cell = f"""
# Added for the CMIP7 GHG manuscript.
# The original run never saved `all_data_with_bins` out,
# but we need it to show the observational network in the manuscript.
from pathlib import Path

manuscript_out_file = Path("{ALL_DATA_WITH_BINS_FILE.as_posix()}")
manuscript_out_file.parent.mkdir(exist_ok=True, parents=True)
all_data_with_bins.to_csv(manuscript_out_file, index=False)
manuscript_out_file
"""

    notebook_name = base_notebook.stem
    ipynb_to_run = bundle_dir / "notebooks-rerun" / f"{notebook_name}.ipynb"
    to_run = write_modified_notebook(
        start_from=start_from,
        out_py=MODIFIED_NOTEBOOKS_DIR / f"{notebook_name}.py",
        out_ipynb=ipynb_to_run,
        extra_cells=[save_cell],
        step_config_id="only",
    )
    run_notebook_from_bundle_dir(
        to_run,
        ipynb_to_run,
        bundle_dir=bundle_dir,
    )

    return pd.read_csv(out_file)


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


def add_network_group(indf: pd.DataFrame) -> pd.DataFrame:
    """
    Add the network group of each observation

    Parameters
    ----------
    indf
        Data to which to add the network group

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


def add_network_legend(ax: matplotlib.axes.Axes, **kwargs: object) -> None:
    """
    Add a legend which identifies the network groups

    Parameters
    ----------
    ax
        Axes to add the legend to

    **kwargs
        Passed on to `ax.legend`
    """
    ax.legend(
        fontsize="x-small",
        framealpha=0.9,
        handletextpad=0.2,
        columnspacing=1.0,
        **kwargs,
    )


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
    normalisation = matplotlib.colors.Normalize(vmin=-90.0, vmax=90.0)

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

    ax.indicate_inset_zoom(ax_inset, edgecolor="black")

    scatter = None
    for group, group_df in indf.groupby(NETWORK_GROUP_COLUMN):
        for axh in [ax, ax_inset]:
            scatter = axh.scatter(
                get_decimal_year(group_df),
                group_df["value"],
                c=group_df["latitude"],
                cmap=LATITUDE_COLOUR_MAP,
                norm=normalisation,
                marker=NETWORK_GROUP_MARKERS[group],
                s=10.0,
                alpha=0.4,
                linewidths=0.0,
            )

    if scatter is None:
        msg = "No data to plot"
        raise AssertionError(msg)

    units = indf["unit"].unique()
    if units.size != 1:
        msg = f"Expected exactly one unit, got {units}"
        raise AssertionError(msg)

    ax.set_ylabel(f"N$_2$O [{units[0]}]", fontsize="small")
    ax.set_xlabel("year", fontsize="small")
    ax.tick_params(labelsize="small")

    # Colour is doing latitude here, so the legend can only speak about shape
    add_network_legend(
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
    # On the panel label's line, rather than over the map (which would cover
    # stations) or below it (which would waste the height the map's fixed
    # aspect already costs us)
    add_network_legend(
        ax,
        loc="center left",
        bbox_to_anchor=(1.05, 0.5),
        ncols=1,
        # frameon=False,
    )


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

    return mesh


def plot_global_mean_from_obs_network(gm: xr.Dataset, ax: matplotlib.axes.Axes) -> None:
    """
    Plot global-mean derived from the observational network
    """
    gm_da = get_only_data_variable(gm)
    ax.scatter(
        gm_da["year"].values.squeeze(),
        gm_da.values.squeeze(),
    )
    ax.set_ylabel(f"[{gm_da.attrs['units']}]", fontsize="small")


def has_fixed_aspect(ax: matplotlib.axes.Axes) -> bool:
    """
    Determine whether an axes' aspect ratio is fixed

    A fixed aspect ratio means the axes' height follows from its width,
    so the axes cannot stretch to fill the row it is in.
    Our maps are the axes for which this is true.

    Parameters
    ----------
    ax
        Axes to check

    Returns
    -------
        `True` if `ax`'s aspect ratio is fixed, `False` otherwise
    """
    return ax.get_aspect() != "auto"


def fit_rows_to_fixed_aspect_panels(
    fig: matplotlib.figure.Figure,
    axes: dict[str, matplotlib.axes.Axes],
    n_iterations: int = 4,
) -> None:
    """
    Shrink each row of the figure onto the fixed-aspect panels it holds

    The maps have a fixed aspect ratio, so their height follows from their
    width and they cannot stretch to fill their row.
    Everything else in the row does stretch,
    which leaves a band of white space under each map
    and stops each map's latitude axis
    from lining up with its neighbours' latitude axes.

    So, we do the opposite: we shrink every row which holds a map
    down onto the height that map wants,
    and let the layout engine hand the height this frees up
    to the rows which can actually use it.
    Once a row is the height of its map,
    every panel in that row has the map's height too,
    which is what lines the latitude axes up.

    Row heights and panel widths depend on each other
    (through the space the layout engine sets aside for labels
    and colour bars), so we iterate rather than solving in one shot.
    A handful of iterations is plenty:
    a map's height is set by its width, which barely moves
    as the rows change height.

    Parameters
    ----------
    fig
        Figure to lay out

        This must be using a layout engine which respects height ratios,
        i.e. the constrained layout engine.

    axes
        The figure's panels

        Panels which span more than one row are ignored,
        because they say nothing about the height any single row needs.

    n_iterations
        Number of times to iterate
    """
    grid_spec = next(iter(axes.values())).get_subplotspec().get_gridspec()
    height_ratios = list(grid_spec.get_height_ratios())

    for _ in range(n_iterations):
        fig.draw_without_rendering()

        heights: dict[int, list[tuple[bool, float]]] = {}
        for ax in axes.values():
            row_span = ax.get_subplotspec().rowspan
            if row_span.stop - row_span.start != 1:
                continue

            heights.setdefault(row_span.start, []).append(
                (has_fixed_aspect(ax), ax.get_position().height)
            )

        for row, row_heights in heights.items():
            fixed_heights = [height for fixed, height in row_heights if fixed]
            if not fixed_heights:
                # Nothing in this row is holding the row's height back
                continue

            # The tallest map is the height the row needs,
            # the tallest panel is the height the row currently has.
            wanted = max(fixed_heights)
            have = max(height for _, height in row_heights)
            if have <= 0.0:
                # The layout engine has given up on this row,
                # so there is nothing to measure. Leave the row as it is:
                # squeezing it further would only make matters worse.
                continue

            height_ratios[row] *= wanted / have

        grid_spec.set_height_ratios(height_ratios)


def add_colour_bar(
    fig: matplotlib.figure.Figure,
    mappable: matplotlib.cm.ScalarMappable,
    ax: matplotlib.axes.Axes | list[matplotlib.axes.Axes],
    label: str,
    ticks: np.typing.ArrayLike | None = None,
    **kwargs: object,
) -> matplotlib.colorbar.Colorbar:
    """
    Add a colour bar to the figure

    Note
    ----
    Every colour bar we add costs the figure space, so each one needs a look
    before it is trusted, and there are two things to look at.

    First, a colour bar takes its width from the panel it belongs to,
    which makes that panel narrower.
    If that panel is a map, the map gets shorter too
    (its aspect ratio is fixed), and its whole row shrinks with it,
    because [fit_rows_to_fixed_aspect_panels][] fits the row to the map.
    A colour bar on a map is therefore much more expensive
    than a colour bar on any other panel.
    If the figure ends up carrying more colour bars than it can,
    the first thing to try is one colour bar shared between the panels
    which share a scale (pass a list of panels as `ax`),
    e.g. one for the two interpolated-coverage maps.

    Second, the colour bar has to be listed in `colour_bars`
    in [generate_n2o_methods_figure][], with the panel it belongs to,
    so that [tuck_colour_bars_against_their_panels][] moves it
    up against its panel.

    Parameters
    ----------
    fig
        Figure to add the colour bar to

    mappable
        Mappable to draw the colour bar for

    ax
        Panel (or panels) the colour bar belongs to

        The space for the colour bar is taken from these panels.

    label
        Label for the colour bar

    ticks
        Ticks to show on the colour bar

        If not supplied, the ticks are left to matplotlib.

    **kwargs
        Passed on to `fig.colorbar`

    Returns
    -------
        The colour bar which was added
    """
    colour_bar = fig.colorbar(mappable, ax=ax, pad=0.02, aspect=18, **kwargs)
    if ticks is not None:
        colour_bar.set_ticks(ticks)

    colour_bar.set_label(label, fontsize="small")
    colour_bar.ax.tick_params(labelsize="small")

    return colour_bar


def tuck_colour_bars_against_their_panels(
    fig: matplotlib.figure.Figure,
    colour_bars: list[tuple[matplotlib.colorbar.Colorbar, matplotlib.axes.Axes]],
    gap: float = 0.1,
) -> None:
    """
    Move each colour bar up against the panel it belongs to

    The layout engine parks a colour bar at the far side of the gap
    between its panel and the next one,
    where it reads as though it belongs to the panel on its right.
    Here we slide it back across the gap, up against its own panel.
    Nothing else moves: the space the colour bar leaves behind
    simply widens the gap before the next panel.

    Parameters
    ----------
    fig
        Figure being laid out

    colour_bars
        The figure's colour bars, each with the panel it belongs to

        Only colour bars which sit to the right of their panel are moved.

    gap
        Gap to leave between each panel and its colour bar, in inches

    Notes
    -----
    This freezes the layout,
    because the layout engine would otherwise undo the move,
    so it has to be the last thing which touches the figure's layout.
    """
    fig.draw_without_rendering()
    fig.set_layout_engine("none")

    gap_in_figure_coords = gap / fig.get_size_inches()[0]
    for colour_bar, panel in colour_bars:
        if colour_bar.orientation != "vertical":
            # Sits under its panel rather than beside it, so it is already home
            continue

        position = colour_bar.ax.get_position()
        colour_bar.ax.set_position(
            (
                panel.get_position().x1 + gap_in_figure_coords,
                position.y0,
                position.width,
                position.height,
            )
        )


def generate_n2o_methods_figure(
    outfile: Path,
    bundle_dir: Path,
    original_run_notebooks_dir: Path = DEFAULT_ORIGINAL_RUN_NOTEBOOKS_DIR,
    force_rerun: bool = False,
) -> Path:
    """
    Generate the N2O methods figure

    Parameters
    ----------
    outfile
        File in which to write the figure

    bundle_dir
        Directory in which to keep the original run's bundle

    original_run_notebooks_dir
        The original run's `notebooks-executed` directory

    force_rerun
        Re-run the original run's notebook even if its output is already there

    Returns
    -------
        `outfile`
    """
    fig, axes = plt.subplot_mosaic(
        MOSAIC,
        figsize=FIGURE_SIZE,
        per_subplot_kw={
            panel: {"projection": ccrs.PlateCarree()} for panel in MAP_PANELS
        },
        layout="constrained",
    )

    all_data_with_bins = add_network_group(
        get_n2o_all_data_with_bins(
            bundle_dir=bundle_dir,
            original_run_notebooks_dir=original_run_notebooks_dir,
            force_rerun=force_rerun,
        )
    )

    timeseries_scatter = plot_station_timeseries(all_data_with_bins, axes["timeseries"])
    plot_station_locations(all_data_with_bins, axes["locations"])
    counts_mesh = plot_observation_counts(all_data_with_bins, axes["counts"])

    latitude_colour_bar = add_colour_bar(
        fig,
        timeseries_scatter,
        ax=axes["timeseries"],
        label=r"latitude [$^{\circ}$N]",
        ticks=LAT_BIN_BOUNDS[::2],
    )
    # The points are drawn see-through so they don't hide each other,
    # but the colour bar should show the colours at full strength
    latitude_colour_bar.solids.set_alpha(1.0)

    counts_colour_bar = add_colour_bar(
        fig,
        counts_mesh,
        ax=axes["counts"],
        label="Number of input data points",
        ticks=np.arange(1, int(counts_mesh.norm.vmax) + 1),
    )

    # Both time axes cover the same period, even though the panels differ in width
    x_limits = (
        get_decimal_year(all_data_with_bins).min() - 1.0,
        get_decimal_year(all_data_with_bins).max() + 1.0,
    )
    for panel in ("timeseries", "counts"):
        axes[panel].set_xlim(x_limits)

    axes["counts"].set_xlabel("year", fontsize="small")

    interpolated_obs_file = (
        bundle_dir / "data/interim/n2o/n2o_observational-network_interpolated.nc"
    )
    interpolated_obs = xr.load_dataset(interpolated_obs_file)
    most_least_coverage = get_interpolated_input_coverage_info(
        all_data_with_bins, interpolated_obs
    )
    coverage_colour_bars_axes = []
    for key in ["most", "least"]:
        ax = axes[f"interpolated-{key}"]
        coverage_mesh = plot_coverage_and_interpolated(
            input_data=all_data_with_bins,
            interpolated=interpolated_obs,
            year_month=most_least_coverage[key],
            ax=ax,
        )
        colour_bar = add_colour_bar(
            fig,
            coverage_mesh,
            ax=ax,
            label=f"[{unit(all_data_with_bins)}]",
            # ticks=LAT_BIN_BOUNDS[::2],
        )

        coverage_colour_bars_axes.append([colour_bar, ax])

    global_mean_from_obs_network = xr.load_dataset(
        bundle_dir / "data/interim/n2o/n2o_observational-network_global-annual-mean.nc"
    )
    plot_global_mean_from_obs_network(global_mean_from_obs_network, axes["gm"])

    for label, (panel, title) in zip(string.ascii_lowercase, PANELS):
        axes[panel].set_title(
            f"$\\bf{{({label})}}$ {title}",
            loc="left",
            fontsize="medium",
            # fontweight="bold",
        )

    # The figure's colour bars, each with the panel it belongs to.
    # Any colour bar we add has to be listed here too,
    # otherwise it is left stranded next to its neighbour's panel.
    colour_bars = [
        (latitude_colour_bar, axes["timeseries"]),
        (counts_colour_bar, axes["counts"]),
        *((cb, ax) for cb, ax in coverage_colour_bars_axes),
    ]

    # These two are last, and in this order,
    # because they need to know how much space everything else has taken up
    # and the second of them freezes the layout.
    fit_rows_to_fixed_aspect_panels(fig, axes)
    tuck_colour_bars_against_their_panels(fig, colour_bars)

    outfile.parent.mkdir(exist_ok=True, parents=True)
    logger.info(f"Writing {outfile}")
    fig.savefig(outfile)
    plt.close(fig)

    return outfile
