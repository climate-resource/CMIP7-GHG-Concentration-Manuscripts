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

ALL_DATA_WITH_BINS_FILE = Path("manuscript-outputs") / "n2o_all-data-with-bins.csv"
"""Where the re-run notebook saves the data we want

Relative to the bundle's root directory,
because that is the notebook's working directory.
"""

LAT_BIN_BOUNDS = np.arange(-90, 91, 15)
"""Bounds of the latitudinal bins used by the original run

This mirrors `local.binning.LAT_BIN_BOUNDS` in the original run.
"""

LON_BIN_BOUNDS = np.arange(-180, 181, 60)
"""Bounds of the longitudinal bins used by the original run

This mirrors `local.binning.LON_BIN_BOUNDS` in the original run.
"""

LAT_AXIS_LIMITS = (-95, 95)
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

LEGEND_MARKER_COLOUR = "0.35"
"""Colour to draw legend markers in where colour means something else

In the timeseries panel colour shows latitude, so the legend
can only speak about the marker shape.
"""


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
        loc="lower right",
        bbox_to_anchor=(1.0, 1.0),
        ncols=2,
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


def align_latitude_axes(
    fig: matplotlib.figure.Figure,
    map_ax: matplotlib.axes.Axes,
    counts_ax: matplotlib.axes.Axes,
    counts_colour_bar_ax: matplotlib.axes.Axes,
) -> None:
    """
    Line the counts panel's latitude axis up with the map's

    Once they line up, a latitude sits at the same height in both panels,
    so the eye can run straight from a station on the map
    to the data it contributes.

    The map's aspect is fixed, so it doesn't fill the height of its cell
    and we only know the height it does fill once the figure has been laid out.
    Hence we lay the figure out, then move the counts panel to match the map.

    Parameters
    ----------
    fig
        Figure being drawn

    map_ax
        Axes holding the map

    counts_ax
        Axes holding the observation counts

    counts_colour_bar_ax
        Axes holding the observation counts' colour bar

        This moves with the counts panel, otherwise it is left standing
        taller than the panel it belongs to.
    """
    fig.draw_without_rendering()
    map_position = map_ax.get_position()

    # Constrained layout would just undo the positions we set below
    fig.set_layout_engine("none")

    for ax in (counts_ax, counts_colour_bar_ax):
        position = ax.get_position()
        ax.set_position(
            (position.x0, map_position.y0, position.width, map_position.height)
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
    all_data_with_bins = add_network_group(
        get_n2o_all_data_with_bins(
            bundle_dir=bundle_dir,
            original_run_notebooks_dir=original_run_notebooks_dir,
            force_rerun=force_rerun,
        )
    )

    panels = ("timeseries", "locations", "counts")
    labels = [string.ascii_lowercase[i] for i in range(len(panels))]
    fig, axes = plt.subplot_mosaic(
        [
            ["timeseries", "timeseries"],
            ["locations", "counts"],
        ],
        figsize=(7.5, 5.6),
        # The map keeps a true 2:1 aspect, so the bottom row is sized to suit it,
        # and given the wider column, rather than leaving it adrift in white space.
        # With the legend above the map rather than below, the row only has to
        # hold the map itself, so the height it used to spare goes to the timeseries.
        height_ratios=[1.0, 0.6],
        width_ratios=[1.28, 1.0],
        per_subplot_kw={"locations": {"projection": ccrs.PlateCarree()}},
        layout="constrained",
    )

    timeseries_scatter = plot_station_timeseries(all_data_with_bins, axes["timeseries"])
    plot_station_locations(all_data_with_bins, axes["locations"])
    counts_mesh = plot_observation_counts(all_data_with_bins, axes["counts"])

    latitude_colour_bar = fig.colorbar(
        timeseries_scatter,
        ax=axes["timeseries"],
        pad=0.02,
        aspect=18,
    )
    latitude_colour_bar.set_ticks(LAT_BIN_BOUNDS[::2])
    latitude_colour_bar.set_label(r"latitude [$^{\circ}$N]", fontsize="small")
    latitude_colour_bar.ax.tick_params(labelsize="small")
    latitude_colour_bar.solids.set_alpha(1.0)

    counts_colour_bar = fig.colorbar(
        counts_mesh,
        ax=axes["counts"],
        pad=0.02,
        aspect=18,
    )
    counts_colour_bar.set_ticks(np.arange(1, int(counts_mesh.norm.vmax) + 1))
    counts_colour_bar.set_label("Number of input data points", fontsize="small")
    counts_colour_bar.ax.tick_params(labelsize="small")

    # Both time axes cover the same period, even though the panels differ in width
    x_limits = (
        get_decimal_year(all_data_with_bins).min() - 1.0,
        get_decimal_year(all_data_with_bins).max() + 1.0,
    )
    for panel in ("timeseries", "counts"):
        axes[panel].set_xlim(x_limits)

    axes["counts"].set_xlabel("year", fontsize="small")

    for label, panel in zip(labels, panels):
        axes[panel].set_title(
            f"({label})", loc="left", fontsize="medium", fontweight="bold"
        )

    # Last, because it freezes the layout
    align_latitude_axes(fig, axes["locations"], axes["counts"], counts_colour_bar.ax)

    outfile.parent.mkdir(exist_ok=True, parents=True)
    logger.info(f"Writing {outfile}")
    fig.savefig(outfile)
    plt.close(fig)

    return outfile
