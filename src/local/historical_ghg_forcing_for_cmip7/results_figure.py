"""
Generation of the results figures

One figure per gas, the same for every gas, showing what we produced
and how it compares to what else is out there:

- our output at its native resolution
- our monthly global- and hemispheric-means
- our yearly global-mean
- the difference between our yearly global-mean and CMIP6's

The observational network (or, for gases without one, the inputs)
is drawn behind everything else, faded, as context:
the methods figures already show it properly.
The comparison datasets (see
[local.historical_ghg_forcing_for_cmip7.comparison_data][])
are what these figures are about, so they are drawn over the top of it.
"""

from __future__ import annotations

import itertools
import string
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path

import matplotlib.axes
import matplotlib.lines
import matplotlib.patheffects
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np
import pandas as pd
import xarray as xr
from loguru import logger

from local.cmip_ghg_generation import DEFAULT_ORIGINAL_RUN_NOTEBOOKS_DIR
from local.historical_ghg_forcing_for_cmip7.c4f10_like_methods_figure import (
    C4F10_LIKE_GASES,
    DROSTE_LABEL,
    get_droste_data,
)
from local.historical_ghg_forcing_for_cmip7.cfc12_like_methods_figure import (
    CFC12_LIKE_GASES,
    get_cfc12_like_all_data_with_bins,
)
from local.historical_ghg_forcing_for_cmip7.ch4_methods_figure import (
    get_ch4_all_data_with_bins,
)
from local.historical_ghg_forcing_for_cmip7.co2_methods_figure import (
    get_co2_all_data_with_bins,
)
from local.historical_ghg_forcing_for_cmip7.comparison_data import (
    LATITUDE_COLUMN,
    NOAA_TRENDS_UNITS,
    TIME_COLUMN,
    VALUE_COLUMN,
    ComparisonTimeseries,
    get_ch4_ice_core_comparisons,
    get_cmip6_comparisons,
    get_noaa_comparisons,
    get_radiative_effect_per_unit,
    get_uci_ch4_comparison,
)
from local.historical_ghg_forcing_for_cmip7.layout import (
    LayoutSettings,
    Panel,
    Row,
    create_figure,
    get_panel_axes_names,
    lay_out_figure,
)
from local.historical_ghg_forcing_for_cmip7.n2o_methods_figure import (
    get_n2o_all_data_with_bins,
)
from local.historical_ghg_forcing_for_cmip7.plotting import (
    LAT_BIN_BOUNDS,
    LAT_BIN_CENTRES,
    LATITUDE_COLOUR_MAP,
    LATITUDE_NORMALISATION,
    LEGEND_MARKER_COLOUR,
    NETWORK_GROUP_COLUMN,
    NETWORK_GROUP_MARKERS,
    OKABE_ITO,
    REGION_COLOURS,
    add_colour_bar,
    add_compact_legend,
    add_network_group,
    clear_ticks_near_break,
    get_decimal_year,
    get_only_data_variable,
    label_name,
    latitude_colour,
    plot_flying_carpet,
)

YEARLY_SEGMENTS: tuple[tuple[float, float | None], ...] = (
    (1.0, 1750.0),
    (1750.0, 2000.0),
    (2000.0, None),
)
"""Pieces the yearly time axis is broken into, as (start, end)

An end of `None` means the end of our data.

Each piece gets its own vertical scale.
The record barely moves for most of its first 1750 years,
then rises through the industrial era,
then keeps rising in the years the observational network covers,
and on one scale the first of those is a flat line
while the last is a steep one.
"""

YEARLY_SEGMENT_SHARES = (0.3, 0.35, 0.35)
"""Share of the yearly panels' width each piece of the time axis gets"""

MONTHLY_SEGMENTS: tuple[tuple[float, float | None], ...] = (
    (2000.0, 2020.0),
    (2020.0, None),
)
"""Pieces the monthly time axis is broken into, as (start, end)

The last few years get a piece of their own,
so their seasonal cycle can be read.
"""

MONTHLY_SEGMENT_SHARES = (0.6, 0.4)
"""Share of the monthly panel's width each piece of the time axis gets"""

SEGMENT_GAP = 0.5
"""Gap between the pieces of a broken time axis, in inches

Each piece has its own vertical scale, hence its own tick labels,
and they need somewhere to go.
"""

CMIP7_LABEL = "CMIP7"
"""What we call our own output in the legend"""

CONTEXT_ALPHA = 0.2
"""How see-through to draw the observational network

It is context, not the subject, so it is drawn faint enough
that anything drawn over it stands out.
"""

RECENT_MONTHLY_CONTEXT_ALPHA = 0.6
"""How see-through to draw the observational network in the last monthly piece

That piece covers only a few years, so the individual measurements
are far enough apart to be read, which is worth drawing them solidly enough for.
"""

CMIP7_LINE_WIDTH = 3.0
"""Width of the lines which show our output

Twice the width of everything else's lines, so our output stands out.
"""

OTHER_LINE_WIDTH = 1.5
"""Width of every other line"""

SHOW_OUTPUT_AT_COMPARISON_LATITUDES = False
"""Whether to draw our output in the latitudinal bins the spatial comparisons fall in

This is the fairer comparison for a site record,
but with more than a couple of sites the lines are too hard to read
(and to explain), so it is off for now.
"""

CONTEXT_MARKER_SIZE = 6.0
"""Marker size to draw the observational network with"""

COMPARISON_MARKER_SIZE = 30.0
"""Marker size to draw spatial comparison datasets with

Much bigger than the observational network's markers, and outlined,
so they stand out from it even where they are the same colour.
"""

DIFFERENCE_COLOUR = OKABE_ITO["black"]
"""Colour to draw the difference from CMIP6 in"""

ZORDERS = {
    "context": 1.0,
    "output-at-latitude": 2.0,
    "output": 3.0,
    "global-mean-comparison": 4.0,
    "spatial-comparison": 5.0,
}
"""Order to draw things in, back to front

The global-mean comparisons (e.g. CMIP6) are drawn over our output:
they are thin dashed lines, so they can be seen on top of our thick ones,
but our thick lines would hide them completely where the two agree.
"""

TITLES = {
    "flying-carpet": "Native resolution",
    "monthly": "Monthly spatial-means",
    "yearly": "Global- annual-mean",
    "yearly-diff": f"Yearly {CMIP7_LABEL} - CMIP6",
}
"""Title of each panel

For a panel split into pieces which each have their own vertical scale,
the title goes on the first piece only.
"""

PANELS_LABELLED_BY_PIECE = ("monthly", "yearly")
"""Panels whose pieces each get their own letter

Each piece of these has its own vertical scale,
so each has to be something the text can point at on its own.
The difference panel's pieces share one vertical scale, so it is one panel.
"""

LEGEND_HEADERS = ("CMIP forcings", "Comparison data", "Input data", "Latitude")
"""Sub-headers of the legend, in the order they are listed"""

LEGEND_COLUMNS = (
    ("CMIP forcings", "Input data"),
    ("Comparison data", "Latitude"),
)
"""Which sub-headers go in which column of the legend

Balanced by eye: the CMIP forcings always have six entries,
the other groups vary by gas.
"""

LEGEND_LATITUDES = LAT_BIN_BOUNDS[::-2]
"""Latitudes to show in the legend's latitude entries, north first"""

LEGEND_PANEL = "legend"
"""Name of the panel which holds the figure's legend"""

ROWS = (
    Row(
        panels=(
            Panel(
                "flying-carpet",
                aspect=1.0,
                projection="3d",
                colour_bar=True,
                colour_bar_height=0.6,
            ),
            Panel(
                "monthly",
                width=2.2,
                broken=True,
                broken_split=MONTHLY_SEGMENT_SHARES,
                broken_gap=SEGMENT_GAP,
            ),
            Panel(LEGEND_PANEL, width=1.0),
        ),
        height=3.0,
    ),
    Row(
        panels=(
            Panel(
                "yearly",
                broken=True,
                broken_split=YEARLY_SEGMENT_SHARES,
                broken_gap=SEGMENT_GAP,
            ),
        ),
        height=3.2,
    ),
    Row(
        panels=(
            Panel(
                "yearly-diff",
                broken=True,
                broken_split=YEARLY_SEGMENT_SHARES,
                broken_gap=SEGMENT_GAP,
            ),
        ),
        height=1.6,
    ),
)
"""Layout of the figure's panels, top to bottom

- Our output: the native resolution in the top-left,
  then the monthly means, and a panel to hold the legend
  for the monthly and yearly panels, which draw the same things the same way.
  The legend also says which colour stands for which latitude,
  so there is no latitude colour bar.
- The yearly global-mean.
- The yearly global-mean's difference from CMIP6,
  directly under the yearly global-mean and on the same time axis.
"""

LAYOUT_SETTINGS = LayoutSettings(align_right=True)
"""Layout settings

The yearly panel and the difference panel share a time axis,
so they have to line up at both ends.
"""


def resolve_segments(
    segments: Sequence[tuple[float, float | None]], last_time: float
) -> tuple[tuple[float, float], ...]:
    """
    Fill in the end of any piece of a time axis which runs to the end of the data

    Parameters
    ----------
    segments
        Pieces of the time axis, as (start, end)

    last_time
        Last time in the data, as a decimal year

    Returns
    -------
    :
        `segments`, with any end of `None` replaced by `last_time`
    """
    return tuple((start, last_time if end is None else end) for start, end in segments)


def in_segment(
    times: np.typing.ArrayLike, segment: tuple[float, float]
) -> np.typing.NDArray[np.bool_]:
    """
    Get which times fall in a piece of a time axis

    Parameters
    ----------
    times
        Times, as decimal years

    segment
        Piece of the time axis, as (start, end)

    Returns
    -------
    :
        Whether each of `times` is in `segment` (both ends included)
    """
    times = np.asarray(times, dtype=float)

    return (times >= segment[0]) & (times <= segment[1])


def get_annual_mean_as_frame(da: xr.DataArray) -> pd.DataFrame:
    """
    Get annual-mean data as a frame, placed at the middle of each year

    Parameters
    ----------
    da
        Data, with a `year` dimension

    Returns
    -------
    :
        Data with a [TIME_COLUMN][] and a [VALUE_COLUMN][]
    """
    return pd.DataFrame(
        {
            TIME_COLUMN: da["year"].values + 0.5,
            VALUE_COLUMN: da.values,
        }
    )


def get_monthly_as_frame(da: xr.DataArray) -> pd.DataFrame:
    """
    Get monthly data as a frame, placed at the middle of each month

    Parameters
    ----------
    da
        Data, with a `year` and a `month` dimension (and nothing else)

    Returns
    -------
    :
        Data with a [TIME_COLUMN][] and a [VALUE_COLUMN][]
    """
    pdf = da.to_dataframe(name=VALUE_COLUMN).reset_index()
    pdf[TIME_COLUMN] = get_decimal_year(pdf)

    return pdf[[TIME_COLUMN, VALUE_COLUMN]].sort_values(TIME_COLUMN)


def get_lat_bin(latitude: float) -> float:
    """
    Get the centre of the latitudinal bin a latitude falls in

    Parameters
    ----------
    latitude
        Latitude of interest

    Returns
    -------
    :
        Centre of the bin `latitude` falls in
    """
    i = np.clip(
        np.searchsorted(LAT_BIN_BOUNDS, latitude) - 1, 0, len(LAT_BIN_CENTRES) - 1
    )

    return float(LAT_BIN_CENTRES[i])


class LegendCollector:
    """
    Collects the handles for the figure's legend

    The monthly and yearly panels draw the same things the same way,
    so they share one legend, which each panel adds to as it goes.
    """

    def __init__(self) -> None:
        self._handles: dict[str, dict[str, matplotlib.artist.Artist]] = {
            header: {} for header in LEGEND_HEADERS
        }

    def add(self, label: str, handle: matplotlib.artist.Artist, group: str) -> None:
        """
        Add a handle, unless there is already one with this label

        Parameters
        ----------
        label
            Label of the handle

        handle
            Handle to add

        group
            Sub-header to list the handle under, one of [LEGEND_HEADERS][]
        """
        if group not in LEGEND_HEADERS:
            raise ValueError(group)

        self._handles[group].setdefault(label, handle)

    def get_columns(
        self,
    ) -> list[list[tuple[str, matplotlib.artist.Artist, bool]]]:
        """
        Get the legend's entries, column by column

        Returns
        -------
        :
            Entries of each column, top to bottom,
            as (label, handle, whether the entry is a sub-header).
            Within a group, handles are in the order they were first added.
            Groups with no entries, and columns with no groups, are left out.
        """
        columns = []
        for headers in LEGEND_COLUMNS:
            column = []
            for header in headers:
                if not self._handles[header]:
                    continue

                column.append(
                    (header, matplotlib.lines.Line2D([], [], linestyle="none"), True)
                )
                column.extend(
                    (label, handle, False)
                    for label, handle in self._handles[header].items()
                )

            if column:
                columns.append(column)

        return columns


def add_legend_with_sub_headers(
    ax: matplotlib.axes.Axes,
    legend: LegendCollector,
    title: str,
    fontsize: str = "x-small",
) -> None:
    """
    Add a legend whose entries are grouped under sub-headers

    Matplotlib's legend has a title but no sub-headers,
    so the sub-headers are entries with an empty handle,
    written in bold and moved over to where the handles start.

    Parameters
    ----------
    ax
        Axes to add the legend to

    legend
        Collector holding the entries

    title
        Title of the whole legend

    fontsize
        Size to draw the legend's text at
    """
    columns = legend.get_columns()
    # Matplotlib fills a legend column by column,
    # so padding each column out to the same length
    # is what puts each group in the column we asked for
    n_rows = max(len(column) for column in columns)
    entries = []
    for column in columns:
        entries.extend(column)
        entries.extend(
            ("", matplotlib.lines.Line2D([], [], linestyle="none"), False)
            for _ in range(n_rows - len(column))
        )

    ax.axis("off")
    add_compact_legend(
        ax,
        fontsize=fontsize,
        handles=[handle for _, handle, _ in entries],
        labels=[label for label, _, _ in entries],
        loc="center left",
        frameon=False,
        ncols=len(columns),
        title=title,
        alignment="left",
    )
    mpl_legend = ax.get_legend()
    mpl_legend.get_title().set_fontsize("small")
    mpl_legend.get_title().set_fontweight("bold")

    header_labels = {label for label, _, is_header in entries if is_header}
    # Relies on how matplotlib packs a legend (one box per column,
    # holding one box per entry, each of which is the handle then the label).
    # This is private, but there is no public way to line a sub-header up
    # with the handles rather than with the labels.
    for column_box in mpl_legend._legend_handle_box.get_children():
        for entry_box in column_box.get_children():
            handle_box, text_box = entry_box.get_children()
            if text_box._text.get_text() in header_labels:
                handle_box.set_width(0.0)
                entry_box.sep = 0.0
                text_box._text.set_fontweight("bold")


def add_latitude_entries(
    legend: LegendCollector, latitudes: Iterable[float] = LEGEND_LATITUDES
) -> None:
    """
    Add entries which say which colour stands for which latitude

    Parameters
    ----------
    legend
        Collector for the figure's legend

    latitudes
        Latitudes to add an entry for, in the order to list them
    """
    for latitude in latitudes:
        legend.add(
            f"{latitude:.0f}" + r"$^{\circ}$N",
            matplotlib.lines.Line2D(
                [],
                [],
                linestyle="none",
                marker="o",
                markersize=5,
                color=latitude_colour(latitude),
            ),
            group="Latitude",
        )


def plot_obs_network_context(
    obs_network: pd.DataFrame,
    axes: Sequence[matplotlib.axes.Axes],
    segments: Sequence[tuple[float, float]],
    legend: LegendCollector,
    alphas: Sequence[float] | None = None,
) -> None:
    """
    Plot the observational network, faded, as context

    This is the same view of it as the methods figures' first panel,
    i.e. coloured by latitude with a marker per network,
    just drawn so that it sits behind everything else.

    Parameters
    ----------
    obs_network
        Observational network data

    axes
        Axes of each piece of the time axis

    segments
        Each piece of the time axis, as (start, end)

    legend
        Collector for the figure's legend

    alphas
        How see-through to draw the network in each piece of the time axis

        If `None`, [CONTEXT_ALPHA][] everywhere.
    """
    if alphas is None:
        alphas = [CONTEXT_ALPHA] * len(axes)

    decimal_year = get_decimal_year(obs_network)
    for ax, segment, alpha in zip(axes, segments, alphas):
        in_seg = obs_network[in_segment(decimal_year, segment)]
        for group, group_df in in_seg.groupby(NETWORK_GROUP_COLUMN):
            ax.scatter(
                get_decimal_year(group_df),
                group_df["value"],
                c=group_df["latitude"],
                cmap=LATITUDE_COLOUR_MAP,
                norm=LATITUDE_NORMALISATION,
                marker=NETWORK_GROUP_MARKERS[group],
                s=CONTEXT_MARKER_SIZE,
                alpha=alpha,
                linewidths=0.0,
                zorder=ZORDERS["context"],
                # Thousands of points which are only there for context
                # make for a very heavy vector graphic
                rasterized=True,
            )

    for group in sorted(obs_network[NETWORK_GROUP_COLUMN].unique()):
        legend.add(
            f"{group} obs. network",
            matplotlib.lines.Line2D(
                [],
                [],
                linestyle="none",
                marker=NETWORK_GROUP_MARKERS[group],
                color=LEGEND_MARKER_COLOUR,
                alpha=0.6,
                markersize=4,
            ),
            group="Input data",
        )


def plot_input_timeseries_context(  # noqa: PLR0913
    inputs: pd.DataFrame,
    axes: Sequence[matplotlib.axes.Axes],
    segments: Sequence[tuple[float, float]],
    legend: LegendCollector,
    label: str,
    alphas: Sequence[float] | None = None,
) -> None:
    """
    Plot the inputs of a gas without an observational network, faded, as context

    Parameters
    ----------
    inputs
        Inputs, with a `lat`, `year` and `value` column

    axes
        Axes of each piece of the time axis

    segments
        Each piece of the time axis, as (start, end)

    legend
        Collector for the figure's legend

    label
        Label to give the inputs in the legend

    alphas
        How see-through to draw the inputs in each piece of the time axis

        If `None`, 0.5 everywhere.
    """
    if alphas is None:
        alphas = [0.5] * len(axes)

    for latitude, lat_df in inputs.groupby("lat"):
        lat_df_sorted = lat_df.sort_values("year")
        # Annual values, so placed at the middle of the year
        times = lat_df_sorted["year"].to_numpy(dtype=float) + 0.5
        for ax, segment, alpha in zip(axes, segments, alphas):
            mask = in_segment(times, segment)
            ax.plot(
                times[mask],
                lat_df_sorted["value"].to_numpy()[mask],
                color=latitude_colour(latitude),
                alpha=alpha,
                linewidth=OTHER_LINE_WIDTH,
                zorder=ZORDERS["context"],
            )

    legend.add(
        label,
        matplotlib.lines.Line2D(
            [], [], color=LEGEND_MARKER_COLOUR, alpha=0.5, linewidth=OTHER_LINE_WIDTH
        ),
        group="Input data",
    )


def plot_line_in_segments(
    pdf: pd.DataFrame,
    axes: Sequence[matplotlib.axes.Axes],
    segments: Sequence[tuple[float, float]],
    **kwargs: object,
) -> matplotlib.lines.Line2D:
    """
    Plot a line on every piece of a broken time axis

    Each piece is only given the part of the line it shows,
    so that nothing outside the piece can affect it.

    Parameters
    ----------
    pdf
        Data, with a [TIME_COLUMN][] and a [VALUE_COLUMN][]

    axes
        Axes of each piece of the time axis

    segments
        Each piece of the time axis, as (start, end)

    **kwargs
        Passed on to `ax.plot`

    Returns
    -------
    :
        The line drawn on the last piece, to use as a legend handle
    """
    line = None
    for ax, segment in zip(axes, segments):
        mask = in_segment(pdf[TIME_COLUMN], segment)
        (line,) = ax.plot(
            pdf.loc[mask, TIME_COLUMN], pdf.loc[mask, VALUE_COLUMN], **kwargs
        )

    if line is None:
        raise AssertionError

    return line


def plot_comparisons(
    comparisons: Iterable[ComparisonTimeseries],
    axes: Sequence[matplotlib.axes.Axes],
    segments: Sequence[tuple[float, float]],
    units: str,
    legend: LegendCollector,
) -> list[pd.DataFrame]:
    """
    Plot comparison datasets

    Datasets with spatial information are drawn as large, outlined markers,
    coloured by latitude like the observational network behind them.
    Global-mean datasets are drawn as lines.
    Both are drawn over the top of the observational network,
    so they stand out from it.

    Parameters
    ----------
    comparisons
        Datasets to plot

    axes
        Axes of each piece of the time axis

    segments
        Each piece of the time axis, as (start, end)

    units
        Units to plot in

    legend
        Collector for the figure's legend

    Returns
    -------
    :
        The data which was plotted, in `units`,
        so it can be taken into account when setting limits
    """
    plotted = []
    for comparison in comparisons:
        comparison_units = comparison.to_units(units)
        pdf = comparison_units.data
        plotted.append(pdf)

        if comparison.is_spatial and comparison.spatial_as_line:
            latitudes = pdf[LATITUDE_COLUMN].unique()
            for latitude in latitudes:
                handle = plot_line_in_segments(
                    pdf[pdf[LATITUDE_COLUMN] == latitude],
                    axes,
                    segments,
                    color=latitude_colour(latitude),
                    linestyle=comparison.linestyle,
                    linewidth=OTHER_LINE_WIDTH,
                    zorder=ZORDERS["spatial-comparison"],
                    # Outlined, like the spatial comparisons' markers,
                    # so it stands out from the observational network
                    # even where it is the same colour
                    path_effects=[
                        matplotlib.patheffects.Stroke(
                            linewidth=OTHER_LINE_WIDTH + 1.2, foreground="k"
                        ),
                        matplotlib.patheffects.Normal(),
                    ],
                )

            latitude_label = (
                f" ({latitudes[0]:.1f}" + r"$^{\circ}$N)" if len(latitudes) == 1 else ""
            )
            legend.add(
                f"{comparison.label}{latitude_label}",
                handle,
                group=comparison.legend_group,
            )

        elif comparison.is_spatial:
            for ax, segment in zip(axes, segments):
                in_seg = pdf[in_segment(pdf[TIME_COLUMN], segment)]
                ax.scatter(
                    in_seg[TIME_COLUMN],
                    in_seg[VALUE_COLUMN],
                    c=in_seg[LATITUDE_COLUMN],
                    cmap=LATITUDE_COLOUR_MAP,
                    norm=LATITUDE_NORMALISATION,
                    marker=comparison.marker,
                    s=COMPARISON_MARKER_SIZE,
                    edgecolors="k",
                    linewidths=0.6,
                    zorder=ZORDERS["spatial-comparison"],
                )

            latitudes = pdf[LATITUDE_COLUMN].unique()
            latitude_label = (
                f" ({latitudes[0]:.1f}" + r"$^{\circ}$N)" if len(latitudes) == 1 else ""
            )
            legend.add(
                f"{comparison.label}{latitude_label}",
                matplotlib.lines.Line2D(
                    [],
                    [],
                    linestyle="none",
                    marker=comparison.marker,
                    markerfacecolor=(
                        latitude_colour(latitudes[0])
                        if len(latitudes) == 1
                        else LEGEND_MARKER_COLOUR
                    ),
                    markeredgecolor="k",
                    markeredgewidth=0.6,
                    markersize=6,
                ),
                group=comparison.legend_group,
            )

        else:
            colour = (
                comparison.colour
                if comparison.colour is not None
                else REGION_COLOURS[comparison.region]
            )
            handle = plot_line_in_segments(
                pdf,
                axes,
                segments,
                color=colour,
                linestyle=comparison.linestyle,
                linewidth=OTHER_LINE_WIDTH,
                zorder=ZORDERS["global-mean-comparison"],
            )
            legend.add(comparison.label, handle, group=comparison.legend_group)

    return plotted


def plot_output_at_comparison_latitudes(
    native_resolution: xr.DataArray,
    comparisons: Iterable[ComparisonTimeseries],
    axes: Sequence[matplotlib.axes.Axes],
    segments: Sequence[tuple[float, float]],
    legend: LegendCollector,
) -> list[pd.DataFrame]:
    """
    Plot our yearly output in the latitudinal bins the spatial comparisons fall in

    A measurement at a site is not a global-mean,
    so this is what it should be read against.

    Parameters
    ----------
    native_resolution
        Our output at its native resolution,
        with a `year`, `month` and `lat` dimension

    comparisons
        Comparison datasets

        Only those with spatial information are used.

    axes
        Axes of each piece of the time axis

    segments
        Each piece of the time axis, as (start, end)

    legend
        Collector for the figure's legend

    Returns
    -------
    :
        The data which was plotted,
        so it can be taken into account when setting limits
    """
    lat_bins = sorted(
        {
            get_lat_bin(latitude)
            for comparison in comparisons
            if comparison.is_spatial
            for latitude in comparison.data[LATITUDE_COLUMN].unique()
        }
    )
    if not lat_bins:
        return []

    annual_mean = native_resolution.mean("month")
    plotted = []
    for lat_bin in lat_bins:
        pdf = get_annual_mean_as_frame(annual_mean.sel(lat=lat_bin))
        plotted.append(pdf)
        plot_line_in_segments(
            pdf,
            axes,
            segments,
            color=latitude_colour(lat_bin),
            linewidth=1.0,
            zorder=ZORDERS["output-at-latitude"],
        )

    legend.add(
        f"{CMIP7_LABEL} at comparison latitudes",
        matplotlib.lines.Line2D([], [], color=LEGEND_MARKER_COLOUR, linewidth=1.0),
        group="CMIP forcings",
    )

    return plotted


def set_segment_y_limits(  # noqa: PLR0913
    axes: Sequence[matplotlib.axes.Axes],
    segments: Sequence[tuple[float, float]],
    foreground: Iterable[pd.DataFrame],
    context: pd.DataFrame | None = None,
    context_percentiles: tuple[float, float] = (1.0, 99.0),
    margin: float = 0.05,
) -> None:
    """
    Give each piece of a broken time axis its own vertical scale

    Parameters
    ----------
    axes
        Axes of each piece of the time axis

    segments
        Each piece of the time axis, as (start, end)

    foreground
        Everything drawn on the pieces which must be on scale,
        each with a [TIME_COLUMN][] and a [VALUE_COLUMN][]

    context
        Observational network data drawn behind everything else

        Its bulk is kept on scale, its outliers are not:
        it is context, so it should not set the scale.

    context_percentiles
        Percentiles of the context data, in each piece, to keep on scale

    margin
        Room to leave above and below the data, as a fraction of its range
    """
    foreground = list(foreground)
    for ax, segment in zip(axes, segments):
        values = [
            pdf.loc[in_segment(pdf[TIME_COLUMN], segment), VALUE_COLUMN].to_numpy(
                dtype=float
            )
            for pdf in foreground
        ]
        if context is not None:
            context_in_segment = context.loc[
                in_segment(get_decimal_year(context), segment), "value"
            ].to_numpy(dtype=float)
            if context_in_segment.size > 0:
                values.append(np.percentile(context_in_segment, context_percentiles))

        all_values = np.concatenate(values)
        all_values = all_values[np.isfinite(all_values)]
        if all_values.size < 1:
            continue

        lowest = all_values.min()
        highest = all_values.max()
        values_range = highest - lowest
        if values_range <= 0.0:
            values_range = max(abs(highest), 1.0)

        ax.set_ylim(lowest - margin * values_range, highest + margin * values_range)


def set_up_segments(  # noqa: PLR0913
    axes: Sequence[matplotlib.axes.Axes],
    segments: Sequence[tuple[float, float]],
    units: str,
    xlabel: str,
    share_y: bool,
    max_x_ticks: int = 5,
) -> None:
    """
    Set up the pieces of a broken time axis

    Parameters
    ----------
    axes
        Axes of each piece of the time axis

    segments
        Each piece of the time axis, as (start, end)

    units
        Units of the values plotted

    xlabel
        Label for the time axis

        Only put on the middle piece, because the pieces are one time axis.

    share_y
        Whether the pieces share one vertical scale

        If they do, only the first piece labels it.

    max_x_ticks
        The most tick labels to put on each piece's time axis
    """
    for i, (ax, segment) in enumerate(zip(axes, segments)):
        ax.set_xlim(segment)
        ax.xaxis.set_major_locator(
            matplotlib.ticker.MaxNLocator(nbins=max_x_ticks, integer=True)
        )
        ax.tick_params(labelsize="small")
        if i < len(axes) - 1:
            # The next piece starts where this one ends,
            # so its first tick label says the same thing
            clear_ticks_near_break(ax, at="right")

        if share_y:
            # The pieces are one axis with bits cut out of it,
            # so only its outer ends are drawn (the broken axis trick)
            if i > 0:
                ax.set_ylim(axes[0].get_ylim())
                ax.spines["left"].set_visible(False)
                ax.tick_params(left=False, labelleft=False)

            if i < len(axes) - 1:
                ax.spines["right"].set_visible(False)

    if share_y:
        add_break_marks(axes)

    axes[0].set_ylabel(label_name(f"[{units}]"), fontsize="small")
    axes[len(axes) // 2].set_xlabel(xlabel, fontsize="small")


def add_break_marks(axes: Sequence[matplotlib.axes.Axes], size: float = 8.0) -> None:
    """
    Mark the breaks between the pieces of a broken axis

    Drawn as markers rather than lines, so that every mark is the same size
    and slant however wide the piece it is on.

    Parameters
    ----------
    axes
        Axes of each piece of the axis, left to right

    size
        Size of the marks, in points
    """
    kwargs = dict(
        marker=[(-1.0, -1.0), (1.0, 1.0)],
        markersize=size,
        linestyle="none",
        color="k",
        markeredgewidth=0.8,
        clip_on=False,
    )
    for left, right in itertools.pairwise(axes):
        left.plot([1.0, 1.0], [0.0, 1.0], transform=left.transAxes, **kwargs)
        right.plot([0.0, 0.0], [0.0, 1.0], transform=right.transAxes, **kwargs)


def get_compact_scalar_formatter() -> matplotlib.ticker.ScalarFormatter:
    """
    Get a formatter which puts very small or very large values' powers of ten aside

    The minor gases' values, and their radiative effects even more so, are tiny,
    so without an offset their tick labels are mostly zeros,
    and wider than the panel can spare.

    Returns
    -------
    :
        Formatter which writes the power of ten once, at the end of the axis
    """
    formatter = matplotlib.ticker.ScalarFormatter(useMathText=True)
    formatter.set_powerlimits((-2, 3))

    return formatter


def add_radiative_effect_axis(
    ax: matplotlib.axes.Axes, gas: str, units: str
) -> matplotlib.axes.Axes:
    """
    Add an axis on the right which shows the approximate radiative effect

    Parameters
    ----------
    ax
        Axes to add the axis to

    gas
        Gas of interest

    units
        Units of the values plotted on `ax`

    Returns
    -------
    :
        The secondary axis
    """
    per_unit = get_radiative_effect_per_unit(gas, units)
    secondary = ax.secondary_yaxis(
        "right",
        functions=(lambda x: x * per_unit, lambda x: x / per_unit),
    )
    secondary.yaxis.set_major_formatter(get_compact_scalar_formatter())
    secondary.set_ylabel(r"approx. radiative effect [W / m$^2$]", fontsize="small")
    secondary.tick_params(labelsize="small")

    return secondary


def plot_yearly(  # noqa: PLR0913
    gm_yearly: xr.DataArray,
    native_resolution: xr.DataArray,
    axes: Sequence[matplotlib.axes.Axes],
    segments: Sequence[tuple[float, float]],
    comparisons: Sequence[ComparisonTimeseries],
    legend: LegendCollector,
    obs_network: pd.DataFrame | None = None,
    inputs: pd.DataFrame | None = None,
    inputs_label: str = "Inputs",
    output_at_comparison_latitudes: bool = SHOW_OUTPUT_AT_COMPARISON_LATITUDES,
) -> None:
    """
    Plot our yearly global-mean against the comparison datasets

    Parameters
    ----------
    gm_yearly
        Our yearly global-mean

    native_resolution
        Our output at its native resolution

    axes
        Axes of each piece of the time axis

    segments
        Each piece of the time axis, as (start, end)

    comparisons
        Datasets to compare against

    legend
        Collector for the figure's legend

    obs_network
        Observational network data, to draw as context

    inputs
        Inputs of a gas without an observational network, to draw as context

    inputs_label
        Label for `inputs` in the legend

    output_at_comparison_latitudes
        Whether to also draw our output in the latitudinal bins
        of the spatial comparison datasets
    """
    units = gm_yearly.attrs["units"]
    if obs_network is not None:
        plot_obs_network_context(obs_network, axes, segments, legend)

    if inputs is not None:
        plot_input_timeseries_context(inputs, axes, segments, legend, inputs_label)

    gm_pdf = get_annual_mean_as_frame(gm_yearly)
    handle = plot_line_in_segments(
        gm_pdf,
        axes,
        segments,
        color=REGION_COLOURS["Global"],
        linewidth=CMIP7_LINE_WIDTH,
        zorder=ZORDERS["output"],
    )
    legend.add(f"{CMIP7_LABEL} global-mean", handle, group="CMIP forcings")

    foreground = [gm_pdf]
    foreground.extend(plot_comparisons(comparisons, axes, segments, units, legend))
    if output_at_comparison_latitudes:
        foreground.extend(
            plot_output_at_comparison_latitudes(
                native_resolution, comparisons, axes, segments, legend
            )
        )

    set_segment_y_limits(axes, segments, foreground, context=obs_network)
    set_up_segments(axes, segments, units=units, xlabel="year", share_y=False)


def plot_monthly(  # noqa: PLR0913
    gm_monthly: xr.DataArray,
    hm_monthly: xr.DataArray,
    axes: Sequence[matplotlib.axes.Axes],
    segments: Sequence[tuple[float, float]],
    comparisons: Sequence[ComparisonTimeseries],
    legend: LegendCollector,
    obs_network: pd.DataFrame | None = None,
    inputs: pd.DataFrame | None = None,
    inputs_label: str = "Inputs",
) -> None:
    """
    Plot our monthly global- and hemispheric-means against the comparison datasets

    Parameters
    ----------
    gm_monthly
        Our monthly global-mean

    hm_monthly
        Our monthly hemispheric-means, with a `lat` dimension
        whose values are the regions' names

    axes
        Axes of each piece of the time axis

    segments
        Each piece of the time axis, as (start, end)

    comparisons
        Datasets to compare against

    legend
        Collector for the figure's legend

    obs_network
        Observational network data, to draw as context

    inputs
        Inputs of a gas without an observational network, to draw as context

    inputs_label
        Label for `inputs` in the legend
    """
    units = gm_monthly.attrs["units"]
    # The last piece is short enough that the individual inputs can be read,
    # so they are drawn more solidly there
    alphas = [CONTEXT_ALPHA] * (len(axes) - 1) + [RECENT_MONTHLY_CONTEXT_ALPHA]
    if obs_network is not None:
        plot_obs_network_context(obs_network, axes, segments, legend, alphas=alphas)

    if inputs is not None:
        plot_input_timeseries_context(
            inputs,
            axes,
            segments,
            legend,
            inputs_label,
            alphas=[0.5] * (len(axes) - 1) + [RECENT_MONTHLY_CONTEXT_ALPHA],
        )

    foreground = []
    # In the same order as the regions are named everywhere else
    for region, da in (
        ("Global", gm_monthly),
        *(
            (region, hm_monthly.sel(lat=region))
            for region in ("Northern hemisphere", "Southern hemisphere")
        ),
    ):
        pdf = get_monthly_as_frame(da)
        foreground.append(pdf)
        handle = plot_line_in_segments(
            pdf,
            axes,
            segments,
            color=REGION_COLOURS[region],
            linewidth=CMIP7_LINE_WIDTH,
            zorder=ZORDERS["output"],
        )
        legend.add(
            f"{CMIP7_LABEL} {'global-mean' if region == 'Global' else region.lower()}",
            handle,
            group="CMIP forcings",
        )

    foreground.extend(plot_comparisons(comparisons, axes, segments, units, legend))

    set_segment_y_limits(axes, segments, foreground, context=obs_network)
    set_up_segments(axes, segments, units=units, xlabel="year", share_y=False)


def plot_difference_from_cmip6(
    gas: str,
    gm_yearly: xr.DataArray,
    cmip6_yearly: ComparisonTimeseries,
    axes: Sequence[matplotlib.axes.Axes],
    segments: Sequence[tuple[float, float]],
) -> None:
    """
    Plot the difference between our yearly global-mean and CMIP6's

    Parameters
    ----------
    gas
        Gas of interest

    gm_yearly
        Our yearly global-mean

    cmip6_yearly
        CMIP6's yearly global-mean

    axes
        Axes of each piece of the time axis

    segments
        Each piece of the time axis, as (start, end)
    """
    units = gm_yearly.attrs["units"]
    ours = get_annual_mean_as_frame(gm_yearly).set_index(TIME_COLUMN)[VALUE_COLUMN]
    cmip6 = cmip6_yearly.to_units(units).data.set_index(TIME_COLUMN)[VALUE_COLUMN]
    # Both are placed at the middle of each year, so the times line up exactly
    diff = (ours - cmip6).dropna().rename(VALUE_COLUMN).reset_index()

    plot_line_in_segments(
        diff,
        axes,
        segments,
        color=DIFFERENCE_COLOUR,
        linewidth=OTHER_LINE_WIDTH,
    )
    for ax in axes:
        ax.axhline(0.0, color="0.6", linewidth=0.8, zorder=0)

    # One scale for every piece, so a difference means the same everywhere
    # and the radiative effect axis on the right can speak for all of them
    set_segment_y_limits(axes[:1], [(-np.inf, np.inf)], [diff])
    set_up_segments(axes, segments, units=units, xlabel="year", share_y=True)
    axes[0].yaxis.set_major_formatter(get_compact_scalar_formatter())
    add_radiative_effect_axis(axes[-1], gas, units)


def label_results_panels(
    axes: Mapping[str, matplotlib.axes.Axes],
) -> dict[str, list[str]]:
    """
    Give each panel its title, labelled in reading order

    Unlike [local.historical_ghg_forcing_for_cmip7.layout.label_panels][],
    each piece of the panels in [PANELS_LABELLED_BY_PIECE][]
    gets its own letter.

    Parameters
    ----------
    axes
        The figure's axes

    Returns
    -------
    :
        Labels given to each panel's pieces, e.g. `{"monthly": ["(b)", "(c)"]}`
    """
    letters = iter(string.ascii_lowercase)
    res = {}
    for row in ROWS:
        for panel in row.panels:
            if panel.name == LEGEND_PANEL:
                continue

            names = get_panel_axes_names(panel)
            if panel.name not in PANELS_LABELLED_BY_PIECE:
                names = names[:1]

            res[panel.name] = []
            for i, name in enumerate(names):
                label = f"({next(letters)})"
                res[panel.name].append(label)
                title = f"$\\bf{{{label}}}$"
                if i < 1:
                    title = f"{title} {TITLES[panel.name]}"

                axes[name].set_title(title, loc="left", fontsize="medium")

    return res


def load_output(
    gas: str, bundle_dir: Path
) -> tuple[xr.DataArray, xr.DataArray, xr.DataArray, xr.DataArray]:
    """
    Load our output for a gas

    Parameters
    ----------
    gas
        Gas of interest

    bundle_dir
        Directory which holds the original run's bundle

    Returns
    -------
    :
        Our native resolution output, yearly global-mean,
        monthly global-mean and monthly hemispheric-means
        (the latter with the regions' names as its `lat` values)
    """
    gas_dir = bundle_dir / "data" / "interim" / gas

    native_resolution = get_only_data_variable(
        xr.load_dataset(gas_dir / f"{gas}_fifteen-degree_monthly.nc")
    )
    gm_yearly = get_only_data_variable(
        xr.load_dataset(gas_dir / f"{gas}_global-mean_annual-mean.nc")
    )
    gm_monthly = get_only_data_variable(
        xr.load_dataset(gas_dir / f"{gas}_global-mean_monthly.nc")
    )
    hm_monthly = get_only_data_variable(
        xr.load_dataset(gas_dir / f"{gas}_hemispheric-mean_monthly.nc")
    )
    sh_lat = -45.0
    hm_monthly = hm_monthly.assign_coords(
        lat=[
            "Southern hemisphere" if v == sh_lat else "Northern hemisphere"
            for v in hm_monthly["lat"].values
        ]
    )

    return native_resolution, gm_yearly, gm_monthly, hm_monthly


def generate_results_figure(  # noqa: PLR0913
    gas: str,
    outfile: Path,
    bundle_dir: Path,
    yearly_comparisons: Sequence[ComparisonTimeseries],
    monthly_comparisons: Sequence[ComparisonTimeseries],
    cmip6_yearly: ComparisonTimeseries,
    obs_network: pd.DataFrame | None = None,
    inputs: pd.DataFrame | None = None,
    inputs_label: str = "Inputs",
    force_rerun: bool = False,
) -> Path:
    """
    Generate a results figure

    Parameters
    ----------
    gas
        Gas of interest

    outfile
        File in which to write the figure

    bundle_dir
        Directory which holds the original run's bundle

    yearly_comparisons
        Datasets to compare our yearly global-mean against

    monthly_comparisons
        Datasets to compare our monthly spatial-means against

    cmip6_yearly
        CMIP6's yearly global-mean, which the difference panel is taken against

    obs_network
        Observational network data, drawn as context

    inputs
        Inputs of a gas without an observational network, drawn as context

    inputs_label
        Label for `inputs` in the legend

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

    native_resolution, gm_yearly, gm_monthly, hm_monthly = load_output(gas, bundle_dir)
    last_time = float(gm_yearly["year"].max()) + 1.0
    yearly_segments = resolve_segments(YEARLY_SEGMENTS, last_time)
    monthly_segments = resolve_segments(MONTHLY_SEGMENTS, last_time)

    fig, axes = create_figure(ROWS, settings=LAYOUT_SETTINGS)

    def panel_axes(name: str) -> list[matplotlib.axes.Axes]:
        (panel,) = (p for row in ROWS for p in row.panels if p.name == name)
        return [axes[n] for n in get_panel_axes_names(panel)]

    legend = LegendCollector()

    max_year = int(native_resolution["year"].max())
    flying_carpet_mesh = plot_flying_carpet(
        native_resolution.sel(year=range(max_year - 9, max_year + 1)).to_dataset(),
        axes["flying-carpet"],
    )
    add_colour_bar(
        fig,
        flying_carpet_mesh,
        cax=axes["flying-carpet-colour-bar"],
        label=label_name(f"{gas} [{native_resolution.attrs['units']}]"),
        label_on_top=True,
    )

    plot_monthly(
        gm_monthly,
        hm_monthly,
        panel_axes("monthly"),
        monthly_segments,
        comparisons=monthly_comparisons,
        legend=legend,
        obs_network=obs_network,
        inputs=inputs,
        inputs_label=inputs_label,
    )

    plot_yearly(
        gm_yearly,
        native_resolution,
        panel_axes("yearly"),
        yearly_segments,
        comparisons=yearly_comparisons,
        legend=legend,
        obs_network=obs_network,
        inputs=inputs,
        inputs_label=inputs_label,
    )

    anything_coloured_by_latitude = (
        obs_network is not None
        or inputs is not None
        or any(c.is_spatial for c in (*yearly_comparisons, *monthly_comparisons))
    )
    # Not for e.g. C8F18, which has neither an observational network nor inputs
    if anything_coloured_by_latitude:
        add_latitude_entries(legend)

    plot_difference_from_cmip6(
        gas, gm_yearly, cmip6_yearly, panel_axes("yearly-diff"), yearly_segments
    )

    piece_labels = label_results_panels(axes)
    add_legend_with_sub_headers(
        axes[LEGEND_PANEL],
        legend,
        title=(
            f"Legend for panels {piece_labels['monthly'][0]} "
            f"to {piece_labels['yearly'][-1]}"
        ),
    )

    # Last, because it needs to know how much room everything takes up
    lay_out_figure(fig, axes, ROWS, settings=LAYOUT_SETTINGS)

    outfile.parent.mkdir(exist_ok=True, parents=True)
    logger.info(f"Writing {outfile}")
    fig.savefig(outfile)
    plt.close(fig)

    return outfile


def get_comparisons(
    gas: str,
) -> tuple[
    tuple[ComparisonTimeseries, ...],
    tuple[ComparisonTimeseries, ...],
    ComparisonTimeseries,
]:
    """
    Get the datasets to compare a gas' output against

    Parameters
    ----------
    gas
        Gas of interest

    Returns
    -------
    :
        Datasets to compare the yearly output against,
        datasets to compare the monthly output against,
        and CMIP6's yearly global-mean
    """
    (cmip6_yearly,) = get_cmip6_comparisons(gas, "yr")
    cmip6_monthly = get_cmip6_comparisons(
        gas,
        "mon",
        regions=("Global", "Northern hemisphere", "Southern hemisphere"),
    )

    yearly_other: list[ComparisonTimeseries] = []
    monthly_other: list[ComparisonTimeseries] = []
    if gas in NOAA_TRENDS_UNITS:
        # The yearly panel is about the trend, not the seasonal cycle,
        # so it gets the records with their seasonal cycle removed
        yearly_other.extend(get_noaa_comparisons(gas, deseasonalised=True))
        monthly_other.extend(get_noaa_comparisons(gas, deseasonalised=False))

    if gas == "ch4":
        ice_cores = get_ch4_ice_core_comparisons()
        yearly_other.extend(ice_cores)
        monthly_other.extend(ice_cores)
        yearly_other.append(get_uci_ch4_comparison(deseasonalised=True))
        monthly_other.append(get_uci_ch4_comparison(deseasonalised=False))

    return (
        (cmip6_yearly, *yearly_other),
        (*cmip6_monthly, *monthly_other),
        cmip6_yearly,
    )


def get_context(
    gas: str,
    bundle_dir: Path,
    original_run_notebooks_dir: Path = DEFAULT_ORIGINAL_RUN_NOTEBOOKS_DIR,
    force_rerun: bool = False,
) -> Mapping[str, object]:
    """
    Get the data to draw behind a gas' output, as context

    Parameters
    ----------
    gas
        Gas of interest

    bundle_dir
        Directory which holds the original run's bundle

    original_run_notebooks_dir
        The original run's `notebooks-executed` directory

        Only used if the observational network data has to be re-generated.

    force_rerun
        Re-generate the observational network data,
        even if it is already there

    Returns
    -------
    :
        Keyword arguments for [generate_results_figure][]
        which give it its context
    """
    obs_network_getters = {
        "co2": get_co2_all_data_with_bins,
        "ch4": get_ch4_all_data_with_bins,
        "n2o": get_n2o_all_data_with_bins,
    }
    if gas in obs_network_getters:
        return {
            "obs_network": add_network_group(
                obs_network_getters[gas](
                    bundle_dir=bundle_dir,
                    original_run_notebooks_dir=original_run_notebooks_dir,
                    force_rerun=force_rerun,
                )
            )
        }

    if gas in CFC12_LIKE_GASES:
        return {
            "obs_network": add_network_group(
                get_cfc12_like_all_data_with_bins(gas, bundle_dir)
            )
        }

    if gas in C4F10_LIKE_GASES:
        return {
            "inputs": get_droste_data(gas, bundle_dir),
            "inputs_label": DROSTE_LABEL,
        }

    # C8F18 is taken from CMIP6 whole, so there is nothing to draw behind it
    return {}


def generate_results_figure_for_gas(
    gas: str,
    outfile: Path,
    bundle_dir: Path,
    original_run_notebooks_dir: Path = DEFAULT_ORIGINAL_RUN_NOTEBOOKS_DIR,
    force_rerun: bool = False,
) -> Path:
    """
    Generate a gas' results figure, with the context and comparisons that gas has

    Parameters
    ----------
    gas
        Gas of interest

    outfile
        File in which to write the figure

    bundle_dir
        Directory which holds the original run's bundle

    original_run_notebooks_dir
        The original run's `notebooks-executed` directory

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

    yearly_comparisons, monthly_comparisons, cmip6_yearly = get_comparisons(gas)

    return generate_results_figure(
        gas,
        outfile,
        bundle_dir=bundle_dir,
        yearly_comparisons=yearly_comparisons,
        monthly_comparisons=monthly_comparisons,
        cmip6_yearly=cmip6_yearly,
        force_rerun=force_rerun,
        **get_context(
            gas,
            bundle_dir,
            original_run_notebooks_dir=original_run_notebooks_dir,
            force_rerun=False,
        ),
    )
