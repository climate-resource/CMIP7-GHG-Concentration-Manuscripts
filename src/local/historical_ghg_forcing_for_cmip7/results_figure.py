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

from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path

import matplotlib.axes
import matplotlib.cm
import matplotlib.lines
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np
import pandas as pd
import xarray as xr
from loguru import logger

from local.cmip_ghg_generation import DEFAULT_ORIGINAL_RUN_NOTEBOOKS_DIR
from local.historical_ghg_forcing_for_cmip7.comparison_data import (
    LATITUDE_COLUMN,
    TIME_COLUMN,
    VALUE_COLUMN,
    ComparisonTimeseries,
    get_ch4_ice_core_comparisons,
    get_cmip6_comparisons,
    get_radiative_effect_per_unit,
)
from local.historical_ghg_forcing_for_cmip7.layout import (
    LayoutSettings,
    Panel,
    Row,
    create_figure,
    get_panel_axes_names,
    label_panels,
    lay_out_figure,
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
    "global-mean-comparison": 3.0,
    "output": 4.0,
    "spatial-comparison": 5.0,
}
"""Order to draw things in, back to front"""

TITLES = {
    "flying-carpet": "Native resolution",
    "monthly": "Monthly spatial-means",
    "yearly": "Yearly global-mean and comparison data",
    "yearly-diff": f"Difference from CMIP6 ({CMIP7_LABEL} - CMIP6)",
}
"""Title of each panel"""

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
                colour_bar=True,
                colour_bar_height=0.8,
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
  then the monthly means, and a panel to hold the figure's legend,
  which serves every panel because they all draw the same things the same way.
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


LEGEND_GROUP_ORDER = (
    "output",
    "output-at-latitude",
    "global-mean-comparison",
    "spatial-comparison",
    "context",
)
"""Order of the groups in the legend

Our output first, because it is the figure's subject,
then what it is compared against, then the context behind it all.
"""


class LegendCollector:
    """
    Collects the handles for the figure's legend

    Every panel draws the same things the same way,
    so the figure has one legend, which each panel adds to as it goes.
    """

    def __init__(self) -> None:
        self._handles: dict[str, tuple[str, matplotlib.artist.Artist]] = {}

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
            Group the handle belongs to, one of [LEGEND_GROUP_ORDER][]
        """
        if group not in LEGEND_GROUP_ORDER:
            raise ValueError(group)

        self._handles.setdefault(label, (group, handle))

    @property
    def handles(self) -> dict[str, matplotlib.artist.Artist]:
        """
        Handles, by label, grouped as [LEGEND_GROUP_ORDER][] says

        Within a group, handles are in the order they were first added.
        """
        ordered = sorted(
            self._handles.items(),
            key=lambda kv: LEGEND_GROUP_ORDER.index(kv[1][0]),
        )

        return {label: handle for label, (_, handle) in ordered}


def plot_obs_network_context(
    obs_network: pd.DataFrame,
    axes: Sequence[matplotlib.axes.Axes],
    segments: Sequence[tuple[float, float]],
    legend: LegendCollector,
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
    """
    decimal_year = get_decimal_year(obs_network)
    for ax, segment in zip(axes, segments):
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
                alpha=CONTEXT_ALPHA,
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
            group="context",
        )


def plot_input_timeseries_context(
    inputs: pd.DataFrame,
    axes: Sequence[matplotlib.axes.Axes],
    segments: Sequence[tuple[float, float]],
    legend: LegendCollector,
    label: str,
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
    """
    for latitude, lat_df in inputs.groupby("lat"):
        lat_df_sorted = lat_df.sort_values("year")
        # Annual values, so placed at the middle of the year
        times = lat_df_sorted["year"].to_numpy(dtype=float) + 0.5
        for ax, segment in zip(axes, segments):
            mask = in_segment(times, segment)
            ax.plot(
                times[mask],
                lat_df_sorted["value"].to_numpy()[mask],
                color=latitude_colour(latitude),
                alpha=0.5,
                linewidth=1.5,
                zorder=ZORDERS["context"],
            )

    legend.add(
        label,
        matplotlib.lines.Line2D(
            [], [], color=LEGEND_MARKER_COLOUR, alpha=0.5, linewidth=1.5
        ),
        group="context",
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

        if comparison.is_spatial:
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
                group="spatial-comparison",
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
                linewidth=1.5,
                zorder=ZORDERS["global-mean-comparison"],
            )
            legend.add(comparison.label, handle, group="global-mean-comparison")

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
        group="output-at-latitude",
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

        if share_y and i > 0:
            ax.set_ylim(axes[0].get_ylim())
            ax.tick_params(labelleft=False)

    axes[0].set_ylabel(label_name(f"[{units}]"), fontsize="small")
    axes[len(axes) // 2].set_xlabel(xlabel, fontsize="small")


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
    output_at_comparison_latitudes: bool = True,
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
        linewidth=2.0,
        zorder=ZORDERS["output"],
    )
    legend.add(CMIP7_LABEL, handle, group="output")

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
    if obs_network is not None:
        plot_obs_network_context(obs_network, axes, segments, legend)

    if inputs is not None:
        plot_input_timeseries_context(inputs, axes, segments, legend, inputs_label)

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
            linewidth=1.5,
            zorder=ZORDERS["output"],
        )
        legend.add(
            CMIP7_LABEL if region == "Global" else f"{CMIP7_LABEL} {region.lower()}",
            handle,
            group="output",
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
        linewidth=1.5,
    )
    for ax in axes:
        ax.axhline(0.0, color="0.6", linewidth=0.8, zorder=0)

    # One scale for every piece, so a difference means the same everywhere
    # and the radiative effect axis on the right can speak for all of them
    set_segment_y_limits(axes[:1], [(-np.inf, np.inf)], [diff])
    set_up_segments(axes, segments, units=units, xlabel="year", share_y=True)
    axes[0].yaxis.set_major_formatter(get_compact_scalar_formatter())
    add_radiative_effect_axis(axes[-1], gas, units)


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
    if anything_coloured_by_latitude:
        latitude_colour_bar = add_colour_bar(
            fig,
            matplotlib.cm.ScalarMappable(
                norm=LATITUDE_NORMALISATION, cmap=LATITUDE_COLOUR_MAP
            ),
            cax=axes["yearly-colour-bar"],
            label=r"latitude [$^{\circ}$N]",
            ticks=LAT_BIN_BOUNDS[::2],
        )
        latitude_colour_bar.solids.set_alpha(1.0)
    else:
        # E.g. C8F18, which has neither an observational network nor inputs
        axes["yearly-colour-bar"].set_visible(False)

    plot_difference_from_cmip6(
        gas, gm_yearly, cmip6_yearly, panel_axes("yearly-diff"), yearly_segments
    )

    legend_ax = axes[LEGEND_PANEL]
    legend_ax.axis("off")
    add_compact_legend(
        legend_ax,
        fontsize="small",
        handles=list(legend.handles.values()),
        labels=list(legend.handles),
        loc="center left",
        frameon=False,
        ncols=1 if len(legend.handles) <= 8 else 2,  # noqa: PLR2004
    )

    label_panels(ROWS, axes, TITLES, unlabelled=(LEGEND_PANEL,))
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

    other: tuple[ComparisonTimeseries, ...] = ()
    if gas == "ch4":
        other = get_ch4_ice_core_comparisons()

    return (cmip6_yearly, *other), (*cmip6_monthly, *other), cmip6_yearly


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
    # Imported here because these modules are about the methods figures,
    # and all we want from them is their data loading
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
    from local.historical_ghg_forcing_for_cmip7.n2o_methods_figure import (
        get_n2o_all_data_with_bins,
    )
    from local.historical_ghg_forcing_for_cmip7.plotting import (
        add_network_group,
    )

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
