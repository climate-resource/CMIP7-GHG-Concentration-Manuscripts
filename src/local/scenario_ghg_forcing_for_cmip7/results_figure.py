"""
Generation of the scenario results figures

One figure per gas, the same for every gas:

- one scenario at native resolution, joined onto the end of the history
- the same scenario's monthly global- and hemispheric-means,
  joined onto the end of the history
- every CMIP7 scenario's yearly global-mean
- the CMIP7 scenarios' yearly global-means next to CMIP6's,
  in groups of high, medium and low scenarios

Built on the historical manuscript's results figure
([local.historical_ghg_forcing_for_cmip7.results_figure][]),
so the two read the same way.
The yearly panels are shared with the emissions figure
([local.scenario_ghg_forcing_for_cmip7.emissions_figure][]),
so the emissions can be read against the concentrations.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path

import matplotlib.axes
import matplotlib.lines
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from loguru import logger

from local.historical_ghg_forcing_for_cmip7.comparison_data import (
    TIME_COLUMN,
    VALUE_COLUMN,
)
from local.historical_ghg_forcing_for_cmip7.layout import (
    LEGEND_FONT_SIZE,
    LayoutSettings,
    Panel,
    Row,
    create_figure,
    get_panel_axes_names,
    label_panels,
    lay_out_figure,
    manuscript_style,
)
from local.historical_ghg_forcing_for_cmip7.plotting import (
    REGION_COLOURS,
    add_colour_bar,
    label_name,
    plot_flying_carpet,
)
from local.historical_ghg_forcing_for_cmip7.results_figure import (
    plot_line_in_segments,
    set_segment_y_limits,
    set_up_segments,
)
from local.scenario_ghg_forcing_for_cmip7.esgf_data import (
    HEMISPHERES,
    convert_units,
    get_history_yearly_global_mean,
    get_monthly_native_resolution,
    get_scenario_yearly_global_mean,
    get_spatial_means,
)
from local.scenario_ghg_forcing_for_cmip7.scenarios import (
    CMIP6_HISTORY_LAST_YEAR,
    CMIP6_LINESTYLE,
    CMIP6_SCENARIOS,
    CMIP7_HISTORY_LAST_YEAR,
    CMIP7_HISTORY_SOURCE_ID,
    CMIP7_LINESTYLE,
    CMIP7_SCENARIOS,
    HISTORY_COLOUR,
    SCENARIO_GROUPS,
    SCENARIOS,
    Scenario,
    get_scenarios_in_group,
)

HIGHLIGHTED_SCENARIO = "l"
"""Scenario shown at native resolution and as monthly spatial-means"""

NATIVE_RESOLUTION_HISTORY_YEARS = 3
"""Number of years of history to show at native resolution

Shown before the scenario, so the reader can see that the two join smoothly.
"""

NATIVE_RESOLUTION_SCENARIO_YEARS = 7
"""Number of years of the scenario to show at native resolution"""

MONTHLY_TIME_AXIS = (2020.0, 2030.0)
"""Time axis of the monthly panel

Short enough to read the seasonal cycle,
long enough to see the scenario carry on from the history.
"""

BREAK_YEAR = 2120.0
"""Year at which the yearly panels' time axes are broken

A little after 2100 (where the extensions start),
so that the 2100 tick and its label sit clear of the break marks.
"""

YEARLY_SEGMENTS = ((2010.0, BREAK_YEAR), (BREAK_YEAR, 2500.0))
"""Pieces the time axis of every yearly panel is broken into, as (start, end)

Starts in 2010, so that CMIP6's history is shown too.
The CMIP7-only panel uses the same time axis, so it can be read against the others.
"""

YEARLY_TICKS = ((2010.0, 2040.0, 2070.0, 2100.0), (2300.0, 2500.0))
"""Labelled ticks of each piece of [YEARLY_SEGMENTS][]

The second piece is too narrow to label every century
(2200 runs into 2100 across the break), so the others get minor ticks,
see [EXTENSION_MINOR_TICKS][].
"""

EXTENSION_MINOR_TICKS = (2200.0, 2400.0)
"""Unlabelled ticks of the extensions' piece of the time axis"""

SEGMENT_SHARES = (0.6, 0.4)
"""Share of a panel's width each piece of its time axis gets

The second piece covers about four times as many years as the first,
but they are smoother and matter less to most readers.
"""

SCENARIO_LINE_WIDTH = 0.8
"""Width of the lines which show the scenarios"""

HISTORY_LINE_WIDTH = 1.0
"""Width of the lines which show the history"""

MONTHLY_LINE_WIDTH = 0.8
"""Width of the lines in the monthly panel"""

SCENARIO_START_COLOUR = "0.5"
"""Colour of the line which marks where the history stops and the scenario starts"""

LEGEND_PANELS = {"CMIP7": "legend-cmip7", "CMIP6": "legend-cmip6"}
"""Name of the panel which holds each CMIP era's legend

The CMIP7 legend goes under the CMIP7 panel,
the CMIP6 one under the comparison panels (the only place CMIP6 is drawn).
"""

TITLES = {
    "native-resolution": (
        f"Native resolution, {SCENARIOS[HIGHLIGHTED_SCENARIO].label}"
    ),
    "monthly": f"Monthly spatial-means, {SCENARIOS[HIGHLIGHTED_SCENARIO].label}",
    "cmip7": "CMIP7 global-, annual-means",
    "high": "High scenarios",
    "medium": "Medium scenarios",
    "low": "Low scenarios",
}
"""Title of each panel"""

YEARLY_ROW_HEIGHT = 0.85
"""Height of each of the CMIP6 comparison panels, in inches"""

YEARLY_ROWS = (
    Row(
        panels=(
            Panel("cmip7", broken=True, broken_split=SEGMENT_SHARES, row_span=3),
            Panel(SCENARIO_GROUPS[0], broken=True, broken_split=SEGMENT_SHARES),
        ),
        height=YEARLY_ROW_HEIGHT,
    ),
    *(
        Row(
            panels=(Panel(group, broken=True, broken_split=SEGMENT_SHARES),),
            height=YEARLY_ROW_HEIGHT,
        )
        for group in SCENARIO_GROUPS[1:]
    ),
    Row(
        panels=(Panel(LEGEND_PANELS["CMIP7"]), Panel(LEGEND_PANELS["CMIP6"])),
        height=0.45,
    ),
)
"""Layout of the yearly panels, top to bottom

- Every CMIP7 scenario on the left, spanning the rows
  which hold the comparisons with CMIP6 on the right, one per group of scenarios.
- The legends for the yearly panels, one under each column.
"""

ROWS = (
    Row(
        panels=(
            Panel(
                "native-resolution",
                aspect=1.0,
                projection="3d",
                colour_bar=True,
                colour_bar_height=0.6,
            ),
            Panel("monthly"),
        ),
        height=1.05,
    ),
    *YEARLY_ROWS,
)
"""Layout of the figure's panels, top to bottom

- One scenario at native resolution, beside its monthly spatial-means.
- The yearly panels, see [YEARLY_ROWS][].
"""

LAYOUT_SETTINGS = LayoutSettings(align_right=True)
"""Layout settings

The comparison panels are stacked on one time axis,
so they have to line up at both ends.
"""


def as_frame(da: xr.DataArray) -> pd.DataFrame:
    """
    Get a timeseries as a frame

    Parameters
    ----------
    da
        Data, with a `time` dimension (as decimal years) and nothing else

    Returns
    -------
    :
        Data with a [TIME_COLUMN][] and a [VALUE_COLUMN][]
    """
    return pd.DataFrame(
        {TIME_COLUMN: da["time"].values, VALUE_COLUMN: da.values}
    ).sort_values(TIME_COLUMN)


def get_history_frame(gas: str, cmip_era: str, units: str) -> pd.DataFrame:
    """
    Get the yearly global-mean history we show

    Parameters
    ----------
    gas
        Gas of interest

    cmip_era
        CMIP era of interest

    units
        Units to return the data in

    Returns
    -------
    :
        History up to the end of the last year we show of it
    """
    last_year = {
        "CMIP7": CMIP7_HISTORY_LAST_YEAR,
        "CMIP6": CMIP6_HISTORY_LAST_YEAR,
    }[cmip_era]
    da = convert_units(get_history_yearly_global_mean(gas, cmip_era), units)

    return as_frame(da.sel(time=da["time"] < last_year + 1))


def branch_off_history(
    scenario: Scenario, scenario_pdf: pd.DataFrame, history: pd.DataFrame
) -> pd.DataFrame:
    """
    Join a scenario onto its history

    The history's last point before the scenario starts
    is put in front of the scenario,
    so the scenario's line branches off the history's.
    For most scenarios, that is the history's last point.
    Some (e.g. VL-CF) start before the end of the history we show,
    so they branch off it part way along.

    Parameters
    ----------
    scenario
        Scenario of interest

    scenario_pdf
        The scenario's data, with a [TIME_COLUMN][] and a [VALUE_COLUMN][]

    history
        The scenario's history, in the same form

    Returns
    -------
    :
        The scenario, starting from the history's last point before it
    """
    branch_point = history[history[TIME_COLUMN] < scenario_pdf[TIME_COLUMN].min()]
    if branch_point.empty:
        msg = f"{scenario.name} starts before the history we have"
        raise AssertionError(msg)

    return pd.concat([branch_point.iloc[-1:], scenario_pdf], ignore_index=True)


def get_scenario_frame(
    gas: str, scenario: Scenario, history: pd.DataFrame, units: str
) -> pd.DataFrame:
    """
    Get a scenario's yearly global-mean, joined onto its history

    Parameters
    ----------
    gas
        Gas of interest

    scenario
        Scenario of interest

    history
        The scenario's history, as returned by [get_history_frame][],
        see [branch_off_history][]

    units
        Units to return the data in

    Returns
    -------
    :
        The scenario, starting from the history's last point before it
    """
    da = convert_units(get_scenario_yearly_global_mean(gas, scenario), units)

    return branch_off_history(scenario, as_frame(da), history)


def plot_yearly_panel(  # noqa: PLR0913
    axes: Sequence[matplotlib.axes.Axes],
    segments: Sequence[tuple[float, float]],
    ticks: Sequence[Sequence[float]],
    scenarios: Mapping[Scenario, pd.DataFrame],
    histories: dict[str, pd.DataFrame],
    units: str,
) -> dict[str, matplotlib.lines.Line2D]:
    """
    Plot scenarios' yearly values, with their histories

    Parameters
    ----------
    axes
        Axes of each piece of the time axis

    segments
        Each piece of the time axis, as (start, end)

    ticks
        Ticks of each piece of the time axis

    scenarios
        Scenarios to plot, with their data, already joined onto their histories
        (see [branch_off_history][])

    histories
        History of each CMIP era to plot, as returned by [get_history_frame][]

    units
        Units the data is in, as written on the vertical axis

    Returns
    -------
    :
        Legend handle for everything drawn, by legend label
    """
    handles = {}
    foreground = []
    for cmip_era, history in histories.items():
        linestyle = CMIP7_LINESTYLE if cmip_era == "CMIP7" else CMIP6_LINESTYLE
        handles[f"{cmip_era} history"] = plot_line_in_segments(
            history,
            axes,
            segments,
            color=HISTORY_COLOUR,
            linestyle=linestyle,
            linewidth=HISTORY_LINE_WIDTH,
            zorder=3.0,
        )
        foreground.append(history)

    for scenario, pdf in scenarios.items():
        handles[scenario.label] = plot_line_in_segments(
            pdf,
            axes,
            segments,
            color=scenario.colour,
            linestyle=scenario.linestyle,
            linewidth=SCENARIO_LINE_WIDTH,
            # CMIP7 on top of CMIP6
            zorder=2.0 if scenario.cmip_era == "CMIP7" else 1.0,
        )
        foreground.append(pdf)

    # The pieces share one scale, so the lines can be followed across the break
    in_view = [
        pdf[
            (pdf[TIME_COLUMN] >= segments[0][0]) & (pdf[TIME_COLUMN] <= segments[-1][1])
        ]
        for pdf in foreground
    ]
    set_segment_y_limits(axes[:1], [(-np.inf, np.inf)], in_view)
    set_up_segments(axes, segments, units=units, xlabel="", share_y=True)
    for ax, segment, piece_ticks in zip(axes, segments, ticks):
        ax.set_xticks(piece_ticks)
        ax.set_xticks(
            [t for t in EXTENSION_MINOR_TICKS if segment[0] < t < segment[1]],
            minor=True,
        )
        # Setting the ticks can move the limits
        ax.set_xlim(segment)
        ax.xaxis.set_major_formatter("{x:.0f}")

    return handles


def get_native_resolution_to_show(gas: str) -> tuple[xr.DataArray, xr.DataArray]:
    """
    Get the native resolution data we show: the end of the history, then the scenario

    Parameters
    ----------
    gas
        Gas of interest

    Returns
    -------
    :
        The end of the history and the start of [HIGHLIGHTED_SCENARIO][],
        as one timeseries,
        and the whole of the history and the scenario (to 2100), likewise
    """
    history = get_monthly_native_resolution(gas, (CMIP7_HISTORY_SOURCE_ID,))
    history = history.sel(time=history["time"] < CMIP7_HISTORY_LAST_YEAR + 1)
    scenario = get_monthly_native_resolution(
        gas, SCENARIOS[HIGHLIGHTED_SCENARIO].source_ids[:1]
    )
    scenario = convert_units(scenario, history.attrs["units"])
    joined = xr.concat([history, scenario], dim="time", combine_attrs="override")

    first_year = CMIP7_HISTORY_LAST_YEAR + 1 - NATIVE_RESOLUTION_HISTORY_YEARS
    last_year = CMIP7_HISTORY_LAST_YEAR + NATIVE_RESOLUTION_SCENARIO_YEARS
    to_show = joined.sel(
        time=(joined["time"] >= first_year) & (joined["time"] < last_year + 1)
    )

    return to_show, joined


def to_year_month(da: xr.DataArray) -> xr.Dataset:
    """
    Split a decimal-year time axis into year and month

    The shape `plot_flying_carpet` expects (the historical output's shape),
    see [local.historical_ghg_forcing_for_cmip7.plotting][].

    Parameters
    ----------
    da
        Monthly data, with `time` as decimal years at the middle of each month

    Returns
    -------
    :
        `da` with `year` and `month` dimensions in place of `time`
    """
    year = np.floor(da["time"].values).astype(int)
    month = np.round((da["time"].values - year) * 12.0 + 0.5).astype(int)
    res = (
        da.assign_coords(year=("time", year), month=("time", month))
        .set_index(time=("year", "month"))
        .unstack("time")
    )

    return res.to_dataset(name=da.name or "value")


def plot_monthly(joined: xr.DataArray, ax: matplotlib.axes.Axes) -> None:
    """
    Plot the monthly global- and hemispheric-means

    Parameters
    ----------
    joined
        The history and the scenario at native resolution, as one timeseries

    ax
        Axes to plot on
    """
    in_view = joined.sel(
        time=(joined["time"] >= MONTHLY_TIME_AXIS[0])
        & (joined["time"] <= MONTHLY_TIME_AXIS[1])
    )
    spatial_means = get_spatial_means(in_view)
    for region in ("Global", *HEMISPHERES):
        pdf = as_frame(spatial_means.sel(region=region))
        ax.plot(
            pdf[TIME_COLUMN],
            pdf[VALUE_COLUMN],
            color=REGION_COLOURS[region],
            linewidth=MONTHLY_LINE_WIDTH,
            label="Global-mean" if region == "Global" else region,
        )

    scenario_start = CMIP7_HISTORY_LAST_YEAR + 1.0
    ax.axvline(
        scenario_start,
        color=SCENARIO_START_COLOUR,
        linestyle=":",
        linewidth=0.6,
        zorder=0,
    )
    for text, x, ha in (
        ("History ", scenario_start, "right"),
        (f" {SCENARIOS[HIGHLIGHTED_SCENARIO].label}", scenario_start, "left"),
    ):
        ax.text(
            x,
            0.97,
            text,
            transform=ax.get_xaxis_transform(),
            ha=ha,
            va="top",
            fontsize=LEGEND_FONT_SIZE,
            color=SCENARIO_START_COLOUR,
        )

    ax.set_xlim(MONTHLY_TIME_AXIS)
    ax.set_xticks(np.arange(MONTHLY_TIME_AXIS[0], MONTHLY_TIME_AXIS[1] + 1, 2))
    ax.xaxis.set_major_formatter("{x:.0f}")
    ax.tick_params(labelsize="small")
    ax.set_ylabel(label_name(f"[{joined.attrs['units']}]"), fontsize="small")
    ax.legend(loc="best", frameon=False, fontsize=LEGEND_FONT_SIZE, handlelength=1.2)


def add_scenario_legend(
    ax: matplotlib.axes.Axes,
    cmip_era: str,
    scenarios: Sequence[Scenario],
    handles: dict[str, matplotlib.lines.Line2D],
    n_rows: int = 3,
) -> None:
    """
    Add the legend for one CMIP era's history and scenarios

    Parameters
    ----------
    ax
        Axes to put the legend in

    cmip_era
        CMIP era of interest

    scenarios
        The era's scenarios, in the order to list them

    handles
        Legend handle for everything drawn, by legend label

    n_rows
        Number of rows to spread the entries over
    """
    entries = [
        ("History", handles[f"{cmip_era} history"]),
        *((s.label, handles[s.label]) for s in scenarios),
    ]
    linestyle = "solid" if cmip_era == "CMIP7" else "dashed"

    ax.axis("off")
    legend = ax.legend(
        handles=[handle for _, handle in entries],
        labels=[label for label, _ in entries],
        loc="center",
        ncols=int(np.ceil(len(entries) / n_rows)),
        frameon=False,
        fontsize=LEGEND_FONT_SIZE,
        title=f"{cmip_era} ({linestyle})",
        title_fontproperties={"weight": "bold", "size": LEGEND_FONT_SIZE},
        handlelength=1.8,
        columnspacing=0.8,
        borderaxespad=0.0,
    )
    # The layout measures the legend against the axes' starting size,
    # which it is much wider than, and would then squeeze every other row to fit.
    # The legend's axes takes up its row, so the legend has room once laid out.
    legend.set_in_layout(False)


def get_axes_of_panel(
    rows: tuple[Row, ...], axes: Mapping[str, matplotlib.axes.Axes], name: str
) -> list[matplotlib.axes.Axes]:
    """
    Get the axes of a panel (more than one if its time axis is broken)

    Parameters
    ----------
    rows
        Rows of the figure, top to bottom

    axes
        The figure's axes

    name
        Name of the panel

    Returns
    -------
    :
        The panel's axes, left to right
    """
    (panel,) = (p for row in rows for p in row.panels if p.name == name)

    return [axes[n] for n in get_panel_axes_names(panel)]


def plot_yearly_panels(
    rows: tuple[Row, ...],
    axes: Mapping[str, matplotlib.axes.Axes],
    histories: dict[str, pd.DataFrame],
    scenarios: Mapping[Scenario, pd.DataFrame],
    units: str,
) -> None:
    """
    Plot the yearly panels and their legends (see [YEARLY_ROWS][])

    Parameters
    ----------
    rows
        Rows of the figure, top to bottom (must include [YEARLY_ROWS][])

    axes
        The figure's axes

    histories
        History of each CMIP era, as returned by [get_history_frame][]

    scenarios
        Every scenario in [CMIP7_SCENARIOS][] and [CMIP6_SCENARIOS][],
        with its data, already joined onto its history (see [branch_off_history][])

    units
        Units the data is in, as written on the vertical axis
    """
    handles = plot_yearly_panel(
        get_axes_of_panel(rows, axes, "cmip7"),
        YEARLY_SEGMENTS,
        YEARLY_TICKS,
        {s: scenarios[s] for s in CMIP7_SCENARIOS},
        {"CMIP7": histories["CMIP7"]},
        units,
    )
    for i, group in enumerate(SCENARIO_GROUPS):
        handles.update(
            plot_yearly_panel(
                get_axes_of_panel(rows, axes, group),
                YEARLY_SEGMENTS,
                YEARLY_TICKS,
                {
                    s: scenarios[s]
                    for s in (
                        *get_scenarios_in_group(group, CMIP7_SCENARIOS),
                        *get_scenarios_in_group(group, CMIP6_SCENARIOS),
                    )
                },
                histories,
                units,
            )
        )
        # The comparison panels share one time axis, only the bottom one labels it
        if i < len(SCENARIO_GROUPS) - 1:
            for ax in get_axes_of_panel(rows, axes, group):
                ax.tick_params(labelbottom=False)

    for cmip_era, era_scenarios in (
        ("CMIP7", CMIP7_SCENARIOS),
        ("CMIP6", CMIP6_SCENARIOS),
    ):
        add_scenario_legend(
            axes[LEGEND_PANELS[cmip_era]], cmip_era, era_scenarios, handles
        )


@manuscript_style
def generate_results_figure(gas: str, outfile: Path) -> Path:
    """
    Generate a scenario results figure

    Parameters
    ----------
    gas
        Gas of interest

    outfile
        File in which to write the figure

    Returns
    -------
    :
        `outfile`
    """
    fig, axes = create_figure(ROWS, settings=LAYOUT_SETTINGS)

    native_resolution, joined = get_native_resolution_to_show(gas)
    units = native_resolution.attrs["units"]
    mesh = plot_flying_carpet(
        to_year_month(native_resolution), axes["native-resolution"]
    )
    add_colour_bar(
        fig,
        mesh,
        cax=axes["native-resolution-colour-bar"],
        label=label_name(f"{gas} [{units}]"),
        label_on_top=True,
    )

    plot_monthly(joined, axes["monthly"])

    histories = {
        cmip_era: get_history_frame(gas, cmip_era, units)
        for cmip_era in ("CMIP7", "CMIP6")
    }
    scenarios = {
        s: get_scenario_frame(gas, s, histories[s.cmip_era], units)
        for s in (*CMIP7_SCENARIOS, *CMIP6_SCENARIOS)
    }
    plot_yearly_panels(ROWS, axes, histories, scenarios, units)

    label_panels(ROWS, axes, TITLES, unlabelled=tuple(LEGEND_PANELS.values()))

    # Last, because it needs to know how much room everything takes up
    lay_out_figure(fig, axes, ROWS, settings=LAYOUT_SETTINGS)

    outfile.parent.mkdir(exist_ok=True, parents=True)
    logger.info(f"Writing {outfile}")
    fig.savefig(outfile)
    plt.close(fig)

    return outfile
