"""
Generation of the scenario manuscript's methods figures

The figures follow the plans in the manuscript's `methods.tex`:

- time axes run to 2500 and are broken a little after 2100
- panels which compare scenarios show all of them, including VL-CF
- the MAGICC ensemble panels (CO2 and the gases processed like N2O)
  show one scenario, L
  (the gases which use a single concentration projection
  are the same in every scenario anyway)
- latitudinal gradient (and seasonality change) PCs are shown over time,
  for all scenarios,
  the regressions behind them are in their own (appendix) figure
- what applies to only one panel is explained in that panel's own legend,
  what applies to more than one panel is explained in the figure's legend panel,
  in groups whose titles say which panels they apply to

Built on the historical manuscript's layout and broken axis helpers,
so all the manuscript's figures read the same way.
"""

from __future__ import annotations

import dataclasses
from collections.abc import Sequence
from pathlib import Path

import matplotlib.artist
import matplotlib.axes
import matplotlib.lines
import matplotlib.patches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

from local.historical_ghg_forcing_for_cmip7.comparison_data import (
    TIME_COLUMN,
    VALUE_COLUMN,
)
from local.historical_ghg_forcing_for_cmip7.layout import (
    LEGEND_FONT_SIZE,
    TITLE_FONT_SIZE,
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
    get_ticks_with_ends,
    label_name,
    plot_pc_timeseries_regression,
)
from local.historical_ghg_forcing_for_cmip7.results_figure import (
    in_segment,
    plot_line_in_segments,
    set_segment_y_limits,
    set_up_segments,
)
from local.scenario_ghg_forcing_for_cmip7.bundle import BUNDLE_DIR
from local.scenario_ghg_forcing_for_cmip7.magicc_output import (
    MAGICC_AR6_VERSION,
    MAGICC_GAS_NAMES,
    MAGICC_VERSION,
)
from local.scenario_ghg_forcing_for_cmip7.methods_data import (
    CO2_SEASONALITY_REGRESSION,
    LAST_HISTORICAL_YEAR,
    get_historical_pcs,
    get_history,
    get_magicc_concentrations,
    get_official_emissions,
    get_projected_co2_seasonality_change_pc,
    get_projected_pcs,
    get_regression,
    get_scenario,
    get_single_projection_inverse_emissions,
    get_units,
    get_western_et_al_2024,
    get_western_et_al_2024_emissions,
    get_wmo_2022,
)
from local.scenario_ghg_forcing_for_cmip7.scenarios import (
    CMIP7_SCENARIOS,
    HISTORY_COLOUR,
    SCENARIOS,
)

SHOWN_SCENARIO = "l"
"""Scenario the single-scenario panels show"""

ALL_SCENARIOS = CMIP7_SCENARIOS
"""Scenarios shown in the panels which compare scenarios, lowest to highest

Including VL-CF, which leaves history in 2015, rather than 2022.
"""

EXTENSIONS_START_YEAR = 2100.0
"""Year in which the extensions (beyond 2100) start"""

BREAK_YEAR = 2120.0
"""Year at which every two-piece time axis is broken

A little after where the extensions start ([EXTENSIONS_START_YEAR][]),
so that the 2100 tick and its label sit clear of the break marks.
"""

END_YEAR = 2500.0
"""Last year shown"""

LONG_SEGMENTS = ((1850.0, BREAK_YEAR), (BREAK_YEAR, END_YEAR))
"""Time axis pieces for panels which show the whole history behind a projection"""

SCENARIO_SEGMENTS = ((2010.0, BREAK_YEAR), (BREAK_YEAR, END_YEAR))
"""Time axis pieces for panels which show the projections, with recent history"""

ONE_BOX_SEGMENTS = ((2000.0, BREAK_YEAR), (BREAK_YEAR, END_YEAR))
"""Time axis pieces for the one-box vs. MAGICC panels"""

CH4_VERSIONS_SEGMENTS = ((2020.0, BREAK_YEAR), (BREAK_YEAR, END_YEAR))
"""Time axis pieces for the comparison of MAGICC versions for CH4"""

ENSEMBLE_SEGMENTS = ((1990.0, BREAK_YEAR), (BREAK_YEAR, END_YEAR))
"""Time axis pieces for the MAGICC ensemble panels"""

ENSEMBLE_SEGMENT_SHARES = (2.0 / 3.0, 1.0 / 3.0)
"""Share of the MAGICC ensemble panels' width given to each piece"""

SEGMENT_SHARES = (0.6, 0.4)
"""Share of a broken panel's width given to each piece"""

HARMONISATION_SEGMENTS = ((2010.0, 2035.0),)
"""Time axis of the panels which focus on the harmonisation (not broken)

Short, so the harmonisation itself can be seen,
and in its own panel, so it gets its own vertical axis.
Starts before 2015, so VL-CF's harmonisation (in 2015) can be seen too.
"""

HARMONISATION_TICKS = (2010.0, 2015.0, 2020.0, 2025.0, 2030.0, 2035.0)
"""Major ticks of [HARMONISATION_SEGMENTS][]"""

PROJECTION_SEGMENTS = ((2015.0, BREAK_YEAR), (BREAK_YEAR, END_YEAR))
"""Time axis pieces for the panels which show the harmonised projections in full

From 2015, so the start of VL-CF isn't cut off.
"""

PROJECTION_TICKS = (2015.0, 2060.0, 2100.0)
"""Major ticks of the first piece of [PROJECTION_SEGMENTS][]

Set explicitly, as the automatic choice can crowd the narrow piece.
"""

HARMONISATION_PANEL_WIDTHS = (0.88, 1.12)
"""Share of the row's width given to the harmonisation and projection panels"""

EXTENSION_TICKS = (2300.0, 2500.0)
"""Major ticks of the piece after the break"""

EXTENSION_MINOR_TICKS = (2200.0, 2400.0)
"""Minor ticks of the piece after the break"""

RANGES_YEARS = (2050, 2100, 2300, 2500)
"""Years shown in the panels of MAGICC's medians and ranges"""

RANGES_GROUP_WIDTH = 0.94
"""Share of the space between each year's ticks which that year's boxes fill"""

RANGES_LEGEND_HEADROOM = 0.2
"""Space above the boxes, as a share of their span, which leaves room for the legend"""

AR6_VERSION_ALPHA = 0.45
"""Opacity of the boxes of the MAGICC version used in AR6"""

ENSEMBLE_RANGE = (0.05, 0.95)
"""Quantiles which give the range of MAGICC's ensemble we show"""

HISTORY_LABEL = "Historical dataset"
"""Legend label of the historical dataset"""

RAW_LINESTYLE = "--"
"""Line style of raw (unharmonised) MAGICC output"""

LINE_WIDTH = 0.8
"""Width of most lines"""

HISTORY_LINE_WIDTH = 1.0
"""Width of the history's line"""

HISTORY_ZORDER = 3.0
"""Z-order of the history's line (over the scenarios which start where it ends)"""

LEAVES_HISTORY_EARLY_ZORDER = 3.5
"""Z-order of scenarios which leave history before it ends (over history)"""

PC_LINESTYLES = ("-", ":")
"""Line style of each PC (first, second)"""

PC_NAMES = ("PC0", "PC1")
"""Name of each PC in legends"""

LEGEND_PANEL = "legend"
"""Name of every figure's legend panel"""

ROW_HEIGHT = 0.95
"""Height of a row of time series panels, in inches"""

LEGEND_ROW_HEIGHT = 0.25
"""Height of the legend row, in inches, if the legend has one line"""

LEGEND_LINE_HEIGHT = 0.15
"""Height each extra line adds to the legend row, in inches"""

ALL_PANELS = "all panels"
"""How a legend group's title says it applies to every panel"""

LAYOUT_SETTINGS = LayoutSettings(align_right=True)
"""Layout settings (panels in a column share time axes, so line up both ends)"""

EXTENSION_SHADE_COLOUR = "0.92"
"""Colour of the shading which marks our extension of Western et al. (2024)"""

WESTERN_ET_AL_2024_GASES = ("hcfc22", "hcfc141b", "hcfc142b")
"""Gases whose projections come from Western et al. (2024), which we extend"""


def broken(
    name: str,
    width: float = 1.0,
    shares: tuple[float, float] = SEGMENT_SHARES,
) -> Panel:
    """Get a panel which is a time axis broken at [BREAK_YEAR][]"""
    return Panel(name, width=width, broken=True, broken_split=shares)


def legend_row(n_lines: int = 1) -> Row:
    """Get the row which holds a figure's legend, which has `n_lines` lines"""
    return Row(
        panels=(Panel(LEGEND_PANEL),),
        height=LEGEND_ROW_HEIGHT + (n_lines - 1) * LEGEND_LINE_HEIGHT,
    )


def to_frame(values: pd.Series) -> pd.DataFrame:
    """Get a series indexed by decimal year as the frame the plotting helpers take"""
    return pd.DataFrame(
        {
            TIME_COLUMN: values.index.to_numpy(dtype=float),
            VALUE_COLUMN: values.to_numpy(dtype=float),
        }
    )


def plot_in_segments(
    values: pd.Series,
    axes: Sequence[matplotlib.axes.Axes],
    segments: Sequence[tuple[float, float]],
    **kwargs: object,
) -> matplotlib.lines.Line2D:
    """Plot a series on every piece of a (possibly broken) time axis"""
    kwargs.setdefault("linewidth", LINE_WIDTH)
    return plot_line_in_segments(to_frame(values), axes, segments, **kwargs)


def fill_in_segments(
    lower: pd.Series,
    upper: pd.Series,
    axes: Sequence[matplotlib.axes.Axes],
    segments: Sequence[tuple[float, float]],
    **kwargs: object,
) -> None:
    """Shade between two series on every piece of a (possibly broken) time axis"""
    for ax, segment in zip(axes, segments):
        mask = in_segment(lower.index.to_numpy(dtype=float), segment)
        ax.fill_between(
            lower.index[mask], lower[mask], upper[mask], linewidth=0.0, **kwargs
        )


def finish_time_axes(
    axes: Sequence[matplotlib.axes.Axes],
    segments: Sequence[tuple[float, float]],
    ylabel: str,
    foreground: Sequence[pd.Series],
) -> None:
    """
    Set the scale, ticks and labels of a (possibly broken) time axis

    Parameters
    ----------
    axes
        Axes of each piece of the time axis

    segments
        Each piece of the time axis, as (start, end)

    ylabel
        Label of the vertical axis

    foreground
        Everything drawn which must be on scale
    """
    frames = [to_frame(values) for values in foreground]
    in_view = [
        df[(df[TIME_COLUMN] >= segments[0][0]) & (df[TIME_COLUMN] <= segments[-1][1])]
        for df in frames
    ]
    # One scale for all the pieces, so lines can be followed across the break
    set_segment_y_limits(axes[:1], [(-np.inf, np.inf)], in_view)
    set_up_segments(axes, segments, units="", xlabel="", share_y=True)
    for i, (ax, segment) in enumerate(zip(axes, segments)):
        if 0 < i and segment[0] < BREAK_YEAR:
            # The previous piece already labels where this one starts
            ax.set_xticks([t for t in ax.get_xticks() if t > segment[0]])
            ax.set_xlim(segment)

        if segment[1] == BREAK_YEAR:
            # Label where the extensions start, rather than the break
            ax.set_xticks(get_ticks_with_ends(segment[0], EXTENSIONS_START_YEAR))
            ax.set_xlim(segment)

        if segment[0] >= BREAK_YEAR:
            ax.set_xticks(EXTENSION_TICKS)
            ax.set_xticks(EXTENSION_MINOR_TICKS, minor=True)
            # Setting the ticks can move the limits
            ax.set_xlim(segment)

        ax.xaxis.set_major_formatter("{x:.0f}")

    axes[0].set_ylabel(ylabel, fontsize="small")


def join_to_history(values: pd.Series, history: pd.Series) -> pd.Series:
    """
    Join a projection which leaves history before history ends onto history

    The published VL-CF data starts the year after it leaves history,
    so without this there would be a gap between history and VL-CF.
    Projections which start at the end of history already overlap it,
    so are returned unchanged.
    """
    previous_year = values.index[0] - 1.0
    if values.index[0] > LAST_HISTORICAL_YEAR or previous_year not in history.index:
        return values

    return pd.concat([history.loc[[previous_year]], values])


def get_zorder(values: pd.Series) -> float:
    """
    Get the z-order of a scenario's line

    Projections which leave history before history ends (VL-CF)
    are drawn over the historical dataset's line,
    so it can be seen where they leave it.
    """
    if values.index[0] < LAST_HISTORICAL_YEAR:
        return LEAVES_HISTORY_EARLY_ZORDER

    return 2.0


def units_label(units: str) -> str:
    """Get the label of a vertical axis which shows values in some units"""
    return label_name(f"[{units}]").replace("dimensionless", "-")


LegendEntry = tuple[str, matplotlib.artist.Artist]
"""Label and handle of a legend entry"""

LegendGroup = tuple[str, Sequence[LegendEntry]]
"""Title and entries of a group of entries in a figure's legend panel"""


def add_legend(
    ax: matplotlib.axes.Axes, lines: Sequence[Sequence[LegendGroup]]
) -> None:
    """
    Put a figure's legend in its legend panel

    The legend panel explains what applies to more than one panel
    (what applies to only one panel is in that panel's legend,
    see [add_panel_legend][]).
    Its entries are in groups, each of which starts with a (bold) title
    that says what the group explains and which panels it applies to.

    Parameters
    ----------
    ax
        The legend panel's axes

    lines
        Groups on each line of the legend, top to bottom
        (the legend panel must be made with `legend_row(len(lines))`)
    """
    ax.axis("off")
    for i, line in enumerate(lines):
        handles: list[matplotlib.artist.Artist] = []
        labels: list[str] = []
        title_indexes = []
        for title, entries in line:
            title_indexes.append(len(labels))
            # A handle which draws nothing, so the title reads as a heading
            handles.append(matplotlib.lines.Line2D([], [], linestyle="none"))
            labels.append(f"{title}:")
            handles.extend(handle for _, handle in entries)
            labels.extend(label for label, _ in entries)

        legend = ax.legend(
            handles=handles,
            labels=labels,
            loc="center",
            # Lines spaced evenly from top to bottom
            bbox_to_anchor=(0.5, 1.0 - (i + 0.5) / len(lines)),
            ncols=len(labels),
            frameon=False,
            fontsize=LEGEND_FONT_SIZE,
            # Tight, so all the datasets fit on one line
            handlelength=1.4,
            handletextpad=0.4,
            columnspacing=0.6,
            borderaxespad=0.0,
        )
        texts = legend.get_texts()
        for j in title_indexes:
            texts[j].set_fontweight("bold")

        # See `add_scenario_legend` in the results figure for why
        legend.set_in_layout(False)
        ax.add_artist(legend)


def add_panel_legend(
    ax: matplotlib.axes.Axes,
    entries: Sequence[LegendEntry],
    loc: str = "best",
    ncols: int = 1,
) -> None:
    """
    Put a legend of what applies to only one panel in that panel

    Parameters
    ----------
    ax
        Axes to put the legend in (for a broken panel, one of its pieces)

    entries
        Label and handle of each entry, in order

    loc
        Where to put the legend
        (`"best"` for panels which are drawn for many gases,
        so the space which is free changes)

    ncols
        Number of columns
    """
    legend = ax.legend(
        handles=[handle for _, handle in entries],
        labels=[label for label, _ in entries],
        loc=loc,
        ncols=ncols,
        frameon=False,
        fontsize=LEGEND_FONT_SIZE,
        handlelength=1.6,
        handletextpad=0.5,
        borderaxespad=0.3,
        labelspacing=0.3,
    )
    legend.set_in_layout(False)


def datasets_group(
    panels: str = ALL_PANELS, scenarios: Sequence = ALL_SCENARIOS
) -> LegendGroup:
    """
    Get the legend group which says which colour each dataset is

    The historical dataset is included,
    so it is clear that the colours distinguish datasets.

    Parameters
    ----------
    panels
        Panels the group applies to, as written in its title

    scenarios
        Scenarios whose colours to explain, in order
    """
    return (
        f"Datasets ({panels})",
        [(HISTORY_LABEL, history_handle()), *scenario_entries(scenarios)],
    )


def line_handle(
    colour: str, linestyle: str = "-", linewidth: float = LINE_WIDTH
) -> matplotlib.lines.Line2D:
    """Get a legend handle for a line"""
    return matplotlib.lines.Line2D(
        [], [], color=colour, linestyle=linestyle, linewidth=linewidth
    )


def history_handle() -> matplotlib.lines.Line2D:
    """Get the legend handle for the historical dataset"""
    return line_handle(HISTORY_COLOUR, linewidth=HISTORY_LINE_WIDTH)


def scenario_entries(
    scenarios: Sequence = ALL_SCENARIOS,
) -> list[tuple[str, matplotlib.artist.Artist]]:
    """Get the legend entries which say which colour each scenario is"""
    return [(s.label, line_handle(s.colour)) for s in scenarios]


def get_panel_axes(
    rows: Sequence[Row], axes: dict[str, matplotlib.axes.Axes], name: str
) -> list[matplotlib.axes.Axes]:
    """Get the axes which make up a panel"""
    (panel,) = (p for row in rows for p in row.panels if p.name == name)
    return [axes[n] for n in get_panel_axes_names(panel)]


def save(
    fig: plt.Figure,
    axes: dict[str, matplotlib.axes.Axes],
    rows: Sequence[Row],
    titles: dict[str, str],
    outfile: Path,
) -> Path:
    """Label, lay out and save a figure, then close it"""
    label_panels(rows, axes, titles, unlabelled=(LEGEND_PANEL,))
    # Last, because it needs to know how much room everything takes up
    lay_out_figure(fig, axes, rows, settings=LAYOUT_SETTINGS)

    outfile.parent.mkdir(exist_ok=True, parents=True)
    fig.savefig(outfile)
    plt.close(fig)

    return outfile


def plot_magicc_ensemble(
    axes: Sequence[matplotlib.axes.Axes], gas: str, summary_file: Path
) -> None:
    """
    Plot MAGICC's ensemble (median and range) for L, with the historical dataset

    The MAGICC output is raw, i.e. before harmonisation.
    """
    units = get_units(gas)
    scenario = SCENARIOS[SHOWN_SCENARIO]
    median = get_magicc_concentrations(summary_file, gas, SHOWN_SCENARIO, 0.5, units)
    lower, upper = (
        get_magicc_concentrations(summary_file, gas, SHOWN_SCENARIO, q, units)
        for q in ENSEMBLE_RANGE
    )
    history = get_history(gas, units)

    fill_in_segments(
        lower, upper, axes, ENSEMBLE_SEGMENTS, color=scenario.colour, alpha=0.25
    )
    plot_in_segments(median, axes, ENSEMBLE_SEGMENTS, color=scenario.colour)
    plot_in_segments(
        history,
        axes,
        ENSEMBLE_SEGMENTS,
        color=HISTORY_COLOUR,
        linewidth=HISTORY_LINE_WIDTH,
    )
    finish_time_axes(
        axes, ENSEMBLE_SEGMENTS, units_label(units), [lower, upper, history]
    )


ENSEMBLE_TITLE = "MAGICC ensemble (raw), L"
"""Title of the MAGICC ensemble panels"""


def ensemble_entries() -> list[LegendEntry]:
    """Get the legend entries for the MAGICC ensemble panel"""
    colour = SCENARIOS[SHOWN_SCENARIO].colour
    return [
        (HISTORY_LABEL, history_handle()),
        ("MAGICC median", line_handle(colour)),
        (
            "MAGICC 5\u201395% range",
            matplotlib.patches.Patch(facecolor=colour, alpha=0.25, linewidth=0.0),
        ),
    ]


def get_ranges_bottom(gas: str, units: str, lowest: float) -> float:
    """
    Get the bottom of the vertical axis of a panel of MAGICC's medians and ranges

    The pre-industrial concentration, rounded down to two significant figures,
    rather than zero, so the boxes get the panel's height,
    unless the boxes go lower than that,
    in which case the bottom of the lowest box (`lowest`).
    """
    pre_industrial = get_history(gas, units).iloc[0]
    if pre_industrial <= 0.0:
        return 0.0

    rounding = 10.0 ** (np.floor(np.log10(pre_industrial)) - 1.0)

    return float(min(np.floor(pre_industrial / rounding) * rounding, lowest))


def plot_magicc_ranges(
    ax: matplotlib.axes.Axes,
    gas: str,
    summary_file: Path,
    climate_models: Sequence[tuple[str, float]],
) -> None:
    """
    Plot MAGICC's raw median and 5-95% range, all scenarios, in [RANGES_YEARS][]

    Each year's boxes fill most of the space between its tick and the next,
    with one box per scenario and MAGICC version.

    Parameters
    ----------
    ax
        Axes to plot on

    gas
        Gas of interest

    summary_file
        File which holds the MAGICC summary

    climate_models
        MAGICC versions to show, with the opacity of their boxes,
        in the order they are shown within each scenario's slot
    """
    units = get_units(gas)
    n_boxes = len(ALL_SCENARIOS) * len(climate_models)
    box_width = RANGES_GROUP_WIDTH / n_boxes
    lowest, highest = np.inf, -np.inf
    for i, year in enumerate(RANGES_YEARS):
        for j, scenario in enumerate(ALL_SCENARIOS):
            for k, (climate_model, alpha) in enumerate(climate_models):
                median, lower, upper = (
                    get_magicc_concentrations(
                        summary_file, gas, scenario.name, q, units, climate_model
                    ).loc[year + 0.5]
                    for q in (0.5, *ENSEMBLE_RANGE)
                )
                x = (
                    i
                    - RANGES_GROUP_WIDTH / 2
                    + (j * len(climate_models) + k + 0.5) * box_width
                )
                # A box from the 5th to the 95th percentile, with the median across it
                ax.bar(
                    x,
                    upper - lower,
                    bottom=lower,
                    width=box_width,
                    color=scenario.colour,
                    alpha=alpha,
                    linewidth=0.0,
                )
                ax.hlines(
                    median,
                    x - box_width / 2,
                    x + box_width / 2,
                    color="0.1",
                    linewidth=0.6,
                )
                lowest = min(lowest, lower)
                highest = max(highest, upper)

    bottom = get_ranges_bottom(gas, units, lowest - 0.03 * (highest - lowest))
    ax.set_ylim(bottom, highest + RANGES_LEGEND_HEADROOM * (highest - bottom))
    ax.set_xticks(range(len(RANGES_YEARS)))
    ax.set_xticklabels([str(y) for y in RANGES_YEARS])
    ax.set_xlim(-0.5, len(RANGES_YEARS) - 0.5)
    ax.tick_params(labelsize="small")
    ax.set_ylabel(units_label(units), fontsize="small")


def ranges_entries(climate_models: Sequence[tuple[str, float]]) -> list[LegendEntry]:
    """Get the legend entries for the panel made by [plot_magicc_ranges][]"""
    boxes = [
        (
            f"{climate_model} 5\u201395% range",
            matplotlib.patches.Patch(facecolor="0.4", alpha=alpha, linewidth=0.0),
        )
        for climate_model, alpha in climate_models
    ]

    return [*boxes, ("MAGICC median", line_handle("0.1", linewidth=0.6))]


def harmonisation_panels(
    harmonisation: str = "harmonisation", projections: str = "projections"
) -> tuple[Panel, Panel]:
    """
    Get the panels which show the harmonisation, for one row

    The first focuses on the harmonisation ([HARMONISATION_SEGMENTS][]),
    the second shows the projections in full ([PROJECTION_SEGMENTS][]).
    """
    harmonisation_width, projections_width = HARMONISATION_PANEL_WIDTHS
    return (
        Panel(harmonisation, width=harmonisation_width),
        broken(projections, width=projections_width),
    )


HARMONISATION_TITLES = {
    "harmonisation": "Global-mean and harmonisation, 2010\u20132035",
    "projections": "Global-mean and harmonisation, 2015\u20132500",
}
"""Titles of the panels made by [harmonisation_panels][]"""


def plot_harmonisation(
    harmonisation_axes: Sequence[matplotlib.axes.Axes],
    projections_axes: Sequence[matplotlib.axes.Axes],
    gas: str,
    summary_file: Path,
) -> None:
    """
    Plot the raw MAGICC median, the harmonised projection and history, all scenarios

    The harmonised projection is the published data.
    It is plotted twice: once focussed on the harmonisation ([HARMONISATION_SEGMENTS][])
    and once in full ([PROJECTION_SEGMENTS][]).
    The panels must be made with [harmonisation_panels][].
    """
    units = get_units(gas)
    history = get_history(gas, units)
    foreground = [history]
    for scenario in ALL_SCENARIOS:
        raw = get_magicc_concentrations(summary_file, gas, scenario.name, 0.5, units)
        harmonised = join_to_history(get_scenario(gas, scenario.name, units), history)
        for axes, segments in (
            (harmonisation_axes, HARMONISATION_SEGMENTS),
            (projections_axes, PROJECTION_SEGMENTS),
        ):
            plot_in_segments(
                raw,
                axes,
                segments,
                color=scenario.colour,
                linestyle=RAW_LINESTYLE,
                linewidth=0.6,
            )
            plot_in_segments(
                harmonised,
                axes,
                segments,
                color=scenario.colour,
                zorder=get_zorder(harmonised),
            )

        foreground.extend([raw, harmonised])

    for axes, segments in (
        (harmonisation_axes, HARMONISATION_SEGMENTS),
        (projections_axes, PROJECTION_SEGMENTS),
    ):
        plot_in_segments(
            history,
            axes,
            segments,
            color=HISTORY_COLOUR,
            linewidth=HISTORY_LINE_WIDTH,
            zorder=HISTORY_ZORDER,
        )
        finish_time_axes(axes, segments, units_label(units), foreground)

    (ax,) = harmonisation_axes
    ax.set_xticks(HARMONISATION_TICKS)
    # Setting the ticks can move the limits
    ax.set_xlim(HARMONISATION_SEGMENTS[0])
    projections_axes[0].set_xticks(PROJECTION_TICKS)
    projections_axes[0].set_xlim(PROJECTION_SEGMENTS[0])


def harmonisation_entries() -> list[LegendEntry]:
    """Get the legend entries for the line styles of the harmonisation panel"""
    return [
        ("CMIP7 forcings (harmonised)", line_handle("0.4")),
        ("Raw MAGICC median", line_handle("0.4", linestyle=RAW_LINESTYLE)),
    ]


def plot_pcs(  # noqa: PLR0913
    axes: Sequence[matplotlib.axes.Axes],
    historical: pd.DataFrame | None,
    projected: pd.DataFrame,
    scenarios: Sequence[str] = (SHOWN_SCENARIO,),
    segments: Sequence[tuple[float, float]] = LONG_SEGMENTS,
    units: str = "-",
) -> None:
    """
    Plot PCs over time: the historical dataset's and the scenarios' projections

    Parameters
    ----------
    axes
        Axes of each piece of the time axis

    historical
        Historical PCs, one column per EOF (`None` if there are none)

    projected
        Projected PCs, with columns `(eof, scenario)` or just `scenario`

    scenarios
        Scenarios to show, each in its colour
        (if the projection is the same in every scenario,
        it has only one scenario, which is shown in L's colour whatever this is)

    segments
        Each piece of the time axis

    units
        Units of the PCs (`"-"` if dimensionless)
    """
    if not isinstance(projected.columns, pd.MultiIndex):
        projected = pd.concat({0: projected}, axis="columns")

    foreground = []
    for eof in projected.columns.get_level_values(0).unique():
        eof_projected = projected[eof]
        if eof_projected.shape[1] == 1:
            to_plot = {SHOWN_SCENARIO: eof_projected.iloc[:, 0]}
        else:
            to_plot = {scenario: eof_projected[scenario] for scenario in scenarios}

        linestyle = PC_LINESTYLES[eof]
        for scenario, values in to_plot.items():
            plot_in_segments(
                values,
                axes,
                segments,
                color=SCENARIOS[scenario].colour,
                linestyle=linestyle,
                zorder=get_zorder(values.dropna()),
            )
            foreground.append(values)

        if historical is not None:
            plot_in_segments(
                historical[eof],
                axes,
                segments,
                color=HISTORY_COLOUR,
                linestyle=linestyle,
                linewidth=HISTORY_LINE_WIDTH,
            )
            foreground.append(historical[eof])

    finish_time_axes(axes, segments, f"PC value [{units}]", foreground)


def add_pc_legend(
    axes: Sequence[matplotlib.axes.Axes], n_eofs: int, loc: str = "best"
) -> None:
    """Put the legend of the PCs' line styles in a PC panel, if it has two PCs"""
    if n_eofs == 1:
        return

    add_panel_legend(
        axes[0],
        [
            (name, line_handle("0.4", linestyle=linestyle))
            for name, linestyle in zip(PC_NAMES[:n_eofs], PC_LINESTYLES)
        ],
        loc=loc,
    )


@manuscript_style
def generate_co2_figure(
    outfile: Path, summary_file: Path, bundle_dir: Path = BUNDLE_DIR
) -> Path:
    """
    Generate the CO2 methods figure

    Parameters
    ----------
    outfile
        File in which to write the figure

    summary_file
        File which holds the MAGICC summary

    bundle_dir
        Directory which holds the bundle (including our re-runs)

    Returns
    -------
    :
        `outfile`
    """
    rows = (
        Row(
            panels=(broken("ensemble", shares=ENSEMBLE_SEGMENT_SHARES),),
            height=ROW_HEIGHT,
        ),
        Row(panels=harmonisation_panels(), height=ROW_HEIGHT),
        Row(panels=(broken("seasonality-pc"),), height=ROW_HEIGHT),
        Row(panels=(broken("lat-gradient-pcs"),), height=ROW_HEIGHT),
        legend_row(2),
    )
    titles = {
        "ensemble": ENSEMBLE_TITLE,
        **HARMONISATION_TITLES,
        "seasonality-pc": "Seasonality change PC",
        "lat-gradient-pcs": "Latitudinal gradient PCs",
    }
    fig, axes = create_figure(rows, settings=LAYOUT_SETTINGS)

    def panel(name: str) -> list[matplotlib.axes.Axes]:
        return get_panel_axes(rows, axes, name)

    plot_magicc_ensemble(panel("ensemble"), "co2", summary_file)
    plot_harmonisation(
        panel("harmonisation"), panel("projections"), "co2", summary_file
    )

    scenario_names = [s.name for s in ALL_SCENARIOS]
    plot_pcs(
        panel("seasonality-pc"),
        get_historical_pcs("co2", kind="seasonality-change"),
        get_projected_co2_seasonality_change_pc(bundle_dir),
        scenario_names,
    )
    plot_pcs(
        panel("lat-gradient-pcs"),
        get_historical_pcs("co2"),
        get_projected_pcs("co2", bundle_dir, include_vl_cf=True),
        scenario_names,
    )

    # After the break, where the plume is narrowest
    add_panel_legend(panel("ensemble")[1], ensemble_entries(), loc="upper right")
    add_pc_legend(panel("lat-gradient-pcs"), 2, loc="lower left")
    add_legend(
        axes[LEGEND_PANEL],
        [[("Panels b, c", harmonisation_entries())], [datasets_group()]],
    )

    return save(fig, axes, rows, titles, outfile)


@manuscript_style
def generate_ch4_figure(
    outfile: Path, summary_file: Path, bundle_dir: Path = BUNDLE_DIR
) -> Path:
    """
    Generate the CH4 methods figure

    Parameters
    ----------
    outfile
        File in which to write the figure

    summary_file
        File which holds the MAGICC summary

    bundle_dir
        Directory which holds the bundle (including our re-runs)

    Returns
    -------
    :
        `outfile`
    """
    rows = (
        Row(panels=(broken("versions"),), height=ROW_HEIGHT),
        Row(panels=(Panel("ranges"),), height=ROW_HEIGHT),
        Row(panels=harmonisation_panels(), height=ROW_HEIGHT),
        Row(panels=(broken("lat-gradient-pcs"),), height=ROW_HEIGHT),
        legend_row(2),
    )
    titles = {
        "versions": (
            f"Projections from {MAGICC_VERSION} compared to {MAGICC_AR6_VERSION}"
        ),
        "ranges": "Median and ranges",
        **HARMONISATION_TITLES,
        "lat-gradient-pcs": "Latitudinal gradient PCs",
    }
    fig, axes = create_figure(rows, settings=LAYOUT_SETTINGS)

    def panel(name: str) -> list[matplotlib.axes.Axes]:
        return get_panel_axes(rows, axes, name)

    units = get_units("ch4")

    # MAGICC versions over time
    foreground = []
    for scenario in ALL_SCENARIOS:
        for climate_model, linestyle in (
            (MAGICC_VERSION, "-"),
            (MAGICC_AR6_VERSION, RAW_LINESTYLE),
        ):
            median = get_magicc_concentrations(
                summary_file, "ch4", scenario.name, 0.5, units, climate_model
            )
            plot_in_segments(
                median,
                panel("versions"),
                CH4_VERSIONS_SEGMENTS,
                color=scenario.colour,
                linestyle=linestyle,
            )
            foreground.append(median)

    finish_time_axes(
        panel("versions"), CH4_VERSIONS_SEGMENTS, units_label(units), foreground
    )

    # MAGICC versions in selected years
    climate_models = ((MAGICC_AR6_VERSION, AR6_VERSION_ALPHA), (MAGICC_VERSION, 1.0))
    (ranges_ax,) = panel("ranges")
    plot_magicc_ranges(ranges_ax, "ch4", summary_file, climate_models)

    plot_harmonisation(
        panel("harmonisation"), panel("projections"), "ch4", summary_file
    )
    plot_pcs(
        panel("lat-gradient-pcs"),
        get_historical_pcs("ch4"),
        get_projected_pcs("ch4", bundle_dir, include_vl_cf=True),
        [s.name for s in ALL_SCENARIOS],
    )

    add_panel_legend(
        panel("versions")[0],
        [
            (MAGICC_VERSION, line_handle("0.4")),
            (MAGICC_AR6_VERSION, line_handle("0.4", linestyle=RAW_LINESTYLE)),
        ],
        loc="upper left",
    )
    # On one line, above the boxes
    add_panel_legend(
        ranges_ax, ranges_entries(climate_models), loc="upper left", ncols=3
    )
    add_pc_legend(panel("lat-gradient-pcs"), 2)
    add_legend(
        axes[LEGEND_PANEL],
        [[("Panels c, d", harmonisation_entries())], [datasets_group()]],
    )

    return save(fig, axes, rows, titles, outfile)


@manuscript_style
def generate_n2o_like_figure(
    gas: str, outfile: Path, summary_file: Path, bundle_dir: Path = BUNDLE_DIR
) -> Path:
    """
    Generate the methods figure of a gas processed like N2O

    Laid out like the CH4 figure ([generate_ch4_figure][]),
    except that the comparison of MAGICC versions is replaced
    by the MAGICC ensemble panel of the CO2 figure ([generate_co2_figure][]).

    Parameters
    ----------
    gas
        Gas of interest

    outfile
        File in which to write the figure

    summary_file
        File which holds the MAGICC summary

    bundle_dir
        Directory which holds the bundle (including our re-runs)

    Returns
    -------
    :
        `outfile`
    """
    rows = (
        Row(
            panels=(broken("ensemble", shares=ENSEMBLE_SEGMENT_SHARES),),
            height=ROW_HEIGHT,
        ),
        Row(panels=harmonisation_panels(), height=ROW_HEIGHT),
        Row(panels=(broken("lat-gradient-pcs"),), height=ROW_HEIGHT),
        legend_row(2),
    )
    projected_pcs = get_projected_pcs(gas, bundle_dir, include_vl_cf=True)
    n_eofs = projected_pcs.columns.get_level_values(0).nunique()
    titles = {
        "ensemble": ENSEMBLE_TITLE,
        **HARMONISATION_TITLES,
        "lat-gradient-pcs": (
            "Latitudinal gradient PCs" if n_eofs > 1 else "Latitudinal gradient PC"
        ),
    }
    fig, axes = create_figure(rows, settings=LAYOUT_SETTINGS)

    def panel(name: str) -> list[matplotlib.axes.Axes]:
        return get_panel_axes(rows, axes, name)

    plot_magicc_ensemble(panel("ensemble"), gas, summary_file)
    plot_harmonisation(panel("harmonisation"), panel("projections"), gas, summary_file)
    plot_pcs(
        panel("lat-gradient-pcs"),
        get_historical_pcs(gas),
        projected_pcs,
        [s.name for s in ALL_SCENARIOS],
    )

    add_panel_legend(panel("ensemble")[0], ensemble_entries())
    add_pc_legend(panel("lat-gradient-pcs"), n_eofs)
    add_legend(
        axes[LEGEND_PANEL],
        [[("Panels b, c", harmonisation_entries())], [datasets_group()]],
    )

    return save(fig, axes, rows, titles, outfile)


def shade_western_et_al_2024_extension(
    gas: str,
    axes: Sequence[matplotlib.axes.Axes],
    segments: Sequence[tuple[float, float]],
) -> None:
    """
    Mark where we extended the Western et al. (2024) projection, if the gas has one
    """
    if gas not in WESTERN_ET_AL_2024_GASES:
        return

    for ax, segment in zip(axes, segments):
        if segment[1] <= EXTENSIONS_START_YEAR:
            continue

        ax.axvspan(
            max(segment[0], EXTENSIONS_START_YEAR),
            segment[1],
            color=EXTENSION_SHADE_COLOUR,
            zorder=0.0,
        )
        if segment[0] < BREAK_YEAR:
            # Only label the shading after the break, where there is room
            continue

        ax.text(
            0.5,
            0.95,
            "Our extension of\nWestern et al. (2024)",
            transform=ax.transAxes,
            ha="center",
            va="top",
            fontsize=LEGEND_FONT_SIZE,
        )


def get_emissions_units(gas: str) -> str:
    """Get the units in which we show a gas' emissions"""
    return f"kt {MAGICC_GAS_NAMES[gas]} / yr"


@manuscript_style
def generate_cfc12_like_figure(
    gas: str, outfile: Path, bundle_dir: Path = BUNDLE_DIR
) -> Path:
    """
    Generate the methods figure of a gas processed like CFC-12

    Parameters
    ----------
    gas
        Gas of interest

    outfile
        File in which to write the figure

    bundle_dir
        Directory which holds the bundle

    Returns
    -------
    :
        `outfile`
    """
    rows = (
        Row(panels=(broken("projection"),), height=ROW_HEIGHT),
        Row(panels=(broken("lat-gradient-pcs"),), height=ROW_HEIGHT),
        legend_row(),
    )
    titles = {
        "projection": "Global-, annual-mean",
        "lat-gradient-pcs": "Latitudinal gradient PC",
    }
    fig, axes = create_figure(rows, settings=LAYOUT_SETTINGS)

    def panel(name: str) -> list[matplotlib.axes.Axes]:
        return get_panel_axes(rows, axes, name)

    units = get_units(gas)
    colour = SCENARIOS[SHOWN_SCENARIO].colour

    projection = get_scenario(gas, SHOWN_SCENARIO, units)
    history = get_history(gas, units)
    shade_western_et_al_2024_extension(gas, panel("projection"), SCENARIO_SEGMENTS)
    plot_in_segments(projection, panel("projection"), SCENARIO_SEGMENTS, color=colour)
    plot_in_segments(
        history,
        panel("projection"),
        SCENARIO_SEGMENTS,
        color=HISTORY_COLOUR,
        linewidth=HISTORY_LINE_WIDTH,
    )
    finish_time_axes(
        panel("projection"),
        SCENARIO_SEGMENTS,
        units_label(units),
        [projection, history],
    )

    # The projection doesn't depend on the scenario,
    # so this is every scenario's PC
    plot_pcs(
        panel("lat-gradient-pcs"),
        get_historical_pcs(gas),
        get_projected_pcs(gas, bundle_dir),
    )

    add_legend(
        axes[LEGEND_PANEL],
        [
            [
                (
                    f"Datasets ({ALL_PANELS})",
                    [
                        (HISTORY_LABEL, history_handle()),
                        (
                            "CMIP7 forcings (the same in all datasets)",
                            line_handle(colour),
                        ),
                    ],
                )
            ]
        ],
    )

    return save(fig, axes, rows, titles, outfile)


@manuscript_style
def generate_hfc125_like_figure(
    gas: str, outfile: Path, summary_file: Path, bundle_dir: Path = BUNDLE_DIR
) -> Path:
    """
    Generate the methods figure of a gas processed like HFC-125

    Parameters
    ----------
    gas
        Gas of interest

    outfile
        File in which to write the figure

    summary_file
        File which holds the MAGICC summary

    bundle_dir
        Directory which holds the bundle

    Returns
    -------
    :
        `outfile`
    """
    rows = (
        Row(panels=(broken("one-box"),), height=ROW_HEIGHT),
        Row(panels=(broken("lat-gradient-pcs"),), height=ROW_HEIGHT),
        legend_row(),
    )
    titles = {
        "one-box": (
            "Global-, annual-mean projections and comparison with MAGICC output"
        ),
        "lat-gradient-pcs": "Latitudinal gradient PC",
    }
    fig, axes = create_figure(rows, settings=LAYOUT_SETTINGS)

    def panel(name: str) -> list[matplotlib.axes.Axes]:
        return get_panel_axes(rows, axes, name)

    units = get_units(gas)
    history = get_history(gas, units)
    foreground = [history]
    for scenario in ALL_SCENARIOS:
        one_box = join_to_history(get_scenario(gas, scenario.name, units), history)
        magicc = get_magicc_concentrations(summary_file, gas, scenario.name, 0.5, units)
        plot_in_segments(
            magicc,
            panel("one-box"),
            ONE_BOX_SEGMENTS,
            color=scenario.colour,
            linestyle=RAW_LINESTYLE,
            linewidth=0.6,
        )
        plot_in_segments(
            one_box,
            panel("one-box"),
            ONE_BOX_SEGMENTS,
            color=scenario.colour,
            zorder=get_zorder(one_box),
        )
        foreground.extend([one_box, magicc])

    plot_in_segments(
        history,
        panel("one-box"),
        ONE_BOX_SEGMENTS,
        color=HISTORY_COLOUR,
        linewidth=HISTORY_LINE_WIDTH,
        zorder=HISTORY_ZORDER,
    )
    finish_time_axes(panel("one-box"), ONE_BOX_SEGMENTS, units_label(units), foreground)

    plot_pcs(
        panel("lat-gradient-pcs"),
        get_historical_pcs(gas),
        get_projected_pcs(gas, bundle_dir, include_vl_cf=True),
        [s.name for s in ALL_SCENARIOS],
    )

    add_panel_legend(
        panel("one-box")[0],
        [
            ("CMIP7 forcings (one-box model)", line_handle("0.4")),
            ("MAGICC median", line_handle("0.4", linestyle=RAW_LINESTYLE)),
        ],
    )
    add_legend(axes[LEGEND_PANEL], [[datasets_group()]])

    return save(fig, axes, rows, titles, outfile)


REGRESSION_PANELS: dict[str, dict[str, tuple[str, str, str]]] = {
    "co2": {
        "seasonality-change": (
            CO2_SEASONALITY_REGRESSION,
            "Seasonality change PC vs. NPP",
            "MAGICC median NPP",
        ),
        "lat-gradient": (
            "co2",
            "Latitudinal gradient PC0 vs. emissions",
            "Fossil co2 emissions",
        ),
    },
    "ch4": {
        "lat-gradient": (
            "ch4",
            "Latitudinal gradient PC0 vs. emissions",
            "Fossil ch4 emissions",
        ),
    },
    "n2o": {
        "lat-gradient": (
            "n2o",
            "Latitudinal gradient PC0 vs. emissions",
            "Total n2o emissions",
        ),
    },
}
"""The panels of each gas' regressions figure

For each panel, the regression it shows (see [get_regression][]),
its title and the name of the quantity the PC is regressed against.
"""

REGRESSION_X_UNITS = {
    "co2": "GtC / yr",
    CO2_SEASONALITY_REGRESSION: "GtC / yr",
    "ch4": "MtCH4 / yr",
    "n2o": "ktN2O / yr",
}
"""Units in which to show the quantity each regression's PC is regressed against

The same as the historical dataset's figures, where they show the same thing.
"""

REGRESSION_ROW_HEIGHT = 1.4
"""Height of the row of a regressions figure, in inches"""


def get_regression_inputs(
    key: str, bundle_dir: Path = BUNDLE_DIR
) -> tuple[xr.Dataset, xr.Dataset, dict[str, tuple[float, str]]]:
    """
    Get a regression we re-derived in the form the historical paper's plotting takes

    See `plot_pc_timeseries_regression`.

    Parameters
    ----------
    key
        Regression of interest (see [get_regression][])

    bundle_dir
        Directory which holds the bundle (including our re-runs)

    Returns
    -------
    :
        PCs (the regression's `y`, as the first EOF's PC),
        the quantity they are regressed against (`x`)
        and the regression's result
    """
    data, info = get_regression(key, bundle_dir)
    years = data["year"].to_numpy(dtype=int)
    pcs = xr.Dataset(
        {
            "principal-components": xr.DataArray(
                data["y"].to_numpy()[np.newaxis, :],
                dims=("eof", "year"),
                coords={"eof": [0], "year": years},
                attrs={"units": info["y_unit"]},
            )
        }
    )
    timeseries = xr.Dataset(
        {
            "x": xr.DataArray(
                data["x"].to_numpy(),
                dims=("year",),
                coords={"year": years},
                attrs={"units": info["x_unit"]},
            )
        }
    )
    regression_info = {
        "m": (info["slope"], info["slope_unit"]),
        "c": (info["intercept"], info["y_unit"]),
    }

    return pcs, timeseries, regression_info


@manuscript_style
def generate_regressions_figure(
    gas: str, outfile: Path, bundle_dir: Path = BUNDLE_DIR
) -> Path:
    """
    Generate the figure of the regressions we re-derived for a gas

    The regressions are plotted exactly as in the historical dataset's figures.
    Figures with one panel are half the usual width.

    Parameters
    ----------
    gas
        Gas of interest (a key of [REGRESSION_PANELS][])

    outfile
        File in which to write the figure

    bundle_dir
        Directory which holds the bundle (including our re-runs)

    Returns
    -------
    :
        `outfile`
    """
    panels = REGRESSION_PANELS[gas]
    rows = (
        Row(
            panels=tuple(Panel(name) for name in panels),
            height=REGRESSION_ROW_HEIGHT,
        ),
    )
    titles = {name: label_name(title) for name, (_, title, _) in panels.items()}
    settings = dataclasses.replace(
        LAYOUT_SETTINGS, width=LAYOUT_SETTINGS.width * len(panels) / 2
    )
    fig, axes = create_figure(rows, settings=settings)

    for name, (key, _, x_name) in panels.items():
        pcs, timeseries, regression_info = get_regression_inputs(key, bundle_dir)
        plot_pc_timeseries_regression(
            pcs,
            timeseries,
            timeseries_name=label_name(x_name),
            regression_info=regression_info,
            ax=axes[name],
            x_unit=REGRESSION_X_UNITS[key],
        )

    if len(panels) == 1:
        # A figure with one panel doesn't need a panel label, just the title
        ((name, title),) = titles.items()
        axes[name].set_title(title, loc="left", fontsize=TITLE_FONT_SIZE)
        titles = {}
        unlabelled = (name,)
    else:
        unlabelled = ()

    return save_without_legend(
        fig, axes, rows, titles, outfile, settings=settings, unlabelled=unlabelled
    )


def save_without_legend(  # noqa: PLR0913
    fig: plt.Figure,
    axes: dict[str, matplotlib.axes.Axes],
    rows: Sequence[Row],
    titles: dict[str, str],
    outfile: Path,
    settings: LayoutSettings = LAYOUT_SETTINGS,
    unlabelled: Sequence[str] = (),
) -> Path:
    """Label, lay out and save a figure which has no legend panel, then close it"""
    label_panels(rows, axes, titles, unlabelled=unlabelled)
    lay_out_figure(fig, axes, rows, settings=settings)

    outfile.parent.mkdir(exist_ok=True, parents=True)
    fig.savefig(outfile)
    plt.close(fig)

    return outfile


@manuscript_style
def generate_western_et_al_2024_figure(
    outfile: Path, bundle_dir: Path = BUNDLE_DIR
) -> Path:
    """
    Generate the figure of our extensions of the Western et al. (2024) projections

    One row per gas: concentrations on the left, emissions on the right.
    The top row is annotated with the steps of the extension:
    invert the projection to get emissions,
    hold the emissions in 2100 constant,
    then run the one-box model forward with them.

    Parameters
    ----------
    outfile
        File in which to write the figure

    bundle_dir
        Directory which holds the bundle

    Returns
    -------
    :
        `outfile`
    """
    rows = (
        *(
            Row(
                panels=(broken(gas), broken(f"{gas}-emissions")),
                height=ROW_HEIGHT,
            )
            for gas in WESTERN_ET_AL_2024_GASES
        ),
        legend_row(2),
    )
    titles = {}
    for gas in WESTERN_ET_AL_2024_GASES:
        titles[gas] = f"{label_name(gas)} concentrations"
        titles[f"{gas}-emissions"] = f"{label_name(gas)} emissions"

    fig, axes = create_figure(rows, settings=LAYOUT_SETTINGS)

    wmo_colour = "0.55"
    western_colour = SCENARIOS[SHOWN_SCENARIO].colour
    for i, gas in enumerate(WESTERN_ET_AL_2024_GASES):
        gas_axes = get_panel_axes(rows, axes, gas)
        units = "ppt"
        raw, extended = get_western_et_al_2024(gas, bundle_dir)
        extension = extended[extended.index >= EXTENSIONS_START_YEAR + 0.5 - 1.0]
        wmo = get_wmo_2022(gas, bundle_dir)
        history = get_history(gas, units)

        plot_in_segments(wmo, gas_axes, SCENARIO_SEGMENTS, color=wmo_colour)
        plot_in_segments(raw, gas_axes, SCENARIO_SEGMENTS, color=western_colour)
        plot_in_segments(
            extension,
            gas_axes,
            SCENARIO_SEGMENTS,
            color=western_colour,
            linestyle=RAW_LINESTYLE,
        )
        plot_in_segments(
            history,
            gas_axes,
            SCENARIO_SEGMENTS,
            color=HISTORY_COLOUR,
            linewidth=HISTORY_LINE_WIDTH,
        )
        finish_time_axes(
            gas_axes,
            SCENARIO_SEGMENTS,
            units_label(units),
            [wmo, raw, extension, history],
        )

        emissions_axes = get_panel_axes(rows, axes, f"{gas}-emissions")
        emissions_units = get_emissions_units(gas)
        inverse, constant = get_western_et_al_2024_emissions(
            gas, emissions_units, bundle_dir
        )
        plot_in_segments(
            inverse, emissions_axes, SCENARIO_SEGMENTS, color=western_colour
        )
        plot_in_segments(
            constant,
            emissions_axes,
            SCENARIO_SEGMENTS,
            color=western_colour,
            linestyle=RAW_LINESTYLE,
        )
        finish_time_axes(
            emissions_axes,
            SCENARIO_SEGMENTS,
            units_label(emissions_units),
            [inverse, constant],
        )

        if i == 0:
            annotate_western_et_al_2024_steps(
                gas_axes, emissions_axes, inverse, constant, extension
            )

    add_legend(
        axes[LEGEND_PANEL],
        [
            [
                (
                    ALL_PANELS.capitalize(),
                    [
                        (
                            "Western et al. (2024) and its inverse emissions",
                            line_handle(western_colour),
                        ),
                        (
                            "Our extension",
                            line_handle(western_colour, linestyle=RAW_LINESTYLE),
                        ),
                    ],
                )
            ],
            [
                (
                    "Panels a, c, e",
                    [
                        (HISTORY_LABEL, history_handle()),
                        ("WMO (2022)", line_handle(wmo_colour)),
                    ],
                )
            ],
        ],
    )

    return save(fig, axes, rows, titles, outfile)


def annotate_step(
    ax: matplotlib.axes.Axes,
    text: str,
    xy: tuple[float, float],
    xytext: tuple[float, float],
) -> None:
    """
    Annotate a step of a process on a plot

    Parameters
    ----------
    ax
        Axes to annotate

    text
        Text of the annotation

    xy
        Point the annotation points to, in data coordinates

    xytext
        Where the text goes, in axes coordinates
    """
    ax.annotate(
        text,
        xy=xy,
        xytext=xytext,
        textcoords="axes fraction",
        fontsize=LEGEND_FONT_SIZE,
        ha="center",
        va="center",
        arrowprops={"arrowstyle": "-|>", "color": "0.3", "linewidth": 0.5},
    )


def annotate_western_et_al_2024_steps(
    concentration_axes: Sequence[matplotlib.axes.Axes],
    emissions_axes: Sequence[matplotlib.axes.Axes],
    inverse: pd.Series,
    constant: pd.Series,
    extension: pd.Series,
) -> None:
    """
    Annotate the steps of our extension of Western et al. (2024) on one gas' panels

    The steps are numbered, so the caption can refer to them.

    Parameters
    ----------
    concentration_axes
        Axes of the gas' concentration panel

    emissions_axes
        Axes of the gas' emissions panel

    inverse
        Emissions implied by the Western et al. (2024) projection

    constant
        Constant emissions which drive our extension

    extension
        Our extension of the concentrations
    """
    before_break, after_break = emissions_axes
    year = 2060.5
    annotate_step(
        before_break,
        "1. Invert\nWestern et al. (2024)",
        (year, inverse.loc[year]),
        (0.65, 0.6),
    )
    year = 2300.5
    annotate_step(
        after_break,
        "2. Hold 2100\nemissions\nconstant",
        (year, constant.loc[year]),
        (0.5, 0.55),
    )
    _, concentrations_after_break = concentration_axes
    annotate_step(
        concentrations_after_break,
        "3. Run one-box\nmodel forward",
        (year, extension.loc[year]),
        (0.5, 0.55),
    )


@manuscript_style
def generate_magicc_mismatch_figure(
    gas: str, outfile: Path, summary_file: Path, bundle_dir: Path = BUNDLE_DIR
) -> Path:
    """
    Generate the figure of the mismatch between our concentrations and MAGICC's

    For a gas which uses a single concentration projection.
    MAGICC is run with the official emissions, including the official extensions,
    which differ from the emissions implied by our concentrations after 2100.

    Parameters
    ----------
    gas
        Gas of interest

    outfile
        File in which to write the figure

    summary_file
        File which holds the MAGICC summary

    bundle_dir
        Directory which holds the bundle

    Returns
    -------
    :
        `outfile`
    """
    rows = (
        Row(panels=(broken("emissions"),), height=ROW_HEIGHT),
        Row(panels=(broken("concentrations"),), height=ROW_HEIGHT),
        legend_row(),
    )
    titles = {
        "emissions": "Emissions",
        "concentrations": "Global-, annual-mean concentrations",
    }
    fig, axes = create_figure(rows, settings=LAYOUT_SETTINGS)

    def panel(name: str) -> list[matplotlib.axes.Axes]:
        return get_panel_axes(rows, axes, name)

    ours_colour = HISTORY_COLOUR
    emissions_units = get_emissions_units(gas)
    inverse_emissions = get_single_projection_inverse_emissions(
        gas, emissions_units, bundle_dir
    )
    units = get_units(gas)
    ours = get_scenario(gas, SHOWN_SCENARIO, units)
    history = get_history(gas, units)

    emissions_foreground = [inverse_emissions]
    concentrations_foreground = [ours, history]
    for scenario in ALL_SCENARIOS:
        official = get_official_emissions(
            gas, scenario.name, emissions_units, bundle_dir
        )
        magicc = get_magicc_concentrations(summary_file, gas, scenario.name, 0.5, units)
        plot_in_segments(
            official,
            panel("emissions"),
            SCENARIO_SEGMENTS,
            color=scenario.colour,
            linestyle=RAW_LINESTYLE,
        )
        plot_in_segments(
            magicc,
            panel("concentrations"),
            SCENARIO_SEGMENTS,
            color=scenario.colour,
            linestyle=RAW_LINESTYLE,
        )
        emissions_foreground.append(official)
        concentrations_foreground.append(magicc)

    plot_in_segments(
        inverse_emissions,
        panel("emissions"),
        SCENARIO_SEGMENTS,
        color=ours_colour,
        linewidth=HISTORY_LINE_WIDTH,
        zorder=HISTORY_ZORDER,
    )
    plot_in_segments(
        ours,
        panel("concentrations"),
        SCENARIO_SEGMENTS,
        color=ours_colour,
        linewidth=HISTORY_LINE_WIDTH,
        zorder=HISTORY_ZORDER,
    )
    plot_in_segments(
        history,
        panel("concentrations"),
        SCENARIO_SEGMENTS,
        color=ours_colour,
        linewidth=HISTORY_LINE_WIDTH,
        zorder=HISTORY_ZORDER,
    )
    finish_time_axes(
        panel("emissions"),
        SCENARIO_SEGMENTS,
        units_label(emissions_units),
        emissions_foreground,
    )
    finish_time_axes(
        panel("concentrations"),
        SCENARIO_SEGMENTS,
        units_label(units),
        concentrations_foreground,
    )

    ours_handle = line_handle(ours_colour, linewidth=HISTORY_LINE_WIDTH)
    add_panel_legend(
        panel("emissions")[0],
        [
            ("Inverse emissions of the CMIP7 forcings", ours_handle),
            ("Official emissions", line_handle("0.4", linestyle=RAW_LINESTYLE)),
        ],
    )
    add_panel_legend(
        panel("concentrations")[0],
        [
            (f"{HISTORY_LABEL} and CMIP7 forcings", ours_handle),
            ("MAGICC median", line_handle("0.4", linestyle=RAW_LINESTYLE)),
        ],
    )
    add_legend(
        axes[LEGEND_PANEL],
        [[(f"Datasets (dashed lines, {ALL_PANELS})", scenario_entries())]],
    )

    return save(fig, axes, rows, titles, outfile)
