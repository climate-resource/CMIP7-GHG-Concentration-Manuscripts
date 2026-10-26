"""
Generation of the scenario manuscript's methods figures

The figures follow the plans in the manuscript's `methods.tex`:

- time axes run to 2500 and are broken a little after 2100
- panels which compare scenarios show all of them, including VL-CF
- everything else shows one scenario, L
  (the gases which use a single concentration projection
  are the same in every scenario anyway)
- latitudinal gradient (and seasonality change) PCs are shown over time,
  the regressions behind them are in their own (appendix) figure

Built on the historical manuscript's layout and broken axis helpers,
so all the manuscript's figures read the same way.
"""

from __future__ import annotations

import dataclasses
from collections.abc import Sequence
from pathlib import Path

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

HARMONISATION_PANEL_WIDTHS = (0.8, 1.2)
"""Share of the row's width given to the harmonisation and projection panels"""

EXTENSION_TICKS = (2300.0, 2500.0)
"""Major ticks of the piece after the break"""

EXTENSION_MINOR_TICKS = (2200.0, 2400.0)
"""Minor ticks of the piece after the break"""

CH4_BAR_YEARS = (2050, 2100, 2300, 2500)
"""Years shown in the box plot comparing MAGICC versions for CH4"""

CH4_BOX_Y_ROUNDING = 100.0
"""The CH4 box plot's axis starts at the pre-industrial value, rounded down to this"""

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

MAX_LEGEND_ENTRIES_PER_LINE = 4
"""Most legend entries which fit on one line across a figure"""

ROW_HEIGHT = 0.95
"""Height of a row of time series panels, in inches"""

LEGEND_ROW_HEIGHT = 0.35
"""Height of the legend row, in inches"""

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


def legend_row(height: float = LEGEND_ROW_HEIGHT) -> Row:
    """Get the row which holds a figure's legend"""
    return Row(panels=(Panel(LEGEND_PANEL),), height=height)


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


def add_legend(
    ax: matplotlib.axes.Axes,
    entries: Sequence[tuple[str, matplotlib.artist.Artist]],
    scenarios: Sequence = (),
) -> None:
    """
    Put a figure's legend in its legend panel

    The entries which explain line styles go on one line,
    the scenarios' colours (if any) on the line below.

    Parameters
    ----------
    ax
        The legend panel's axes

    entries
        Label and handle of each entry which explains a line style, in order

    scenarios
        Scenarios whose colours to explain, in order
    """
    entries = list(entries)
    # As few lines as fit, with the entries shared evenly between them
    n_lines = int(np.ceil(len(entries) / MAX_LEGEND_ENTRIES_PER_LINE))
    per_line = int(np.ceil(len(entries) / max(n_lines, 1)))
    lines = [entries[i : i + per_line] for i in range(0, len(entries), per_line)]
    if scenarios:
        lines.append(scenario_entries(scenarios))

    ax.axis("off")
    for i, line in enumerate(lines):
        legend = ax.legend(
            handles=[handle for _, handle in line],
            labels=[label for label, _ in line],
            loc="center",
            # Lines spaced evenly from top to bottom
            bbox_to_anchor=(0.5, 1.0 - (i + 0.5) / len(lines)),
            ncols=len(line),
            frameon=False,
            fontsize=LEGEND_FONT_SIZE,
            handlelength=1.8,
            columnspacing=0.8,
            borderaxespad=0.0,
        )
        # See `add_scenario_legend` in the results figure for why
        legend.set_in_layout(False)
        ax.add_artist(legend)


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


def ensemble_entries() -> list[tuple[str, matplotlib.artist.Artist]]:
    """Get the legend entries for the MAGICC ensemble panel"""
    colour = SCENARIOS[SHOWN_SCENARIO].colour
    return [
        ("MAGICC median, L", line_handle(colour)),
        (
            "MAGICC 5-95%, L",
            matplotlib.patches.Patch(facecolor=colour, alpha=0.25, linewidth=0.0),
        ),
    ]


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
    "harmonisation": "Harmonisation, all scenarios, 2010\u20132035",
    "projections": "Harmonisation, all scenarios, 2015\u20132500",
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


def harmonisation_entries() -> list[tuple[str, matplotlib.artist.Artist]]:
    """Get the legend entries for the line styles of the harmonisation panel"""
    return [
        ("Raw MAGICC median", line_handle("0.4", linestyle=RAW_LINESTYLE)),
        ("Harmonised (output)", line_handle("0.4")),
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


def pc_entries(n_eofs: int) -> list[tuple[str, matplotlib.artist.Artist]]:
    """Get the legend entries for the line styles of the PC panels"""
    if n_eofs == 1:
        return []

    return [
        (name, line_handle("0.4", linestyle=linestyle))
        for name, linestyle in zip(PC_NAMES[:n_eofs], PC_LINESTYLES)
    ]


def history_and_l_entries(
    projection_label: str = "L",
) -> list[tuple[str, matplotlib.artist.Artist]]:
    """Get the legend entries for panels which show history and L"""
    return [
        (HISTORY_LABEL, history_handle()),
        (projection_label, line_handle(SCENARIOS[SHOWN_SCENARIO].colour)),
    ]


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
        legend_row(height=0.55),
    )
    titles = {
        "ensemble": "MAGICC ensemble (raw), L",
        **HARMONISATION_TITLES,
        "seasonality-pc": "Seasonality change PC, all scenarios",
        "lat-gradient-pcs": "Latitudinal gradient PCs, all scenarios",
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

    add_legend(
        axes[LEGEND_PANEL],
        [
            (HISTORY_LABEL, history_handle()),
            *ensemble_entries(),
            *harmonisation_entries(),
            *pc_entries(2),
        ],
        scenarios=ALL_SCENARIOS,
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
        Row(panels=(Panel("versions-bars"),), height=ROW_HEIGHT),
        Row(panels=harmonisation_panels(), height=ROW_HEIGHT),
        Row(panels=(broken("lat-gradient-pcs"),), height=ROW_HEIGHT),
        legend_row(height=0.55),
    )
    titles = {
        "versions": (
            f"Raw MAGICC median, {MAGICC_VERSION} (solid), "
            f"{MAGICC_AR6_VERSION} (dashed)"
        ),
        "versions-bars": (
            "Raw MAGICC median (line) and 5-95% range (box)\n"
            f"{MAGICC_AR6_VERSION} (lighter, left), {MAGICC_VERSION} (right)"
        ),
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
    (ax,) = panel("versions-bars")
    n_scenarios = len(ALL_SCENARIOS)
    bar_width = 0.8 / (2 * n_scenarios)
    highest = -np.inf
    for i, year in enumerate(CH4_BAR_YEARS):
        for j, scenario in enumerate(ALL_SCENARIOS):
            for k, (climate_model, alpha) in enumerate(
                ((MAGICC_AR6_VERSION, 0.45), (MAGICC_VERSION, 1.0))
            ):
                median, lower, upper = (
                    get_magicc_concentrations(
                        summary_file, "ch4", scenario.name, q, units, climate_model
                    ).loc[year + 0.5]
                    for q in (0.5, *ENSEMBLE_RANGE)
                )
                x = i - 0.4 + (2 * j + k + 0.5) * bar_width
                # A box from the 5th to the 95th percentile, with the median across it
                ax.bar(
                    x,
                    upper - lower,
                    bottom=lower,
                    width=bar_width,
                    color=scenario.colour,
                    alpha=alpha,
                    linewidth=0.0,
                )
                ax.hlines(
                    median,
                    x - bar_width / 2,
                    x + bar_width / 2,
                    color="0.1",
                    linewidth=0.6,
                )
                highest = max(highest, upper)

    # Start at (roughly) the pre-industrial concentration, rather than zero,
    # so the boxes get the panel's height
    pre_industrial = get_history("ch4", units).iloc[0]
    ax.set_ylim(
        np.floor(pre_industrial / CH4_BOX_Y_ROUNDING) * CH4_BOX_Y_ROUNDING,
        highest + 0.03 * (highest - pre_industrial),
    )
    ax.set_xticks(range(len(CH4_BAR_YEARS)))
    ax.set_xticklabels([str(y) for y in CH4_BAR_YEARS])
    ax.tick_params(labelsize="small")
    ax.set_ylabel(units_label(units), fontsize="small")

    plot_harmonisation(
        panel("harmonisation"), panel("projections"), "ch4", summary_file
    )
    plot_pcs(
        panel("lat-gradient-pcs"),
        get_historical_pcs("ch4"),
        get_projected_pcs("ch4", bundle_dir),
    )

    add_legend(
        axes[LEGEND_PANEL],
        [
            (HISTORY_LABEL, history_handle()),
            *harmonisation_entries(),
            *pc_entries(2),
        ],
        scenarios=ALL_SCENARIOS,
    )

    return save(fig, axes, rows, titles, outfile)


@manuscript_style
def generate_n2o_like_figure(
    gas: str, outfile: Path, summary_file: Path, bundle_dir: Path = BUNDLE_DIR
) -> Path:
    """
    Generate the methods figure of a gas processed like N2O

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
            panels=(
                broken("lat-gradient-pcs"),
                broken("ensemble", shares=ENSEMBLE_SEGMENT_SHARES),
            ),
            height=ROW_HEIGHT,
        ),
        Row(panels=harmonisation_panels(), height=ROW_HEIGHT),
        legend_row(height=0.55),
    )
    projected_pcs = get_projected_pcs(gas, bundle_dir)
    n_eofs = projected_pcs.columns.get_level_values(0).nunique()
    titles = {
        "ensemble": "MAGICC ensemble (raw), L",
        **HARMONISATION_TITLES,
        "lat-gradient-pcs": "Latitudinal gradient PCs"
        if n_eofs > 1
        else "Latitudinal gradient PC",
    }
    fig, axes = create_figure(rows, settings=LAYOUT_SETTINGS)

    def panel(name: str) -> list[matplotlib.axes.Axes]:
        return get_panel_axes(rows, axes, name)

    plot_pcs(
        panel("lat-gradient-pcs"),
        get_historical_pcs(gas),
        projected_pcs,
    )
    plot_magicc_ensemble(panel("ensemble"), gas, summary_file)
    plot_harmonisation(panel("harmonisation"), panel("projections"), gas, summary_file)

    add_legend(
        axes[LEGEND_PANEL],
        [
            (HISTORY_LABEL, history_handle()),
            *pc_entries(n_eofs),
            *ensemble_entries(),
            *harmonisation_entries(),
        ],
        scenarios=ALL_SCENARIOS,
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
        Row(
            panels=(broken("projection"), broken("inverse-emissions")),
            height=ROW_HEIGHT,
        ),
        Row(panels=(broken("lat-gradient-pcs"),), height=ROW_HEIGHT),
        legend_row(),
    )
    titles = {
        "projection": "Global-, annual-mean",
        "inverse-emissions": "Inverse emissions",
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

    emissions_units = get_emissions_units(gas)
    inverse_emissions = get_single_projection_inverse_emissions(
        gas, emissions_units, bundle_dir
    )
    shade_western_et_al_2024_extension(
        gas, panel("inverse-emissions"), SCENARIO_SEGMENTS
    )
    plot_in_segments(
        inverse_emissions, panel("inverse-emissions"), SCENARIO_SEGMENTS, color=colour
    )
    finish_time_axes(
        panel("inverse-emissions"),
        SCENARIO_SEGMENTS,
        units_label(emissions_units),
        [inverse_emissions],
    )

    plot_pcs(
        panel("lat-gradient-pcs"),
        get_historical_pcs(gas),
        get_projected_pcs(gas, bundle_dir),
    )

    add_legend(
        axes[LEGEND_PANEL],
        history_and_l_entries("Projection (the same in all scenarios)"),
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
        legend_row(height=0.5),
    )
    titles = {
        "one-box": (
            "One-box model (our output) and MAGICC median (not used), all scenarios"
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
        get_projected_pcs(gas, bundle_dir),
    )

    add_legend(
        axes[LEGEND_PANEL],
        [
            (HISTORY_LABEL, history_handle()),
            ("One-box model (output)", line_handle("0.4")),
            ("MAGICC median", line_handle("0.4", linestyle=RAW_LINESTYLE)),
        ],
        scenarios=ALL_SCENARIOS,
    )

    return save(fig, axes, rows, titles, outfile)


REGRESSION_PANELS: dict[str, dict[str, tuple[str, str, str]]] = {
    "co2": {
        "seasonality-change": (
            CO2_SEASONALITY_REGRESSION,
            "Seasonality change PC\nvs. NPP",
            "MAGICC median NPP",
        ),
        "lat-gradient": (
            "co2",
            "Lat. gradient PC0\nvs. emissions",
            "Fossil co2 emissions",
        ),
    },
    "ch4": {
        "lat-gradient": (
            "ch4",
            "Lat. gradient PC0\nvs. emissions",
            "Fossil ch4 emissions",
        ),
    },
    "n2o": {
        "lat-gradient": (
            "n2o",
            "Lat. gradient PC0\nvs. emissions",
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
            Row(panels=(broken(gas),), height=ROW_HEIGHT)
            for gas in WESTERN_ET_AL_2024_GASES
        ),
        legend_row(),
    )
    titles = {gas: label_name(gas) for gas in WESTERN_ET_AL_2024_GASES}
    fig, axes = create_figure(rows, settings=LAYOUT_SETTINGS)

    wmo_colour = "0.55"
    western_colour = SCENARIOS[SHOWN_SCENARIO].colour
    for gas in WESTERN_ET_AL_2024_GASES:
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

    add_legend(
        axes[LEGEND_PANEL],
        [
            (HISTORY_LABEL, history_handle()),
            ("Western et al. (2024)", line_handle(western_colour)),
            ("Our extension", line_handle(western_colour, linestyle=RAW_LINESTYLE)),
            ("WMO (2022)", line_handle(wmo_colour)),
        ],
    )

    return save(fig, axes, rows, titles, outfile)


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
        legend_row(height=0.5),
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

    add_legend(
        axes[LEGEND_PANEL],
        [
            (
                "Ours (inverse emissions; historical dataset and our concentrations)",
                line_handle(ours_colour, linewidth=HISTORY_LINE_WIDTH),
            ),
            (
                "Official emissions; MAGICC median",
                line_handle("0.4", linestyle=RAW_LINESTYLE),
            ),
        ],
        scenarios=ALL_SCENARIOS,
    )

    return save(fig, axes, rows, titles, outfile)
