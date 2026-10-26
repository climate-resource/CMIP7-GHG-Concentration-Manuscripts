"""
Generation of the scenario manuscript's methods figures

The figures follow the plans in the manuscript's `methods.tex`:

- time axes run to 2500 and are broken a little after 2100
- harmonisation and one-box vs. MAGICC panels show all the marker scenarios
- everything else shows one scenario, L
  (the gases which use a single concentration projection
  are the same in every scenario anyway)
- latitudinal gradient (and seasonality change) PCs are shown over time,
  the regressions behind them are in their own (appendix) figure

Built on the historical manuscript's layout and broken axis helpers,
so all the manuscript's figures read the same way.
"""

from __future__ import annotations

import re
from collections.abc import Sequence
from pathlib import Path

import matplotlib.axes
import matplotlib.lines
import matplotlib.patches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

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
    add_colour_bar,
    get_ticks_with_ends,
    label_name,
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
    NPP_VARIABLE,
)
from local.scenario_ghg_forcing_for_cmip7.methods_data import (
    CO2_SEASONALITY_REGRESSION,
    LAST_HISTORICAL_YEAR,
    get_historical_pcs,
    get_history,
    get_magicc,
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

MARKER_SCENARIOS = tuple(s for s in CMIP7_SCENARIOS if s.name != "vl-cf")
"""Marker scenarios, lowest to highest

VL-CF is left out of the harmonisation and one-box panels:
it is harmonised in 2015, rather than 2022, so isn't comparable.
"""

ALL_MAGICC_SCENARIOS = CMIP7_SCENARIOS
"""Scenarios MAGICC was run for, lowest to highest (including VL-CF)"""

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

ENSEMBLE_SEGMENTS = ((1990.0, EXTENSIONS_START_YEAR),)
"""Time axis for the MAGICC ensemble panels (not broken)"""

SEGMENT_SHARES = (0.6, 0.4)
"""Share of a broken panel's width given to each piece"""

HARMONISATION_DETAIL_SEGMENTS = (
    (2010.0, 2050.0),
    (2050.0, BREAK_YEAR),
    (BREAK_YEAR, END_YEAR),
)
"""Time axis pieces for harmonisation panels which show the harmonisation in detail

The first piece, where the harmonisation happens, gets half the panel's width
(see [HARMONISATION_DETAIL_SHARES][]).
The breaks are a little after round years (2040 and 2100),
so those years' tick labels sit clear of the break marks.
"""

HARMONISATION_DETAIL_SHARES = (0.5, 0.25, 0.25)
"""Share of the panel's width given to each of [HARMONISATION_DETAIL_SEGMENTS][]"""

EXTENSION_TICKS = (2300.0, 2500.0)
"""Major ticks of the piece after the break"""

HARMONISATION_DETAIL_TICKS = (
    (2010.0, 2020.0, 2030.0, 2040.0),
    (2060.0, 2080.0, 2100.0),
    EXTENSION_TICKS,
)
"""Major ticks of each of [HARMONISATION_DETAIL_SEGMENTS][]"""

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

PC_LINESTYLES = ("-", ":")
"""Line style of each PC (first, second)"""

PC_NAMES = ("1st PC", "2nd PC")
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


def broken(name: str, width: float = 1.0) -> Panel:
    """Get a panel which is a time axis broken at [BREAK_YEAR][]"""
    return Panel(name, width=width, broken=True, broken_split=SEGMENT_SHARES)


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
    scenarios: Sequence = MARKER_SCENARIOS,
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


def harmonisation_panel(name: str = "harmonisation") -> Panel:
    """Get a harmonisation panel (broken into [HARMONISATION_DETAIL_SEGMENTS][])"""
    return Panel(name, broken=True, broken_split=HARMONISATION_DETAIL_SHARES)


def plot_harmonisation(
    axes: Sequence[matplotlib.axes.Axes], gas: str, summary_file: Path
) -> None:
    """
    Plot the raw MAGICC median, the harmonised projection and history, all scenarios

    The harmonised projection is the published data.
    The time axis is broken into [HARMONISATION_DETAIL_SEGMENTS][],
    so the harmonisation itself gets half the panel's width.
    The panel must be made with [harmonisation_panel][].
    """
    segments = HARMONISATION_DETAIL_SEGMENTS
    units = get_units(gas)
    history = get_history(gas, units)
    foreground = [history]
    for scenario in MARKER_SCENARIOS:
        raw = get_magicc_concentrations(summary_file, gas, scenario.name, 0.5, units)
        harmonised = get_scenario(gas, scenario.name, units)
        plot_in_segments(
            raw,
            axes,
            segments,
            color=scenario.colour,
            linestyle=RAW_LINESTYLE,
            linewidth=0.6,
        )
        plot_in_segments(harmonised, axes, segments, color=scenario.colour)
        foreground.extend([raw, harmonised])

    plot_in_segments(
        history,
        axes,
        segments,
        color=HISTORY_COLOUR,
        linewidth=HISTORY_LINE_WIDTH,
        zorder=3.0,
    )
    finish_time_axes(axes, segments, units_label(units), foreground)
    for ax, segment, ticks in zip(axes, segments, HARMONISATION_DETAIL_TICKS):
        ax.set_xticks(ticks)
        # Setting the ticks can move the limits
        ax.set_xlim(segment)


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
    scenario: str,
    segments: Sequence[tuple[float, float]] = LONG_SEGMENTS,
    units: str = "-",
) -> None:
    """
    Plot PCs over time: the historical dataset's and one scenario's projection

    Parameters
    ----------
    axes
        Axes of each piece of the time axis

    historical
        Historical PCs, one column per EOF (`None` if there are none)

    projected
        Projected PCs, with columns `(eof, scenario)` or just `scenario`

    scenario
        Scenario to show (if the projection is the same in every scenario,
        it has only one scenario, which is shown whatever this is)

    segments
        Each piece of the time axis

    units
        Units of the PCs (`"-"` if dimensionless)
    """
    if not isinstance(projected.columns, pd.MultiIndex):
        projected = pd.concat({0: projected}, axis="columns")

    colour = SCENARIOS[SHOWN_SCENARIO].colour
    foreground = []
    for eof in projected.columns.get_level_values(0).unique():
        eof_projected = projected[eof]
        if eof_projected.shape[1] == 1:
            eof_projected = eof_projected.iloc[:, 0]
        else:
            eof_projected = eof_projected[scenario]

        linestyle = PC_LINESTYLES[eof]
        plot_in_segments(
            eof_projected, axes, segments, color=colour, linestyle=linestyle
        )
        foreground.append(eof_projected)
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
        Row(panels=(Panel("ensemble"),), height=ROW_HEIGHT),
        Row(panels=(harmonisation_panel(),), height=ROW_HEIGHT),
        Row(panels=(broken("npp"), broken("seasonality-pc")), height=ROW_HEIGHT),
        Row(panels=(broken("lat-gradient-pcs"),), height=ROW_HEIGHT),
        legend_row(height=0.55),
    )
    titles = {
        "ensemble": "MAGICC ensemble (raw), L",
        "harmonisation": "Harmonisation, all scenarios",
        "npp": "MAGICC median NPP",
        "seasonality-pc": "Seasonality change PC",
        "lat-gradient-pcs": "Latitudinal gradient PCs",
    }
    fig, axes = create_figure(rows, settings=LAYOUT_SETTINGS)

    def panel(name: str) -> list[matplotlib.axes.Axes]:
        return get_panel_axes(rows, axes, name)

    plot_magicc_ensemble(panel("ensemble"), "co2", summary_file)
    plot_harmonisation(panel("harmonisation"), "co2", summary_file)

    # NPP is the same in every scenario up to the end of the history
    npp = get_magicc(summary_file, NPP_VARIABLE, SHOWN_SCENARIO, 0.5)
    npp_units = str(npp.name)
    npp_history = npp[npp.index < LAST_HISTORICAL_YEAR + 1]
    npp_projected = npp[npp.index >= LAST_HISTORICAL_YEAR]
    plot_in_segments(
        npp_projected,
        panel("npp"),
        LONG_SEGMENTS,
        color=SCENARIOS[SHOWN_SCENARIO].colour,
    )
    plot_in_segments(
        npp_history,
        panel("npp"),
        LONG_SEGMENTS,
        color=HISTORY_COLOUR,
        linewidth=HISTORY_LINE_WIDTH,
    )
    finish_time_axes(
        panel("npp"),
        LONG_SEGMENTS,
        units_label(npp_units),
        [npp_history, npp_projected],
    )

    plot_pcs(
        panel("seasonality-pc"),
        get_historical_pcs("co2", kind="seasonality-change"),
        get_projected_co2_seasonality_change_pc(bundle_dir),
        SHOWN_SCENARIO,
    )
    plot_pcs(
        panel("lat-gradient-pcs"),
        get_historical_pcs("co2"),
        get_projected_pcs("co2", bundle_dir),
        SHOWN_SCENARIO,
    )

    add_legend(
        axes[LEGEND_PANEL],
        [
            (HISTORY_LABEL, history_handle()),
            *ensemble_entries(),
            *harmonisation_entries(),
            *pc_entries(2),
        ],
        scenarios=MARKER_SCENARIOS,
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
        Row(panels=(harmonisation_panel(),), height=ROW_HEIGHT),
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
        "harmonisation": "Harmonisation, all scenarios",
        "lat-gradient-pcs": "Latitudinal gradient PCs",
    }
    fig, axes = create_figure(rows, settings=LAYOUT_SETTINGS)

    def panel(name: str) -> list[matplotlib.axes.Axes]:
        return get_panel_axes(rows, axes, name)

    units = get_units("ch4")

    # MAGICC versions over time
    foreground = []
    for scenario in ALL_MAGICC_SCENARIOS:
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
    n_scenarios = len(ALL_MAGICC_SCENARIOS)
    bar_width = 0.8 / (2 * n_scenarios)
    highest = -np.inf
    for i, year in enumerate(CH4_BAR_YEARS):
        for j, scenario in enumerate(ALL_MAGICC_SCENARIOS):
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

    plot_harmonisation(panel("harmonisation"), "ch4", summary_file)
    plot_pcs(
        panel("lat-gradient-pcs"),
        get_historical_pcs("ch4"),
        get_projected_pcs("ch4", bundle_dir),
        SHOWN_SCENARIO,
    )

    add_legend(
        axes[LEGEND_PANEL],
        [
            (HISTORY_LABEL, history_handle()),
            *harmonisation_entries(),
            *pc_entries(2),
        ],
        scenarios=ALL_MAGICC_SCENARIOS,
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
            panels=(broken("lat-gradient-pcs", width=1.25), Panel("ensemble")),
            height=ROW_HEIGHT,
        ),
        Row(panels=(harmonisation_panel(),), height=ROW_HEIGHT),
        legend_row(height=0.55),
    )
    projected_pcs = get_projected_pcs(gas, bundle_dir)
    n_eofs = projected_pcs.columns.get_level_values(0).nunique()
    titles = {
        "ensemble": "MAGICC ensemble (raw), L",
        "harmonisation": "Harmonisation, all scenarios",
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
        SHOWN_SCENARIO,
    )
    plot_magicc_ensemble(panel("ensemble"), gas, summary_file)
    plot_harmonisation(panel("harmonisation"), gas, summary_file)

    add_legend(
        axes[LEGEND_PANEL],
        [
            (HISTORY_LABEL, history_handle()),
            *pc_entries(n_eofs),
            *ensemble_entries(),
            *harmonisation_entries(),
        ],
        scenarios=MARKER_SCENARIOS,
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
        SHOWN_SCENARIO,
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
    for scenario in MARKER_SCENARIOS:
        one_box = get_scenario(gas, scenario.name, units)
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
            one_box, panel("one-box"), ONE_BOX_SEGMENTS, color=scenario.colour
        )
        foreground.extend([one_box, magicc])

    plot_in_segments(
        history,
        panel("one-box"),
        ONE_BOX_SEGMENTS,
        color=HISTORY_COLOUR,
        linewidth=HISTORY_LINE_WIDTH,
        zorder=3.0,
    )
    finish_time_axes(panel("one-box"), ONE_BOX_SEGMENTS, units_label(units), foreground)

    plot_pcs(
        panel("lat-gradient-pcs"),
        get_historical_pcs(gas),
        get_projected_pcs(gas, bundle_dir),
        SHOWN_SCENARIO,
    )

    add_legend(
        axes[LEGEND_PANEL],
        [
            (HISTORY_LABEL, history_handle()),
            ("One-box model (output)", line_handle("0.4")),
            ("MAGICC median", line_handle("0.4", linestyle=RAW_LINESTYLE)),
        ],
        scenarios=MARKER_SCENARIOS,
    )

    return save(fig, axes, rows, titles, outfile)


REGRESSION_PANELS = {
    "co2": ("co2", "CO2, first lat. gradient PC"),
    "ch4": ("ch4", "CH4, first lat. gradient PC"),
    "n2o": ("n2o", "N2O, first lat. gradient PC"),
    "co2-seasonality": (CO2_SEASONALITY_REGRESSION, "CO2, seasonality change PC"),
}
"""Regression shown in each panel of the regressions figure, and its title"""


@manuscript_style
def generate_regressions_figure(outfile: Path, bundle_dir: Path = BUNDLE_DIR) -> Path:
    """
    Generate the figure of the regressions we re-derived

    Parameters
    ----------
    outfile
        File in which to write the figure

    bundle_dir
        Directory which holds the bundle (including our re-runs)

    Returns
    -------
    :
        `outfile`
    """
    names = list(REGRESSION_PANELS)
    rows = tuple(
        Row(
            panels=(Panel(left, colour_bar=True), Panel(right, colour_bar=True)),
            height=1.3,
        )
        for left, right in zip(names[::2], names[1::2])
    )
    titles = {name: label_name(title) for name, (_, title) in REGRESSION_PANELS.items()}
    fig, axes = create_figure(rows, settings=LAYOUT_SETTINGS)

    for name, (key, _) in REGRESSION_PANELS.items():
        ax = axes[name]
        data, info = get_regression(key, bundle_dir)
        scatter = ax.scatter(
            data["x"], data["y"], c=data["year"], s=3.0, cmap="viridis", linewidths=0.0
        )
        x = np.linspace(data["x"].min(), data["x"].max(), 2)
        ax.plot(
            x, info["slope"] * x + info["intercept"], color="k", linewidth=LINE_WIDTH
        )
        ax.set_xlabel(
            label_name(f"[{get_short_units(str(info['x_unit']))}]"), fontsize="small"
        )
        ax.set_ylabel("PC value [-]", fontsize="small")
        ax.tick_params(labelsize="small")
        add_colour_bar(fig, scatter, cax=axes[f"{name}-colour-bar"], label="Year")

    return save_without_legend(fig, axes, rows, titles, outfile)


UNIT_ABBREVIATIONS = {
    "gigametric_ton": "Gt",
    "megametric_ton": "Mt",
    "kilometric_ton": "kt",
    "metric_ton": "t",
    "gigatC": "GtC",
}
"""Short form of the units pint writes out in full"""


def get_short_units(units: str) -> str:
    """
    Get units as written for a reader

    For example, "Mt CO2 / yr" for "CO2 * megametric_ton / yr".

    Parameters
    ----------
    units
        Units, as pint writes them

    Returns
    -------
    :
        Units, as a reader expects them
    """
    for long, short in UNIT_ABBREVIATIONS.items():
        units = re.sub(rf"\b{long}\b", short, units)

    # pint puts the species first
    match = re.fullmatch(r"(\w+) \* (\w+) / (\w+)", units)
    if match:
        species, mass, time = match.groups()
        units = f"{mass} {species} / {time}"

    return units


def save_without_legend(
    fig: plt.Figure,
    axes: dict[str, matplotlib.axes.Axes],
    rows: Sequence[Row],
    titles: dict[str, str],
    outfile: Path,
) -> Path:
    """Label, lay out and save a figure which has no legend panel, then close it"""
    label_panels(rows, axes, titles)
    lay_out_figure(fig, axes, rows, settings=LAYOUT_SETTINGS)

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
    for scenario in MARKER_SCENARIOS:
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
        zorder=3.0,
    )
    plot_in_segments(
        ours,
        panel("concentrations"),
        SCENARIO_SEGMENTS,
        color=ours_colour,
        linewidth=HISTORY_LINE_WIDTH,
        zorder=3.0,
    )
    plot_in_segments(
        history,
        panel("concentrations"),
        SCENARIO_SEGMENTS,
        color=ours_colour,
        linewidth=HISTORY_LINE_WIDTH,
        zorder=3.0,
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
        scenarios=MARKER_SCENARIOS,
    )

    return save(fig, axes, rows, titles, outfile)
