"""
Layout of the historical GHG manuscript's methods figures

The methods figures have panels whose shape is not up for negotiation
(the maps have a fixed aspect ratio, the flying carpet is square)
sat alongside panels which can be any shape at all.
Matplotlib's constrained layout engine ties every column of the figure
to every other one, which the fixed panels then over-constrain.

The figures are drawn at the size they are printed at
([FIGURE_WIDTH][], included in the latex at the same width),
so that a font size or a line width here is the one the reader sees.
Drawn any bigger and shrunk to fit the page,
every point size in the figure is shrunk with it,
and text which reads fine on screen comes out unreadable in print.
[manuscript_style][] sets the sizes everything is drawn at to match.

So we do the layout ourselves, in inches, one row at a time.
Each row is independent of the others,
so a row can have as many panels as we like, of whatever widths we like.
Within a row, the fixed-shape panels are sized first
and the flexible panels share out whatever width is left.

The only thing we need matplotlib for is to tell us
how much room each panel's decorations
(tick labels, axis labels, titles, colour bars) take up.
We measure that, lay the figure out again with the room it asked for,
and repeat until the measurements stop moving.
"""

from __future__ import annotations

import functools
import itertools
import string
from collections.abc import Callable, Collection, Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import ParamSpec, TypeVar

import matplotlib.axes
import matplotlib.axis
import matplotlib.figure
import matplotlib.pyplot as plt
import matplotlib.text
import matplotlib.ticker
from loguru import logger
from matplotlib.transforms import Bbox

P = ParamSpec("P")
T = TypeVar("T")

CM = 1.0 / 2.54
"""One centimetre, in inches"""

FIGURE_WIDTH = 12.0 * CM
"""Width of every figure, in inches

The width Copernicus' template asks full-width figures to be included at
(`\\includegraphics[width=12cm]{...}`).
The latex includes the figures at this width too,
so they are printed at exactly the size they are drawn at.
"""

TITLE_FONT_SIZE = 6.5
"""Size of the panels' titles, in points

Set outright rather than relative to the base size:
a panel is only a few centimetres wide,
and a title much bigger than this is wider than its panel.
"""

LEGEND_FONT_SIZE = 5.5
"""Size of the text in legends, in points

Set outright rather than relative to the base size,
because it is already about as small as print can go:
shrinking the axis labels to make room should not shrink the legends too.
"""

RC_PARAMS = {
    # Axis and tick labels are "small", which with this base size is 6 pt
    "font.size": 7.2,
    "legend.fontsize": LEGEND_FONT_SIZE,
    "legend.title_fontsize": LEGEND_FONT_SIZE,
    "axes.linewidth": 0.6,
    "axes.titlepad": 3.0,
    "axes.labelpad": 2.0,
    "lines.linewidth": 1.0,
    "lines.markersize": 3.0,
    "patch.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.minor.width": 0.4,
    "ytick.minor.width": 0.4,
    "xtick.major.size": 2.5,
    "ytick.major.size": 2.5,
    "xtick.minor.size": 1.5,
    "ytick.minor.size": 1.5,
    "xtick.major.pad": 2.0,
    "ytick.major.pad": 2.0,
    "legend.handlelength": 1.5,
    "grid.linewidth": 0.5,
    "savefig.dpi": 300,
}
"""Matplotlib settings for figures drawn at the size they are printed at

Matplotlib's defaults are for a figure of about 16 by 12 cm seen on screen.
Ours are 12 cm wide with a dozen panels in them,
so everything drawn in points (text, lines, markers, ticks)
is scaled down to match.
"""


def manuscript_style(func: Callable[P, T]) -> Callable[P, T]:
    """
    Draw everything `func` draws with [RC_PARAMS][]

    Font sizes are fixed when a piece of text is created,
    so the settings have to be in place for all of the drawing,
    not only when the figure is created or saved.

    Parameters
    ----------
    func
        Function which draws a figure

    Returns
    -------
    :
        `func`, drawing with [RC_PARAMS][]
    """

    @functools.wraps(func)
    def wrapper(*args: P.args, **kwargs: P.kwargs) -> T:
        with plt.rc_context(RC_PARAMS):
            return func(*args, **kwargs)

    return wrapper


@dataclass(frozen=True)
class Panel:
    """
    A panel in the figure
    """

    name: str
    """Name of the panel, used to look it up in the axes"""

    width: float = 1.0
    """
    Share of the row's flexible width this panel takes

    Ignored if `aspect` is set.
    """

    aspect: float | None = None
    """
    Height of the panel's data box divided by its width, if fixed

    If `None`, the panel takes the height of its row
    and its share of the row's flexible width.
    """

    projection: object = None
    """Projection to create the panel's axes with"""

    broken: bool = False
    """
    Whether the panel is a broken axis

    If `True`, the panel is made of axes side by side with a small gap between.
    If `broken_split` is a float, there are two of them,
    `f"{name}-l"` and `f"{name}-r"`.
    If it is a tuple, there is one per element,
    `f"{name}-0"`, `f"{name}-1"` etc.
    """

    broken_split: float | tuple[float, ...] = 0.5
    """
    Share of a broken panel's width given to its left half

    If a tuple, the share of the panel's width given to each of its pieces,
    left to right, which is how a panel is broken into more than two pieces.

    Ignored if `broken` is `False`.

    The two halves of a broken axis rarely cover spans
    which deserve the same amount of room:
    the left half of an extended timeseries
    is usually a millennium of flat line,
    while the right half is where everything happens.
    Giving the left half less than half the panel
    spends the panel's width where the reader needs it.
    """

    broken_gap: float | None = None
    """
    Gap between the pieces of a broken axis, in inches

    If `None`, [LayoutSettings.broken_gap][] is used.
    Pieces which each carry their own vertical axis need a gap wide enough
    to hold that axis' tick labels, which the default is not.
    """

    colour_bar: bool = False
    """
    Whether the panel has a colour bar beside it

    If `True`, an axes for the colour bar is created as `f"{name}-colour-bar"`.
    """

    colour_bar_height: float = 1.0
    """Height of the colour bar, as a fraction of the panel's height"""


@dataclass(frozen=True)
class Row:
    """
    A row of panels in the figure
    """

    panels: tuple[Panel, ...]
    """Panels in the row, left to right"""

    height: float | None = None
    """
    Height of the row's data boxes, in inches

    If `None`, every panel in the row must have a fixed aspect ratio
    and the row is as tall as the panels are when they fill the row's width.
    """


@dataclass
class Pads:
    """
    Room a panel's decorations take up around its data box, in inches
    """

    left: float = 0.35
    right: float = 0.1
    top: float = 0.15
    bottom: float = 0.25

    title_right: float = 0.0
    """
    How far the panel's title reaches past the right of its data box

    Kept apart from `right` because a title sits above its panel,
    so it can reach out over the space the panel's own decorations
    (e.g. a colour bar) leave above them without running into anything.
    All it has to keep clear of is the next panel's title,
    which starts where the next panel's left-hand decorations do,
    see [place_titles][].
    """


@dataclass(frozen=True)
class LayoutSettings:
    """
    Spacing settings for the layout, in inches
    """

    width: float = FIGURE_WIDTH
    """Width of the figure"""

    margin: float = 0.02
    """Space around the outside of the figure's decorations"""

    wspace: float = 0.1
    """Space between the decorations of neighbouring panels in a row"""

    hspace: float = 0.06
    """Space between the decorations of neighbouring rows"""

    broken_gap: float = 0.04
    """Gap between the two halves of a broken axis"""

    colour_bar_gap: float = 0.05
    """Gap between a panel and its colour bar"""

    colour_bar_width: float = 0.07
    """Width of a colour bar"""

    title_gap: float = 0.08
    """Smallest space between a title and the title of the next panel in its row"""

    align_right: bool = False
    """
    Whether every row's data should end at the same place on the right

    By default, each row runs all the way to the edge of the figure.
    Rows whose panels share a time axis have to line up at both ends,
    in which case every row stops where the row which needs the most room
    for its right-hand decorations does.
    """


def get_panel_axes_names(panel: Panel) -> tuple[str, ...]:
    """
    Get the names of the axes which make up a panel's data box

    Parameters
    ----------
    panel
        Panel of interest

    Returns
    -------
    :
        Names of the axes, left to right (colour bar not included)
    """
    if panel.broken:
        if isinstance(panel.broken_split, tuple):
            return tuple(f"{panel.name}-{i}" for i in range(len(panel.broken_split)))

        return (f"{panel.name}-l", f"{panel.name}-r")

    return (panel.name,)


def get_broken_shares(panel: Panel) -> tuple[float, ...]:
    """
    Get the share of a broken panel's width which each of its pieces takes

    Parameters
    ----------
    panel
        Panel of interest

    Returns
    -------
    :
        Share of the panel's width each piece takes, left to right
    """
    if isinstance(panel.broken_split, tuple):
        total = sum(panel.broken_split)
        return tuple(share / total for share in panel.broken_split)

    return (panel.broken_split, 1.0 - panel.broken_split)


def get_colour_bar_axes_name(panel: Panel) -> str:
    """
    Get the name of a panel's colour bar axes
    """
    return f"{panel.name}-colour-bar"


def get_titles(ax: matplotlib.axes.Axes) -> tuple[matplotlib.text.Text, ...]:
    """
    Get the titles of an axes which have something in them

    Parameters
    ----------
    ax
        Axes of interest

    Returns
    -------
    :
        The axes' centre, left and right titles, if they have any text
    """
    return tuple(
        title
        for title in (ax.title, ax._left_title, ax._right_title)
        if title.get_text()
    )


def merge_axes(
    *axes: Mapping[str, matplotlib.axes.Axes],
) -> dict[str, matplotlib.axes.Axes]:
    """
    Merge the axes of several figures into one lookup

    For a function which draws panels that are split between figures:
    each panel is drawn into whichever figure it is in
    without the drawing code having to know which that is.

    Parameters
    ----------
    *axes
        Axes of each figure, as returned by [create_figure][]

    Returns
    -------
    :
        Every figure's axes

    Raises
    ------
    AssertionError
        Two figures have an axes with the same name
    """
    res: dict[str, matplotlib.axes.Axes] = {}
    for figure_axes in axes:
        clashes = set(res).intersection(figure_axes)
        if clashes:
            msg = f"Axes names are used in more than one figure: {sorted(clashes)}"
            raise AssertionError(msg)

        res.update(figure_axes)

    return res


def save_figure(  # noqa: PLR0913
    fig: matplotlib.figure.Figure,
    axes: Mapping[str, matplotlib.axes.Axes],
    rows: tuple[Row, ...],
    titles: Mapping[str, str],
    outfile: Path,
    settings: LayoutSettings = LayoutSettings(),
) -> Path:
    """
    Label, lay out and save a figure, then close it

    Parameters
    ----------
    fig
        Figure to save

    axes
        The figure's axes (only this figure's, see [lay_out_figure][])

    rows
        Rows of the figure, top to bottom

    titles
        Title of each panel, see [label_panels][]

    outfile
        File in which to write the figure

    settings
        Layout settings

    Returns
    -------
    :
        `outfile`
    """
    label_panels(rows, axes, titles)
    # Last, because it needs to know how much room everything takes up
    lay_out_figure(fig, axes, rows, settings=settings)

    outfile.parent.mkdir(exist_ok=True, parents=True)
    logger.info(f"Writing {outfile}")
    fig.savefig(outfile)
    plt.close(fig)

    return outfile


def create_figure(
    rows: tuple[Row, ...], settings: LayoutSettings = LayoutSettings()
) -> tuple[matplotlib.figure.Figure, dict[str, matplotlib.axes.Axes]]:
    """
    Create the figure and all its axes

    The axes are not in their final positions,
    see [lay_out_figure][] for that.

    Parameters
    ----------
    rows
        Rows of the figure, top to bottom

    settings
        Layout settings

    Returns
    -------
    :
        Figure and its axes
    """
    # Height is set properly once the layout is done
    fig = plt.figure(figsize=(settings.width, settings.width))
    axes = {}
    for row in rows:
        for panel in row.panels:
            for name in get_panel_axes_names(panel):
                axes[name] = fig.add_axes(
                    (0.0, 0.0, 0.1, 0.1), projection=panel.projection
                )

            if panel.colour_bar:
                axes[get_colour_bar_axes_name(panel)] = fig.add_axes(
                    (0.0, 0.0, 0.01, 0.1)
                )

    if len(axes) != len(fig.axes):
        msg = "Panel names are not unique"
        raise AssertionError(msg)

    return fig, axes


def label_panels(
    rows: tuple[Row, ...],
    axes: Mapping[str, matplotlib.axes.Axes],
    titles: Mapping[str, str],
    unlabelled: Collection[str] = (),
) -> None:
    """
    Give each panel its title, labelled in reading order

    Parameters
    ----------
    rows
        Rows of the figure, top to bottom

    axes
        The figure's axes

    titles
        Title of each panel

    unlabelled
        Panels which get neither a title nor a label

        For panels which hold something other than data, e.g. a legend.
    """
    panels = [
        panel for row in rows for panel in row.panels if panel.name not in unlabelled
    ]
    if set(titles) != {panel.name for panel in panels}:
        msg = f"Titles don't match the panels: {set(titles)=}"
        raise AssertionError(msg)

    for label, panel in zip(string.ascii_lowercase, panels):
        # Titles go on the leftmost axes of the panel
        ax = axes[get_panel_axes_names(panel)[0]]
        ax.set_title(
            f"$\\bf{{({label})}}$ {titles[panel.name]}",
            loc="left",
            fontsize=TITLE_FONT_SIZE,
        )


def _place(  # noqa: PLR0912, PLR0915
    fig: matplotlib.figure.Figure,
    axes: Mapping[str, matplotlib.axes.Axes],
    rows: tuple[Row, ...],
    pads: Mapping[str, Pads],
    settings: LayoutSettings,
) -> dict[str, tuple[float, float, float, float]]:
    """
    Place every axes, given the room each panel's decorations need

    Returns
    -------
    :
        Data box of each panel, as (x0, y0, width, height) in inches,
        with y measured from the bottom of the figure
    """

    def right_pad(panel: Panel) -> float:
        # Nothing comes after the last panel in a row,
        # so its title has to fit inside the figure like everything else
        return max(pads[panel.name].right, pads[panel.name].title_right)

    def space_between(left: Panel, right: Panel) -> float:
        # Whichever needs more room: the decorations of the two panels,
        # or the left panel's title, which has to clear the right one's.
        # The right one's title starts where its left-hand decorations do,
        # see place_titles.
        return max(
            pads[left.name].right + settings.wspace + pads[right.name].left,
            pads[left.name].title_right + settings.title_gap + pads[right.name].left,
        )

    # Every row's data starts at the same place on the left,
    # so the vertical axes line up.
    # On the right, each row runs all the way to the edge,
    # otherwise a colour bar in one row would leave a gap at the end of all the others,
    # unless we've been asked to line the rows up on the right too.
    left_edge = settings.margin + max(pads[row.panels[0].name].left for row in rows)
    shared_right_edge = (
        settings.width
        - settings.margin
        - max(right_pad(row.panels[-1]) for row in rows)
    )

    # Positions from the top first, because we don't know the figure's height yet
    boxes_from_top = {}
    y = settings.margin
    for row in rows:
        panels = row.panels
        if settings.align_right:
            right_edge = shared_right_edge
        else:
            right_edge = settings.width - settings.margin - right_pad(panels[-1])

        between = sum(
            space_between(left, right) for left, right in itertools.pairwise(panels)
        )
        available = right_edge - left_edge - between

        if row.height is None:
            if any(panel.aspect is None for panel in panels):
                msg = "A row without a height can only hold fixed-aspect panels"
                raise AssertionError(msg)

            height = available / sum(1.0 / panel.aspect for panel in panels)

        else:
            height = row.height

        fixed_width = sum(
            height / panel.aspect for panel in panels if panel.aspect is not None
        )
        flexible_width = available - fixed_width
        flexible_weight = sum(panel.width for panel in panels if panel.aspect is None)
        if flexible_weight > 0.0 and flexible_width <= 0.0:
            msg = (
                f"No width left for the flexible panels in {row}. "
                "The panels' decorations need more width than the row has. "
                "This is usually titles (or legends) which are wider than their panel: "
                "a panel's decorations are measured against its data box, "
                "so they take up more room the narrower the panel gets, "
                "until there is nothing left for the panels themselves. "
                "Shorten or wrap the titles, or move a panel to another row."
            )
            raise AssertionError(msg)

        y += max(pads[panel.name].top for panel in panels)
        x = left_edge
        for i, panel in enumerate(panels):
            if panel.aspect is None:
                width = flexible_width * panel.width / flexible_weight
            else:
                width = height / panel.aspect

            boxes_from_top[panel.name] = (x, y, width, height)
            if i < len(panels) - 1:
                x += width + space_between(panel, panels[i + 1])

        y += height + max(pads[panel.name].bottom for panel in panels)
        y += settings.hspace

    figure_height = y - settings.hspace + settings.margin
    fig.set_size_inches(settings.width, figure_height)

    def to_figure_coords(x0, y0_from_bottom, w, h):
        return (
            x0 / settings.width,
            y0_from_bottom / figure_height,
            w / settings.width,
            h / figure_height,
        )

    boxes = {}
    for row in rows:
        for panel in row.panels:
            x0, y_top, width, height = boxes_from_top[panel.name]
            y0 = figure_height - y_top - height
            boxes[panel.name] = (x0, y0, width, height)

            names = get_panel_axes_names(panel)
            if panel.broken:
                gap = (
                    settings.broken_gap
                    if panel.broken_gap is None
                    else panel.broken_gap
                )
                usable = width - gap * (len(names) - 1)
                x_piece = x0
                for name, share in zip(names, get_broken_shares(panel)):
                    piece_width = usable * share
                    axes[name].set_position(
                        to_figure_coords(x_piece, y0, piece_width, height)
                    )
                    x_piece += piece_width + gap
            else:
                axes[names[0]].set_position(to_figure_coords(x0, y0, width, height))

            if panel.colour_bar:
                colour_bar_height = height * panel.colour_bar_height
                axes[get_colour_bar_axes_name(panel)].set_position(
                    to_figure_coords(
                        x0 + width + settings.colour_bar_gap,
                        y0 + (height - colour_bar_height) / 2.0,
                        settings.colour_bar_width,
                        colour_bar_height,
                    )
                )

    return boxes


def place_titles(
    axes: Mapping[str, matplotlib.axes.Axes],
    rows: tuple[Row, ...],
    pads: Mapping[str, Pads],
    boxes: Mapping[str, tuple[float, float, float, float]],
    settings: LayoutSettings,
) -> None:
    """
    Move each panel's title to the left of its panel

    Matplotlib lines a left-hand title up with the left of the data box,
    which puts the panel's label somewhere in the middle of the panel
    once its vertical axis labels are counted.
    Lined up with the left of the panel's decorations instead,
    the labels sit at the panels' top-left corners,
    and the first panel in every row has its label hard up against
    the left of the figure, so a column of them reads straight down.

    The pieces of a broken panel which have titles of their own
    (they have their own vertical scale)
    have them moved to the start of the gap before them,
    which is where their decorations start.

    Parameters
    ----------
    axes
        The figure's axes

    rows
        Rows of the figure, top to bottom

    pads
        Room each panel's decorations take up

    boxes
        Data box of each panel, as returned by [_place][]

    settings
        Layout settings
    """
    for row in rows:
        for i, panel in enumerate(row.panels):
            x0, _, _, _ = boxes[panel.name]
            title_left = settings.margin if i == 0 else x0 - pads[panel.name].left

            names = get_panel_axes_names(panel)
            gap = settings.broken_gap if panel.broken_gap is None else panel.broken_gap
            for j, name in enumerate(names):
                ax = axes[name]
                # The axes' position, in inches
                fig_width = ax.get_figure().get_size_inches()[0]
                position = ax.get_position(original=True)
                ax_x0 = position.x0 * fig_width
                ax_width = position.width * fig_width
                if j > 0:
                    title_left = ax_x0 - gap

                ax._left_title.set_x((title_left - ax_x0) / ax_width)


def _measure(
    fig: matplotlib.figure.Figure,
    axes: Mapping[str, matplotlib.axes.Axes],
    rows: tuple[Row, ...],
    boxes: Mapping[str, tuple[float, float, float, float]],
) -> dict[str, Pads]:
    """
    Measure the room each panel's decorations take up around its data box
    """
    fig.draw_without_rendering()
    renderer = fig.canvas.get_renderer()

    pads = {}
    for row in rows:
        for panel in row.panels:
            data_names = get_panel_axes_names(panel)
            names = list(data_names)
            if panel.colour_bar:
                names.append(get_colour_bar_axes_name(panel))

            # How far the titles reach to the right is measured on its own,
            # see Pads.title_right. `for_layout_only` is matplotlib's way
            # of leaving a title's width out of the box (but not its height).
            # Only the panel's own titles: a colour bar's label can be drawn
            # as its title, but it is the colour bar's, not the panel's,
            # so it is measured in full like the rest of the colour bar.
            titles = [
                title
                for name in data_names
                for title in get_titles(axes[name])
                if title.get_visible()
            ]

            # An axes which has been hidden has no box, and takes no room
            tight_boxes = [
                box
                for name in names
                if (
                    box := axes[name].get_tightbbox(
                        renderer, for_layout_only=name in data_names
                    )
                )
                is not None
            ]

            # `for_layout_only` also squashes the axis labels to a sliver,
            # which is right for a label under the middle of a row
            # (it can reach over its neighbours' space),
            # but not for one at the end of a row, which then runs off the figure.
            # So they go in at their full size.
            tight_boxes.extend(
                axis.label.get_window_extent(renderer)
                for name in names
                for axis in (axes[name].xaxis, axes[name].yaxis)
                if axis.get_visible()
                and axis.label.get_visible()
                and axis.label.get_text()
            )

            # A 3D axes' box leaves out its axes' labels,
            # which hang out below it, so they are added by hand
            tight_boxes.extend(
                box
                for name in names
                if axes[name].name == "3d"
                for axis in (axes[name].xaxis, axes[name].yaxis, axes[name].zaxis)
                if axis.get_visible()
                and (box := axis.get_tightbbox(renderer)) is not None
            )
            to_inches = fig.dpi_scale_trans.inverted()
            tight = Bbox.union(tight_boxes).transformed(to_inches)

            x0, y0, width, height = boxes[panel.name]
            pads[panel.name] = Pads(
                left=max(x0 - tight.x0, 0.0),
                right=max(tight.x1 - (x0 + width), 0.0),
                top=max(tight.y1 - (y0 + height), 0.0),
                bottom=max(y0 - tight.y0, 0.0),
            )
            if titles:
                title_box = Bbox.union(
                    [title.get_window_extent(renderer) for title in titles]
                ).transformed(to_inches)
                pads[panel.name].top = max(
                    pads[panel.name].top, title_box.y1 - (y0 + height)
                )
                pads[panel.name].title_right = max(title_box.x1 - (x0 + width), 0.0)

    return pads


def _get_tick_labels_in_view(
    axis: matplotlib.axis.Axis,
) -> tuple[list[float], list[matplotlib.text.Text]]:
    """
    Get the major ticks an axis draws, and their labels

    Parameters
    ----------
    axis
        Axis of interest

    Returns
    -------
    :
        Location of each tick in view with a label, and that label
    """
    lower, upper = sorted(axis.get_view_interval())
    tolerance = 1e-9 * max(upper - lower, 1.0)
    in_view = [
        (loc, label)
        for loc, label in zip(axis.get_majorticklocs(), axis.get_majorticklabels())
        if lower - tolerance <= loc <= upper + tolerance
        and label.get_visible()
        and label.get_text()
    ]

    return [loc for loc, _ in in_view], [label for _, label in in_view]


def thin_keeping_ends(locs: list[float]) -> list[float]:
    """
    Take every other tick, keeping the first and the last

    The ends of an axis are often labelled on purpose (to show where it stops),
    so they are the last ticks to go.

    Parameters
    ----------
    locs
        Ticks to thin, in order

    Returns
    -------
    :
        Every other tick of `locs`, always including the first,
        and the last too if there were more than two
    """
    keep = list(locs[::2])
    if len(locs) > 2 and keep[-1] != locs[-1]:  # noqa: PLR2004
        # Swapped for the tick before it, so the gap stays at least two ticks wide
        keep[-1] = locs[-1]

    return keep


def _thin_crowded_tick_labels(
    fig: matplotlib.figure.Figure,
    axes: Mapping[str, matplotlib.axes.Axes],
    min_gap: float = 2.0,
    max_rounds: int = 5,
) -> None:
    """
    Take every other tick off any axis whose tick labels run into each other

    Matplotlib spaces ticks by the font size alone,
    not by how long the labels are,
    so on an axis only a couple of centimetres long
    labels like "1500" or "-2.5" end up touching.
    Which axes that happens to depends on how wide each panel ends up,
    which isn't known until the figure is laid out,
    so this is done as part of the layout rather than panel by panel.

    Parameters
    ----------
    fig
        Figure of interest

    axes
        The figure's axes

    min_gap
        Smallest gap to leave between neighbouring labels, in points

    max_rounds
        Most times to thin any one axis
    """
    renderer = fig.canvas.get_renderer()
    min_gap_pixels = min_gap * fig.dpi / 72.0
    to_check = [
        (axis, horizontal)
        for ax in axes.values()
        # 3D axes place their own ticks around a box, not along a line
        if ax.name != "3d" and ax.get_visible()
        for axis, horizontal in ((ax.xaxis, True), (ax.yaxis, False))
        if axis.get_visible()
    ]
    for _ in range(max_rounds):
        # Once per round rather than once per axis: drawing is the slow part
        fig.draw_without_rendering()
        crowded_axes = []
        for axis, horizontal in to_check:
            locs, labels = _get_tick_labels_in_view(axis)
            extents = sorted(
                (box.x0, box.x1) if horizontal else (box.y0, box.y1)
                for box in (label.get_window_extent(renderer) for label in labels)
            )
            if any(
                following[0] < preceding[1] + min_gap_pixels
                for preceding, following in itertools.pairwise(extents)
            ):
                crowded_axes.append((axis, locs))

        if not crowded_axes:
            break

        for axis, locs in crowded_axes:
            limits = axis.get_view_interval()
            axis.set_major_locator(
                matplotlib.ticker.FixedLocator(thin_keeping_ends(locs))
            )
            # Fixing the ticks can reset the limits they were chosen for
            axis.set_view_interval(*limits, ignore=True)

        # Only the axes which were crowded can still be
        to_check = [
            (axis, horizontal)
            for axis, horizontal in to_check
            if any(axis is crowded for crowded, _ in crowded_axes)
        ]


def lay_out_figure(
    fig: matplotlib.figure.Figure,
    axes: Mapping[str, matplotlib.axes.Axes],
    rows: tuple[Row, ...],
    settings: LayoutSettings = LayoutSettings(),
    n_iterations: int = 5,
) -> None:
    """
    Put every panel in its place

    This has to come after everything has been plotted,
    because it depends on how much room each panel's decorations take up.

    Parameters
    ----------
    fig
        Figure to lay out

    axes
        The figure's axes, as returned by [create_figure][]

        Only this figure's: every axes in here is drawn and measured
        as though it were part of `fig`.

    rows
        Rows of the figure, top to bottom

    settings
        Layout settings

    n_iterations
        Number of times to measure the decorations and re-place the panels

        The decorations barely change size as the panels move,
        so a handful of iterations is plenty.
    """
    pads = {panel.name: Pads() for row in rows for panel in row.panels}
    for _ in range(n_iterations):
        boxes = _place(fig, axes, rows, pads, settings)
        place_titles(axes, rows, pads, boxes, settings)
        _thin_crowded_tick_labels(fig, axes)
        pads = _measure(fig, axes, rows, boxes)

    boxes = _place(fig, axes, rows, pads, settings)
    place_titles(axes, rows, pads, boxes, settings)
    # The panels move a little in the last placement, so check again
    _thin_crowded_tick_labels(fig, axes)
    _warn_about_crowded_titles(fig, axes, rows, pads, settings)


def _warn_about_crowded_titles(  # noqa: PLR0913
    fig: matplotlib.figure.Figure,
    axes: Mapping[str, matplotlib.axes.Axes],
    rows: tuple[Row, ...],
    pads: Mapping[str, Pads],
    settings: LayoutSettings,
    max_cost: float = 0.2,
) -> None:
    """
    Warn about titles which run into each other or cost their panel width

    Neither stops the layout, but both are easy to miss.
    A title which reaches further past its panel
    than the next panel's decorations do
    pushes the next panel along, so the row's panels all get narrower
    (see the error [_place][] raises when that leaves them no width at all).
    The titles on the pieces of a broken panel
    are not kept apart by the layout at all.

    Parameters
    ----------
    fig
        Figure of interest

    axes
        The figure's axes

    rows
        Rows of the figure, top to bottom

    pads
        Room each panel's decorations take up, as measured by [_measure][]

    settings
        Layout settings

    max_cost
        Most width a title may cost its row before we warn, in inches

        A little is a fair price for a title which reads well,
        so this only catches the titles which cost a noticeable share
        of a panel's width.
    """
    fig.draw_without_rendering()
    renderer = fig.canvas.get_renderer()
    to_inches = fig.dpi_scale_trans.inverted()

    titles = {
        (name, title.get_text()): title.get_window_extent(renderer).transformed(
            to_inches
        )
        for name, ax in axes.items()
        for title in get_titles(ax)
        if title.get_visible()
    }
    for (key, box), (other_key, other_box) in itertools.combinations(titles.items(), 2):
        if box.overlaps(other_box):
            logger.warning(
                f"Titles overlap: {key} and {other_key}. Shorten or wrap one."
            )

    for row in rows:
        for i, panel in enumerate(row.panels):
            if i < len(row.panels) - 1:
                without_title = (
                    pads[panel.name].right
                    + settings.wspace
                    + pads[row.panels[i + 1].name].left
                )
                cost = pads[panel.name].title_right + settings.title_gap - without_title
            else:
                cost = pads[panel.name].title_right - pads[panel.name].right

            if cost > max_cost:
                logger.warning(
                    f"The title of {panel.name!r} costs its row {cost:.2f} inches "
                    "of width. Shorten or wrap it."
                )
