"""
Layout of the historical GHG manuscript's methods figures

The methods figures have panels whose shape is not up for negotiation
(the maps have a fixed aspect ratio, the flying carpet is square)
sat alongside panels which can be any shape at all.
Matplotlib's constrained layout engine ties every column of the figure
to every other one, which the fixed panels then over-constrain.

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

import itertools
import string
from collections.abc import Mapping
from dataclasses import dataclass

import matplotlib.axes
import matplotlib.figure
import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox


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

    If `True`, the panel is made of two axes,
    `f"{name}-l"` and `f"{name}-r"`, side by side with a small gap between.
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

    left: float = 0.7
    right: float = 0.2
    top: float = 0.35
    bottom: float = 0.5


@dataclass(frozen=True)
class LayoutSettings:
    """
    Spacing settings for the layout, in inches
    """

    width: float = 16.0
    """Width of the figure"""

    margin: float = 0.05
    """Space around the outside of the figure's decorations"""

    wspace: float = 0.2
    """Space between the decorations of neighbouring panels in a row"""

    hspace: float = 0.15
    """Space between the decorations of neighbouring rows"""

    broken_gap: float = 0.08
    """Gap between the two halves of a broken axis"""

    colour_bar_gap: float = 0.1
    """Gap between a panel and its colour bar"""

    colour_bar_width: float = 0.15
    """Width of a colour bar"""


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
        return (f"{panel.name}-l", f"{panel.name}-r")

    return (panel.name,)


def get_colour_bar_axes_name(panel: Panel) -> str:
    """
    Get the name of a panel's colour bar axes
    """
    return f"{panel.name}-colour-bar"


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
    """
    panels = [panel for row in rows for panel in row.panels]
    if set(titles) != {panel.name for panel in panels}:
        msg = f"Titles don't match the panels: {set(titles)=}"
        raise AssertionError(msg)

    for label, panel in zip(string.ascii_lowercase, panels):
        # Titles go on the leftmost axes of the panel
        ax = axes[get_panel_axes_names(panel)[0]]
        ax.set_title(
            f"$\\bf{{({label})}}$ {titles[panel.name]}",
            loc="left",
            fontsize="medium",
        )


def _place(  # noqa: PLR0912
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
    # Every row's data starts at the same place on the left,
    # so the vertical axes line up.
    # On the right, each row runs all the way to the edge,
    # otherwise a colour bar in one row would leave a gap at the end of all the others.
    left_edge = settings.margin + max(pads[row.panels[0].name].left for row in rows)

    # Positions from the top first, because we don't know the figure's height yet
    boxes_from_top = {}
    y = settings.margin
    for row in rows:
        panels = row.panels
        right_edge = settings.width - settings.margin - pads[panels[-1].name].right
        between = sum(
            pads[left.name].right + settings.wspace + pads[right.name].left
            for left, right in itertools.pairwise(panels)
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
            msg = f"No width left for the flexible panels in {row}"
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
                x += (
                    width
                    + pads[panel.name].right
                    + settings.wspace
                    + pads[panels[i + 1].name].left
                )

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
                half = (width - settings.broken_gap) / 2.0
                axes[names[0]].set_position(to_figure_coords(x0, y0, half, height))
                axes[names[1]].set_position(
                    to_figure_coords(x0 + half + settings.broken_gap, y0, half, height)
                )
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
            names = list(get_panel_axes_names(panel))
            if panel.colour_bar:
                names.append(get_colour_bar_axes_name(panel))

            tight = Bbox.union(
                [axes[name].get_tightbbox(renderer) for name in names]
            ).transformed(fig.dpi_scale_trans.inverted())

            x0, y0, width, height = boxes[panel.name]
            pads[panel.name] = Pads(
                left=max(x0 - tight.x0, 0.0),
                right=max(tight.x1 - (x0 + width), 0.0),
                top=max(tight.y1 - (y0 + height), 0.0),
                bottom=max(y0 - tight.y0, 0.0),
            )

    return pads


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
        pads = _measure(fig, axes, rows, boxes)

    _place(fig, axes, rows, pads, settings)
