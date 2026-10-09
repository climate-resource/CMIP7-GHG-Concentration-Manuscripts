"""
Generation of the appendix figures for the gases which share a figure layout

The gases processed like CFC-12 (and those processed like C4F10)
each have the same figures.
Rather than write these out by hand for thirty-odd gases,
we write them here and the build inlines them in the appendices
in place of their tags.
Their captions point to the figure of the same kind
which is in the text (e.g. the CFC-12 methods figure),
rather than repeat its caption.

The gas names are written the way the manuscript's source text writes them
(e.g. `C_2F_6`, `CFC12`), so the manuscript's replacements
format them exactly as they are formatted in the text around them.
"""

from __future__ import annotations

import re
from collections.abc import Callable, Iterable
from dataclasses import dataclass
from pathlib import Path

from loguru import logger

from local.historical_ghg_forcing_for_cmip7.cfc12_like_tables import (
    MANUSCRIPT_GAS_NAMES as CFC12_LIKE_MANUSCRIPT_GAS_NAMES,
)

C4F10_LIKE_MANUSCRIPT_GAS_NAMES = {
    "c4f10": "C_4F_10",
    "c5f12": "C_5F_12",
    "c6f14": "C_6F_14",
    "c7f16": "C_7F_16",
    "cc4f8": "cC_4F_8",
}
"""How the manuscript's source text writes the name of each gas processed like C4F10

These are keys of the manuscript's `replacements.yaml`,
which is where they are turned into tex.
"""

MANUSCRIPT_GAS_NAMES = {
    **CFC12_LIKE_MANUSCRIPT_GAS_NAMES,
    **C4F10_LIKE_MANUSCRIPT_GAS_NAMES,
}
"""How the manuscript's source text writes each gas' name"""

FIGURE_TAG = re.compile(r"\\includegraphics(?:\[[^\]]*\])?\{(@[^}@]*@)\}")
"""A figure included in the latex, which is still a tag"""

COMMENT = re.compile(r"(?<!\\)%.*")
"""A latex comment (`%` up to the end of the line, unless the `%` is escaped)"""


@dataclass(frozen=True)
class AppendixFigure:
    """
    A figure to put in an appendix
    """

    tag: str
    """
    Tag which the build replaces with the figure's file
    """

    caption: str
    """
    Caption of the figure, as latex
    """

    label: str
    """
    Label of the figure
    """

    def to_tex(self, placement: str) -> str:
        """
        Get the figure as latex

        Parameters
        ----------
        placement
            Float placement specifier e.g. `tp`

        Returns
        -------
        :
            The figure's latex
        """
        return "\n".join(
            [
                rf"\begin{{figure}}[{placement}]",
                rf"\includegraphics[width=12cm]{{{self.tag}}}",
                r"\caption{",
                f"    {self.caption}",
                "}",
                rf"\label{{{self.label}}}",
                r"\end{figure}",
            ]
        )


def get_like_caption(reference_label: str, gas: str) -> str:
    """
    Get a caption which says a figure is like another figure, but for another gas

    Parameters
    ----------
    reference_label
        Label of the figure this figure is like

    gas
        Gas this figure is for

    Returns
    -------
    :
        The caption
    """
    return (
        rf"Like Figure~\ref{{{reference_label}}}, "
        f"except for {MANUSCRIPT_GAS_NAMES[gas]}."
    )


def get_cfc12_like_methods_figures(gas: str) -> tuple[AppendixFigure, ...]:
    """Get the methods figures of a gas processed like CFC-12"""
    return (
        AppendixFigure(
            tag=f"@{gas}-methods-figure@",
            caption=get_like_caption("fig:methods-cfc12", gas),
            label=f"fig:methods-{gas}",
        ),
        AppendixFigure(
            tag=f"@{gas}-methods-appendix-figure@",
            caption=get_like_caption("fig:methods-appendix-cfc12", gas),
            label=f"fig:methods-appendix-{gas}",
        ),
    )


def get_c4f10_like_methods_figures(gas: str) -> tuple[AppendixFigure, ...]:
    """Get the methods figures of a gas processed like C4F10"""
    return (
        AppendixFigure(
            tag=f"@{gas}-methods-figure@",
            caption=get_like_caption("fig:methods-c4f10", gas),
            label=f"fig:methods-{gas}",
        ),
    )


def get_cfc12_like_results_figures(gas: str) -> tuple[AppendixFigure, ...]:
    """Get the results figures of a gas processed like CFC-12"""
    return (
        AppendixFigure(
            tag=f"@{gas}-results-figure@",
            caption=get_like_caption("fig:results-cfc12", gas),
            label=f"fig:results-{gas}",
        ),
    )


def get_c4f10_like_results_figures(gas: str) -> tuple[AppendixFigure, ...]:
    """Get the results figures of a gas processed like C4F10"""
    return (
        AppendixFigure(
            tag=f"@{gas}-results-figure@",
            caption=get_like_caption("fig:results-c4f10", gas),
            label=f"fig:results-{gas}",
        ),
    )


APPENDIX_FIGURE_GETTERS: dict[str, Callable[[str], tuple[AppendixFigure, ...]]] = {
    "cfc12-like-methods": get_cfc12_like_methods_figures,
    "c4f10-like-methods": get_c4f10_like_methods_figures,
    "cfc12-like-results": get_cfc12_like_results_figures,
    "c4f10-like-results": get_c4f10_like_results_figures,
}
"""Function which gets a gas' figures, for each group of appendix figures"""


def get_figure_tags_in_text(texts: Iterable[str]) -> set[str]:
    """
    Get the tags of the figures which are included in some latex

    Parameters
    ----------
    texts
        Latex to search

    Returns
    -------
    :
        Tags of the figures included in `texts` (ignoring comments)
    """
    return {
        match.group(1)
        for text in texts
        for line in text.splitlines()
        for match in FIGURE_TAG.finditer(COMMENT.sub("", line))
    }


def generate_appendix_figures(
    out_file: Path,
    group: str,
    gases: Iterable[str],
    figure_tags_in_text: set[str],
) -> Path:
    """
    Generate the latex for a group of appendix figures

    Unlike the tables, this is always re-generated:
    it doesn't depend on any data, so it is quick
    and can't go stale when the text changes.

    Parameters
    ----------
    out_file
        File in which to write the latex

    group
        Group of figures to generate, a key of [APPENDIX_FIGURE_GETTERS][]

    gases
        Gases to generate the figures for, in the order they should appear

    figure_tags_in_text
        Tags of the figures which the text already includes.

        These are left out, so they aren't in the manuscript twice.
        For example, CFC-12's methods figures are in the methods,
        so they are not repeated with the other gases processed like CFC-12.

    Returns
    -------
    :
        `out_file`
    """
    get_figures = APPENDIX_FIGURE_GETTERS[group]

    gases_tex = []
    for gas in gases:
        figures = [
            figure
            for figure in get_figures(gas)
            if figure.tag not in figure_tags_in_text
        ]
        figures_tex = [
            # Get latex to keep the first figure after the heading of the appendix,
            # or at the bottom of this page or at the top of the next
            # (with `tp` it goes at the top of the page, above the heading).
            # The ! because it is bigger than latex normally allows.
            figure.to_tex(placement="hbt!" if not gases_tex and j == 0 else "tp")
            for j, figure in enumerate(figures)
        ]
        if figures_tex:
            gases_tex.append("\n\n".join(figures_tex))

    out_file.parent.mkdir(parents=True, exist_ok=True)
    out_file.write_text(
        # A new page for each gas.
        # Otherwise latex runs out of room for the figures it hasn't placed yet
        # ("Too many unprocessed floats"), and it keeps each gas' figures together.
        "\n\n\\clearpage\n\n".join(gases_tex) + "\n"
    )
    logger.info(f"Wrote {out_file}")

    return out_file
