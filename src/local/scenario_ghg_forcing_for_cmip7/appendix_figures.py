"""
Generation of the scenario manuscript's per-gas appendix figures

Each group of gases (processed like N2O, CFC-12 or HFC-125) has the same figure,
so rather than write them out by hand, we write them here
and the build inlines them in the appendices in place of their tags.
Their captions point to the figure of the same kind in the methods,
rather than repeat its caption.
The same is done for the comparison of MAGICC's concentrations with ours
for the gases which use a single concentration projection.

Gas names are written the way the manuscript's source text writes them
(e.g. `C_2F_6`, `CFC-12`), so the manuscript's replacements
format them exactly as they are formatted in the text around them.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable
from pathlib import Path

from loguru import logger

from local.historical_ghg_forcing_for_cmip7.appendix_figures import (
    AppendixFigure,
    get_figure_tags_in_text,
)
from local.scenario_ghg_forcing_for_cmip7.rerun_projections import (
    ONE_EOF_GASES,
    SINGLE_PROJECTION_GASES,
)

MANUSCRIPT_GAS_NAMES = {
    "co2": "CO_2",
    "ch4": "CH_4",
    "n2o": "N_2O",
    "c2f6": "C_2F_6",
    "c3f8": "C_3F_8",
    "c4f10": "C_4F_10",
    "c5f12": "C_5F_12",
    "c6f14": "C_6F_14",
    "c7f16": "C_7F_16",
    "c8f18": "C_8F_18",
    "cc4f8": "cC_4F_8",
    "ccl4": "CCl_4",
    "cf4": "CF_4",
    "cfc11": "CFC-11",
    "cfc113": "CFC-113",
    "cfc114": "CFC-114",
    "cfc115": "CFC-115",
    "cfc12": "CFC-12",
    "ch2cl2": "CH_2Cl_2",
    "ch3br": "CH_3Br",
    "ch3ccl3": "CH_3CCl_3",
    "ch3cl": "CH_3Cl",
    "chcl3": "CHCl_3",
    "halon1211": "Halon-1211",
    "halon1301": "Halon-1301",
    "halon2402": "Halon-2402",
    "hcfc141b": "HCFC-141b",
    "hcfc142b": "HCFC-142b",
    "hcfc22": "HCFC-22",
    "hfc125": "HFC-125",
    "hfc134a": "HFC-134a",
    "hfc143a": "HFC-143a",
    "hfc152a": "HFC-152a",
    "hfc227ea": "HFC-227ea",
    "hfc23": "HFC-23",
    "hfc236fa": "HFC-236fa",
    "hfc245fa": "HFC-245fa",
    "hfc32": "HFC-32",
    "hfc365mfc": "HFC-365mfc",
    "hfc4310mee": "HFC-43-10mee",
    "nf3": "NF_3",
    "sf6": "SF_6",
    "so2f2": "SO_2F_2",
}
"""How the manuscript's source text writes each gas' name

These are keys of the manuscript's `replacements.yaml`.
"""

N2O_LIKE_GASES = ONE_EOF_GASES
"""Gases processed like N2O, other than N2O, in the order the methods list them"""

CFC12_LIKE_GASES = tuple(g for g in SINGLE_PROJECTION_GASES if g != "cfc12")
"""Gases processed like CFC-12, other than CFC-12, in the order the methods list them"""

HFC125_LIKE_GASES = (
    "c2f6",
    "c3f8",
    "c4f10",
    "c5f12",
    "c6f14",
    "c7f16",
    "c8f18",
    "cc4f8",
    "cf4",
    "ch2cl2",
    "chcl3",
    "hfc143a",
    "hfc227ea",
    "hfc236fa",
    "hfc4310mee",
    "nf3",
    "sf6",
    "so2f2",
)
"""Gases processed like HFC-125, other than HFC-125, in the methods' order"""

MAIN_TEXT_GASES = ("co2", "ch4", "n2o", "cfc12", "hfc125")
"""Gases whose methods figure is in the main text"""

MAGICC_MISMATCH_GASES = SINGLE_PROJECTION_GASES
"""Gases whose concentrations we compare with MAGICC's"""


def get_methods_figure_tag(gas: str) -> str:
    """Get the tag of a gas' methods figure"""
    return f"@{gas}-methods-figure@"


def get_methods_figure_label(gas: str) -> str:
    """Get the label of a gas' methods figure"""
    return f"fig:methods-{gas}"


def get_magicc_mismatch_figure_tag(gas: str) -> str:
    """Get the tag of a gas' comparison of MAGICC's concentrations with ours"""
    return f"@{gas}-magicc-mismatch-figure@"


REGRESSION_GASES = ("co2", "ch4", "n2o")
"""Gases whose regressions we re-derive (so have a regressions figure)"""


def get_regressions_figure_tag(gas: str) -> str:
    """Get the tag of the figure of the regressions we re-derived for a gas"""
    return f"@{gas}-regressions-figure@"


WESTERN_ET_AL_2024_FIGURE_TAG = "@western-et-al-2024-extensions-figure@"
"""Tag of the figure of our extensions of Western et al. (2024)"""


def get_like_caption(reference_gas: str, gas: str, extra: str = "") -> str:
    """Get a caption which says a figure is like a methods figure, for another gas"""
    res = (
        rf"Like Fig.~\ref{{{get_methods_figure_label(reference_gas)}}}, "
        f"except for {MANUSCRIPT_GAS_NAMES[gas]}."
    )
    if extra:
        res = f"{res} {extra}"

    return res


def get_n2o_like_figures(gas: str) -> tuple[AppendixFigure, ...]:
    """Get the appendix figures of a gas processed like N2O"""
    return (
        AppendixFigure(
            tag=get_methods_figure_tag(gas),
            caption=get_like_caption(
                "n2o",
                gas,
                "Unlike N_2O, only one latitudinal gradient EOF is retained, "
                "so there is only one PC in d).",
            ),
            label=get_methods_figure_label(gas),
        ),
    )


def get_cfc12_like_figures(gas: str) -> tuple[AppendixFigure, ...]:
    """Get the appendix figures of a gas processed like CFC-12"""
    extra = ""
    if gas in ("hcfc22", "hcfc141b", "hcfc142b"):
        extra = (
            r"The projection comes from \citet{western_hcfcs_2024}, "
            "which ends in 2100. "
            "The shading marks our extension of it "
            r"(see also Fig.~\ref{fig:methods-western-et-al-2024})."
        )

    return (
        AppendixFigure(
            tag=get_methods_figure_tag(gas),
            caption=get_like_caption("cfc12", gas, extra),
            label=get_methods_figure_label(gas),
        ),
    )


def get_hfc125_like_figures(gas: str) -> tuple[AppendixFigure, ...]:
    """Get the appendix figures of a gas processed like HFC-125"""
    extra = ""
    if gas == "c8f18":
        extra = (
            "In the historical dataset, the C_8F_18 latitudinal gradient "
            "has an assumed shape, scaled to reproduce the CMIP6 hemispheric-means "
            "[TODO cite historical paper]. "
            "There is no regression between it and emissions, "
            "so in the projections it is held constant "
            "at its value in December 2022, "
            "i.e. the EOF is the December 2022 latitudinal gradient "
            "and the PC is one. "
            "The historical PC shown is the historical latitudinal gradient "
            "expressed relative to this EOF."
        )

    return (
        AppendixFigure(
            tag=get_methods_figure_tag(gas),
            caption=get_like_caption("hfc125", gas, extra),
            label=get_methods_figure_label(gas),
        ),
    )


def get_magicc_mismatch_figures(gas: str) -> tuple[AppendixFigure, ...]:
    """Get the comparison of MAGICC's concentrations of a gas with ours"""
    name = MANUSCRIPT_GAS_NAMES[gas]
    return (
        AppendixFigure(
            tag=get_magicc_mismatch_figure_tag(gas),
            caption=(
                f"{name} emissions and concentrations: "
                "the CMIP7 forcings compared with those used in, "
                "and calculated by, MAGICC. "
                f"a) The emissions implied by the CMIP7 {name} forcings "
                "(inverse emissions) "
                "and the official emissions "
                "(including the official extensions beyond 2100), "
                "with which MAGICC is run, in each dataset. "
                f"b) The CMIP7 {name} forcings (the same in all datasets), "
                "with the historical dataset, "
                "and the median of MAGICC's concentrations in each dataset."
            ),
            label=f"fig:magicc-mismatch-{gas}",
        ),
    )


APPENDIX_FIGURE_GROUPS: dict[
    str, tuple[Callable[[str], tuple[AppendixFigure, ...]], tuple[str, ...]]
] = {
    "n2o-like-methods": (get_n2o_like_figures, N2O_LIKE_GASES),
    "cfc12-like-methods": (get_cfc12_like_figures, CFC12_LIKE_GASES),
    "hfc125-like-methods": (get_hfc125_like_figures, HFC125_LIKE_GASES),
    "magicc-mismatch": (get_magicc_mismatch_figures, MAGICC_MISMATCH_GASES),
}
"""Function which gets a gas' figures, and the gases, for each group of figures"""


def get_appendix_figures_tag(group: str) -> str:
    """Get the tag the build replaces with a group's appendix figures"""
    return f"@{group}-appendix-figures@"


def get_generated_figure_tags(group: str) -> tuple[str, ...]:
    """Get the tags of all the figures a group of appendix figures can include"""
    get_figures, gases = APPENDIX_FIGURE_GROUPS[group]
    return tuple(figure.tag for gas in gases for figure in get_figures(gas))


def generate_appendix_figures(
    out_file: Path,
    group: str,
    figure_tags_in_text: set[str],
) -> Path:
    """
    Generate the latex for a group of appendix figures

    Parameters
    ----------
    out_file
        File in which to write the latex

    group
        Group of figures to generate, a key of [APPENDIX_FIGURE_GROUPS][]

    figure_tags_in_text
        Tags of the figures which the text already includes,
        which are left out so they aren't in the manuscript twice

    Returns
    -------
    :
        `out_file`
    """
    get_figures, gases = APPENDIX_FIGURE_GROUPS[group]

    gases_tex = []
    for gas in gases:
        figures = [f for f in get_figures(gas) if f.tag not in figure_tags_in_text]
        figures_tex = [
            # The first figure stays after the heading of its appendix
            figure.to_tex(placement="hbt!" if not gases_tex and j == 0 else "tp")
            for j, figure in enumerate(figures)
        ]
        if figures_tex:
            gases_tex.append("\n\n".join(figures_tex))

    out_file.parent.mkdir(parents=True, exist_ok=True)
    # A new page every few figures,
    # otherwise latex runs out of room for the figures it hasn't placed yet
    # ("Too many unprocessed floats").
    out_file.write_text("\n\n\\clearpage\n\n".join(gases_tex) + "\n")
    logger.info(f"Wrote {out_file}")

    return out_file


def get_tags_in_texts(texts: Iterable[Path]) -> set[str]:
    """Get the tags of the figures included in some latex files"""
    return get_figure_tags_in_text(t.read_text() for t in texts)
