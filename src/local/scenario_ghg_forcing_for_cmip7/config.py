"""
Paths used by the scenario GHG manuscript's build

The build is driven by pytask (see `tasks/scenario_ghg_forcing_for_cmip7`).
These live here, rather than next to the tasks,
because pytask doesn't put the task modules' directories on the import path,
so the task modules can't import a module which sits next to them.
"""

from __future__ import annotations

from pathlib import Path

from local.cmip_ghg_generation import DEFAULT_BUNDLE_DIR as HISTORICAL_BUNDLE_DIR
from local.paths import FIGURES_DIR, REPO_ROOT
from local.scenario_ghg_forcing_for_cmip7 import rerun_vl_cf
from local.scenario_ghg_forcing_for_cmip7.appendix_figures import (
    APPENDIX_FIGURE_GROUPS,
    CFC12_LIKE_GASES,
    HFC125_LIKE_GASES,
    MAGICC_MISMATCH_GASES,
    MAIN_TEXT_GASES,
    N2O_LIKE_GASES,
    REGRESSION_GASES,
    WESTERN_ET_AL_2024_FIGURE_TAG,
    get_appendix_figures_tag,
    get_magicc_mismatch_figure_tag,
    get_methods_figure_tag,
    get_regressions_figure_tag,
)
from local.scenario_ghg_forcing_for_cmip7.bundle import (
    BUNDLE_DIR,
    BUNDLE_TARBALLS,
    INVERSE_EMISSIONS_FILE,
    INVERSE_EMISSIONS_INTERMEDIATE_FILE,
    INVERSE_EMISSIONS_NOTEBOOK,
    MODIFIED_NOTEBOOKS_DIR,
    MONTHLY_MEAN_FILE,
)
from local.scenario_ghg_forcing_for_cmip7.rerun_projections import (
    CO2_SEASONALITY_PC_FILE,
    CO2_SEASONALITY_REGRESSION_KEY,
    FOSSIL_BIO_SPLIT_FILE,
    ONE_EOF_GASES,
    SINGLE_PROJECTION_GASES,
    TWO_EOF_GASES,
    get_eof_scaling_emissions_file,
    get_lat_gradient_info_file,
    get_pcs_file,
    get_regression_file,
    get_regression_info_file,
)

MANUSCRIPT_NAME = "scenario-ghg-forcing-for-cmip7"
"""Name of the manuscript, used for its directories and the output PDF"""

# Bundle
BUNDLE_PYPROJECT = BUNDLE_DIR / "pyproject.toml"
"""The bundle's `pyproject.toml`"""

BUNDLE_PIXI_LOCK = BUNDLE_DIR / "pixi.lock"
"""The bundle's lock file, which defines the environment the notebooks run in"""

BUNDLE_MONTHLY_MEAN = BUNDLE_DIR / MONTHLY_MEAN_FILE
"""Input to the inverse emissions notebook"""

BUNDLE_INVERSE_EMISSIONS = BUNDLE_DIR / INVERSE_EMISSIONS_FILE
"""The re-run inverse emissions notebook's usual output"""

BUNDLE_INVERSE_EMISSIONS_INTERMEDIATE = BUNDLE_DIR / INVERSE_EMISSIONS_INTERMEDIATE_FILE
"""Intermediates saved by the re-run inverse emissions notebook"""

EXECUTED_INVERSE_EMISSIONS_NOTEBOOK = INVERSE_EMISSIONS_NOTEBOOK
"""The inverse emissions notebook, as executed in the original run"""

MODIFIED_INVERSE_EMISSIONS_NOTEBOOK = (
    MODIFIED_NOTEBOOKS_DIR / f"{INVERSE_EMISSIONS_NOTEBOOK.stem}.py"
)
"""The modified inverse emissions notebook, in `py:percent` format"""

BUNDLE_EXTRACTED_MARKERS = [
    BUNDLE_DIR / "zenodo-tarballs" / f"{tarball}.extracted"
    for tarball in BUNDLE_TARBALLS
]
"""Files which mark that each of the bundle's tarballs has been extracted"""

HISTORICAL_BUNDLE_TARBALLS = ("data--interim.tar.gz", "data--processed.tar.gz")
"""Tarballs we need from the historical run's bundle

Its interim data (the historical PCs etc.)
and its copy of the published historical files.
"""

HISTORICAL_BUNDLE_EXTRACTED_MARKERS = [
    HISTORICAL_BUNDLE_DIR / "zenodo-tarballs" / f"{tarball}.extracted"
    for tarball in HISTORICAL_BUNDLE_TARBALLS
]
"""Files which mark that each of the historical bundle's tarballs has been extracted"""

MAGICC_SUMMARY = BUNDLE_DIR / "manuscript-outputs" / "magicc-summary.feather"
"""Summary of the original run's MAGICC output"""

RERUN_EMISSIONS_FILES = [
    BUNDLE_DIR / FOSSIL_BIO_SPLIT_FILE,
    *(
        BUNDLE_DIR / get_eof_scaling_emissions_file(g)
        for g in (*TWO_EOF_GASES, *ONE_EOF_GASES)
    ),
]
"""Emissions extracted for the re-runs"""

RERUN_CO2_SEASONALITY_FILES = [
    BUNDLE_DIR / CO2_SEASONALITY_PC_FILE,
    BUNDLE_DIR / get_regression_file(CO2_SEASONALITY_REGRESSION_KEY),
    BUNDLE_DIR / get_regression_info_file(CO2_SEASONALITY_REGRESSION_KEY),
]
"""What the re-run of the CO2 seasonality projection writes that we use"""

RERUN_CO2_SEASONALITY_VL_CF_FILES = [
    BUNDLE_DIR / f for f in rerun_vl_cf.RERUN_CO2_SEASONALITY_VL_CF_FILES
]
"""What the re-run of the CO2 seasonality projection for VL-CF writes that we use"""


def get_rerun_lat_gradient_files(ghg: str) -> list[Path]:
    """Get the outputs we use of the re-run of a gas' latitudinal gradient projection"""
    res = [BUNDLE_DIR / get_lat_gradient_info_file(ghg), BUNDLE_DIR / get_pcs_file(ghg)]
    if ghg in TWO_EOF_GASES:
        res.extend(
            [
                BUNDLE_DIR / get_regression_file(ghg),
                BUNDLE_DIR / get_regression_info_file(ghg),
            ]
        )

    return res


RERUN_LAT_GRADIENT_GASES = (*TWO_EOF_GASES, *ONE_EOF_GASES, *SINGLE_PROJECTION_GASES)
"""Gases whose latitudinal gradient projections we re-run"""

# Manuscript sources
MANUSCRIPT_DIR = REPO_ROOT / "manuscripts" / MANUSCRIPT_NAME
"""Directory which holds the manuscript's source"""

ABSTRACT = MANUSCRIPT_DIR / "abstract.tex"
INTRODUCTION = MANUSCRIPT_DIR / "introduction.tex"
SECTIONS = (
    MANUSCRIPT_DIR / "output-requirements.tex",
    MANUSCRIPT_DIR / "methods.tex",
    MANUSCRIPT_DIR / "discussion.tex",
)
"""Body sections, in the order they appear"""

APPENDICES = (MANUSCRIPT_DIR / "appendices.tex",)
"""Appendices, in the order they appear"""

CONCLUSION = MANUSCRIPT_DIR / "conclusion.tex"
CODE_AND_DATA_AVAILABILITY = MANUSCRIPT_DIR / "code-and-data-availability.tex"
AUTHOR_CONTRIBUTION = MANUSCRIPT_DIR / "author-contribution.tex"
COMPETING_INTERESTS = MANUSCRIPT_DIR / "competing-interests.tex"
ACKNOWLEDGEMENTS = MANUSCRIPT_DIR / "acknowledgements.tex"
METADATA = MANUSCRIPT_DIR / "metadata.toml"
REPLACEMENTS = MANUSCRIPT_DIR / "replacements.yaml"

REFERENCES_BIB = REPO_ROOT / "references" / "references.bib"

COPERNICUS_TEMPLATE_DIR = REPO_ROOT / "copernicus-latex-package"
CLEAN_COPERNICUS_TEMPLATE_FILENAME = "template_clean.tex"

# Generated tex inputs
FIGURES_OUT_DIR = FIGURES_DIR / MANUSCRIPT_NAME
TABLES_OUT_DIR = REPO_ROOT / "tables" / MANUSCRIPT_NAME
BUILD_DIR = REPO_ROOT / "build" / MANUSCRIPT_NAME

AUTHOR_CONTRIBUTIONS = BUILD_DIR / "author-contributions.tex"
"""The author contributions statement, generated from the metadata file"""

METHODS_FIGURES_DIR = FIGURES_OUT_DIR / "methods"
"""Directory which holds the methods figures"""


def get_methods_figure_file(gas: str) -> Path:
    """Get the file which holds a gas' methods figure"""
    return METHODS_FIGURES_DIR / f"{gas}_methods.pdf"


def get_magicc_mismatch_figure_file(gas: str) -> Path:
    """Get the file of a gas' comparison of MAGICC's concentrations with ours"""
    return METHODS_FIGURES_DIR / f"{gas}_magicc-mismatch.pdf"


def get_regressions_figure_file(gas: str) -> Path:
    """Get the file of the figure of the regressions we re-derived for a gas"""
    return METHODS_FIGURES_DIR / f"{gas}_regressions.pdf"


WESTERN_ET_AL_2024_FIGURE = METHODS_FIGURES_DIR / "western-et-al-2024-extensions.pdf"
"""The figure of our extensions of the Western et al. (2024) projections"""

METHODS_FIGURE_GASES = (
    *MAIN_TEXT_GASES,
    *N2O_LIKE_GASES,
    *CFC12_LIKE_GASES,
    *HFC125_LIKE_GASES,
)
"""Gases which have a methods figure"""

FIGURE_FILES: dict[str, Path] = {
    **{
        get_methods_figure_tag(g): get_methods_figure_file(g)
        for g in METHODS_FIGURE_GASES
    },
    **{
        get_magicc_mismatch_figure_tag(g): get_magicc_mismatch_figure_file(g)
        for g in MAGICC_MISMATCH_GASES
    },
    **{
        get_regressions_figure_tag(g): get_regressions_figure_file(g)
        for g in REGRESSION_GASES
    },
    WESTERN_ET_AL_2024_FIGURE_TAG: WESTERN_ET_AL_2024_FIGURE,
}
"""Tag in the latex which each figure replaces

The main text's figures are only included in the build if their tags are in the text.
"""

APPENDIX_FIGURES_FILES = {
    group: BUILD_DIR / f"{group}-appendix-figures.tex"
    for group in APPENDIX_FIGURE_GROUPS
}
"""File which holds the latex of each group of appendix figures"""

FILES_TO_INLINE = {
    "@author-contributions@": AUTHOR_CONTRIBUTIONS,
    **{
        get_appendix_figures_tag(group): f
        for group, f in APPENDIX_FIGURES_FILES.items()
    },
}
"""Tag in the latex which each file to inline replaces"""

# Outputs
LATEX_DIR = BUILD_DIR / "latex"
LATEX_MAIN = LATEX_DIR / "main.tex"
OUTPUT_PDF = REPO_ROOT / "compiled-manuscripts" / f"{MANUSCRIPT_NAME}.pdf"
