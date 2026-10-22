"""
Paths used by the scenario GHG manuscript's build

The build is driven by pytask (see `tasks/scenario_ghg_forcing_for_cmip7`).
These live here, rather than next to the tasks,
because pytask doesn't put the task modules' directories on the import path,
so the task modules can't import a module which sits next to them.
"""

from __future__ import annotations

from local.paths import FIGURES_DIR, REPO_ROOT
from local.scenario_ghg_forcing_for_cmip7.bundle import (
    BUNDLE_DIR,
    INVERSE_EMISSIONS_FILE,
    INVERSE_EMISSIONS_INTERMEDIATE_FILE,
    INVERSE_EMISSIONS_NOTEBOOK,
    MODIFIED_NOTEBOOKS_DIR,
    MONTHLY_MEAN_FILE,
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

# Manuscript sources
MANUSCRIPT_DIR = REPO_ROOT / "manuscripts" / MANUSCRIPT_NAME
"""Directory which holds the manuscript's source"""

ABSTRACT = MANUSCRIPT_DIR / "abstract.tex"
INTRODUCTION = MANUSCRIPT_DIR / "introduction.tex"
SECTIONS = (
    MANUSCRIPT_DIR / "output-requirements.tex",
    MANUSCRIPT_DIR / "methods.tex",
)
"""Body sections, in the order they appear"""

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

INVERSE_EMISSIONS_FIGURE = FIGURES_OUT_DIR / "ccl4_inverse_emissions.pdf"
SUMMARY_TABLE = TABLES_OUT_DIR / "ccl4_summary.tex"
AUTHOR_CONTRIBUTIONS = BUILD_DIR / "author-contributions.tex"
"""The author contributions statement, generated from the metadata file"""

FIGURE_FILES = {
    "@ccl4-inverse-emissions-figure@": INVERSE_EMISSIONS_FIGURE,
}
"""Tag in the latex which each figure replaces"""

FILES_TO_INLINE = {
    "@ccl4-summary-table@": SUMMARY_TABLE,
    "@author-contributions@": AUTHOR_CONTRIBUTIONS,
}
"""Tag in the latex which each file to inline replaces"""

# Outputs
LATEX_DIR = BUILD_DIR / "latex"
LATEX_MAIN = LATEX_DIR / "main.tex"
OUTPUT_PDF = REPO_ROOT / "compiled-manuscripts" / f"{MANUSCRIPT_NAME}.pdf"
