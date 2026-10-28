"""
Tasks which assemble the manuscript and compile it to PDF
"""

import shutil
import sys
from pathlib import Path
from typing import Annotated

from pytask import Product

import local.gmd_latex
import local.scenario_ghg_forcing_for_cmip7.value_checks
from local.gmd_latex import (
    COPERNICUS_TEMPLATE_FILES,
    assemble_main_tex,
    compile_latex,
)
from local.scenario_ghg_forcing_for_cmip7.appendix_figures import (
    APPENDIX_FIGURE_GROUPS,
    get_generated_figure_tags,
    get_tags_in_texts,
)
from local.scenario_ghg_forcing_for_cmip7.bundle import BUNDLE_DIR
from local.scenario_ghg_forcing_for_cmip7.config import (
    ABSTRACT,
    ACKNOWLEDGEMENTS,
    APPENDICES,
    AUTHOR_CONTRIBUTION,
    BUILD_DIR,
    CLEAN_COPERNICUS_TEMPLATE_FILENAME,
    CODE_AND_DATA_AVAILABILITY,
    COMPETING_INTERESTS,
    CONCLUSION,
    COPERNICUS_TEMPLATE_DIR,
    FIGURE_FILES,
    FILES_TO_INLINE,
    INTRODUCTION,
    LATEX_DIR,
    LATEX_MAIN,
    MAGICC_SUMMARY,
    METADATA,
    OUTPUT_PDF,
    REFERENCES_BIB,
    REPLACEMENTS,
    SECTIONS,
)
from local.scenario_ghg_forcing_for_cmip7.magicc_jump_summary import (
    MAGICC_JUMP_SUMMARY,
)
from local.scenario_ghg_forcing_for_cmip7.value_checks import get_value_checks
from local.value_checks import run_value_checks

MANUSCRIPT_FILES = {
    "abstract": ABSTRACT,
    "introduction": INTRODUCTION,
    "sections": list(SECTIONS),
    "appendices": list(APPENDICES),
    "conclusion": CONCLUSION,
    "code_and_data_availability": CODE_AND_DATA_AVAILABILITY,
    "author_contribution": AUTHOR_CONTRIBUTION,
    "competing_interests": COMPETING_INTERESTS,
    "acknowledgements": ACKNOWLEDGEMENTS,
}
"""The manuscript's text, by the part of the template it goes into"""

COPERNICUS_TEMPLATE = [
    COPERNICUS_TEMPLATE_DIR / fn
    for fn in (CLEAN_COPERNICUS_TEMPLATE_FILENAME, *COPERNICUS_TEMPLATE_FILES)
]
"""The template files the assembly reads or copies"""

FIGURE_TAGS_IN_TEXT = get_tags_in_texts([*SECTIONS, *APPENDICES]) | {
    tag for group in APPENDIX_FIGURE_GROUPS for tag in get_generated_figure_tags(group)
}
"""Tags of the figures the manuscript includes

Those included by hand, plus those the generated appendix figures include
(which is all of each group's figures, except those included by hand).
Figures which aren't included aren't copied into the build,
so they are only products of the assembly if they are included.
"""

FIGURE_FILES_IN_TEXT = {
    tag: path for tag, path in FIGURE_FILES.items() if tag in FIGURE_TAGS_IN_TEXT
}
"""Figures the manuscript includes, by tag"""

LATEX_SUPPORT_FILES = [
    *(LATEX_DIR / fn for fn in COPERNICUS_TEMPLATE_FILES),
    LATEX_DIR / REFERENCES_BIB.name,
    *(LATEX_DIR / figure.name for figure in FIGURE_FILES_IN_TEXT.values()),
]
"""Files the assembly puts next to `main.tex`, which latex needs to compile it"""


def task_assemble_latex(  # noqa: PLR0913
    manuscript_files: dict[str, Path | list[Path]] = MANUSCRIPT_FILES,
    metadata: Path = METADATA,
    replacements: Path = REPLACEMENTS,
    references_bib: Path = REFERENCES_BIB,
    copernicus_template: list[Path] = COPERNICUS_TEMPLATE,
    figure_files: dict[str, Path] = FIGURE_FILES_IN_TEXT,
    files_to_inline: dict[str, Path] = FILES_TO_INLINE,
    value_check_data: Path = MAGICC_SUMMARY,
    magicc_jump_summary: Path = BUNDLE_DIR / MAGICC_JUMP_SUMMARY,
    value_checks_source: Path = Path(
        local.scenario_ghg_forcing_for_cmip7.value_checks.__file__
    ),
    gmd_latex_source: Path = Path(local.gmd_latex.__file__),
    main_tex: Annotated[Path, Product] = LATEX_MAIN,
    latex_support_files: Annotated[list[Path], Product] = LATEX_SUPPORT_FILES,
) -> None:
    """
    Assemble `main.tex` and check the values behind its statements

    Fails if any value check comment's statement no longer holds
    or has no check.
    """
    latex_main, collector = assemble_main_tex(
        **manuscript_files,
        replacements=replacements,
        metadata=metadata,
        references_bib_file=references_bib,
        copernicus_template_dir=COPERNICUS_TEMPLATE_DIR,
        clean_copernicus_template_filename=CLEAN_COPERNICUS_TEMPLATE_FILENAME,
        build_dir=BUILD_DIR,
        figure_files=list(figure_files.items()),
        files_to_inline=list(files_to_inline.items()),
    )
    if latex_main != main_tex:
        msg = f"Assembly wrote {latex_main}, expected {main_tex}"
        raise AssertionError(msg)

    value_check_report = run_value_checks(
        collector.specs, get_value_checks(bundle_dir=BUNDLE_DIR)
    )
    print(value_check_report.to_str())
    value_check_report.raise_if_not_ok()


def task_compile_pdf(
    main_tex: Path = LATEX_MAIN,
    latex_support_files: list[Path] = LATEX_SUPPORT_FILES,
    output_pdf: Annotated[Path, Product] = OUTPUT_PDF,
) -> None:
    """Compile `main.tex` to PDF"""
    # stderr, as in scripts/compile-gmd-template-based-latex.py
    built_pdf = compile_latex(main_tex, output_sink=sys.stderr)

    output_pdf.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(built_pdf, output_pdf)
