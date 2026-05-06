"""
Compile a GMD-templated based latex document to PDF
"""

import copy
import shutil
import subprocess
import tomllib
from pathlib import Path
from typing import Annotated, Any

import typer

REPO_ROOT = Path(__file__).parents[1]


def do_basic_replacements(
    start: str, metadata: dict[str, Any], references_bib_stem: str
) -> str:
    """
    Do basic replacements
    """
    res = copy.deepcopy(start)
    for old, new in (
        (r"\title{TEXT}", rf"\title{{{metadata['title_info']['title']}}}"),
        (
            r"\runningtitle{TEXT}",
            rf"\runningtitle{{{metadata['title_info']['running_title']}}}",
        ),
        (
            r"\documentclass[journal abbreviation, manuscript]{copernicus}",
            r"\documentclass[gmd, manuscript]{copernicus}",
        ),
        (
            r"\bibliography{example.bib}",
            rf"\bibliography{{{references_bib_stem}}}",
        ),
    ):
        res = res.replace(old, new)

    return res


def copy_template_files(template_dir: Path, build_dir: Path) -> None:
    """Copy template files into the latex build dir"""
    for fn in [
        "copernicus.bst",
        "copernicus.cls",
        "copernicus.cfg",
        "pdfscreen.sty",
        "pdfscreencop.sty",
    ]:
        shutil.copy2(template_dir / fn, build_dir / fn)


def aux_file_has_citations(aux_file: Path) -> bool:
    """Check whether BibTeX has citation entries to process"""
    if not aux_file.exists():
        return False

    return any(
        line.startswith(r"\citation")
        for line in aux_file.read_text(encoding="utf-8").splitlines()
    )


def remove_stale_latex_artifacts(latex_main: Path) -> None:
    """Remove outputs that can otherwise be reused from an earlier build."""
    for suffix in [".aux", ".bbl", ".bcf", ".blg", ".log", ".out", ".pdf", ".run.xml"]:
        latex_main.with_suffix(suffix).unlink(missing_ok=True)


def main(  # noqa: PLR0913
    metadata: Annotated[
        Path,
        typer.Option(
            help="Path to the metadata file",
            dir_okay=False,
            file_okay=True,
        ),
    ],
    references_bib_file: Annotated[
        Path, typer.Option(help="Bibtex references file to use")
    ],
    copernicus_template_dir: Annotated[
        Path,
        typer.Option(
            help="Path in which the copernicus latex template was extracted",
            dir_okay=True,
            file_okay=False,
        ),
    ],
    clean_copernicus_template_filename: Annotated[
        str,
        typer.Option(
            help=(
                "Name of the clean copernicus template file "
                "(must be in `copernicus_template_dir`)"
            )
        ),
    ],
    build_dir: Annotated[
        Path,
        typer.Option(
            help="Path in which to the build is being done",
            dir_okay=True,
            file_okay=False,
        ),
    ],
    output: Annotated[
        Path,
        typer.Option(
            help="Path in which to write the output PDF",
            dir_okay=False,
            file_okay=True,
        ),
    ],
) -> None:
    """
    Compile the PDF
    """
    with open(metadata, "rb") as fh:
        metadata_values = tomllib.load(fh)

    with open(copernicus_template_dir / clean_copernicus_template_filename) as fh:
        raw = fh.read()

    res = do_basic_replacements(
        raw, metadata_values, references_bib_stem=references_bib_file.stem
    )

    # Insert a dummy citation for now
    res = res.replace(r"\maketitle", (r"\maketitle" "\n\\citet{dunne2025evolving}"))

    latex_dir = build_dir / "latex"
    latex_dir.mkdir(exist_ok=True, parents=True)

    latex_main = latex_dir / "main.tex"
    with open(latex_main, "w") as fh:
        fh.write(res)

    remove_stale_latex_artifacts(latex_main)

    copy_template_files(template_dir=copernicus_template_dir, build_dir=latex_dir)
    shutil.copy2(references_bib_file, latex_main.parent / references_bib_file.name)

    subprocess.run(
        ["pdflatex", "-interaction=nonstopmode", "-halt-on-error", latex_main.name],
        cwd=latex_main.parent,
        check=True,
    )

    if aux_file_has_citations(latex_main.with_suffix(".aux")):
        subprocess.run(
            ["bibtex", latex_main.stem],
            cwd=latex_main.parent,
            check=True,
        )

    subprocess.run(
        ["pdflatex", "-interaction=nonstopmode", "-halt-on-error", latex_main.name],
        cwd=latex_main.parent,
        check=True,
    )
    subprocess.run(
        ["pdflatex", "-interaction=nonstopmode", "-halt-on-error", latex_main.name],
        cwd=latex_main.parent,
        check=True,
    )

    built_pdf = latex_main.parent / f"{latex_main.stem}.pdf"
    if not built_pdf.exists():
        raise FileNotFoundError(built_pdf)

    output.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(built_pdf, output)


if __name__ == "__main__":
    typer.run(main)
