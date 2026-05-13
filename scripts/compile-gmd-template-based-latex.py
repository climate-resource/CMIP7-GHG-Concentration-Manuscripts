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
import yaml

REPO_ROOT = Path(__file__).parents[1]


def do_basic_replacements(
    start: str, metadata: dict[str, Any], references_bib_stem: str
) -> str:
    """
    Do basic replacements
    """
    res = copy.deepcopy(start)
    for old, new in (
        (
            r"\documentclass[journal abbreviation, manuscript]{copernicus}",
            r"\documentclass[gmd, manuscript]{copernicus}",
        ),
        (r"\title{TEXT}", rf"\title{{{metadata['title_info']['title']}}}"),
        (
            r"\runningtitle{TEXT}",
            rf"\runningtitle{{{metadata['title_info']['running_title']}}}",
        ),
        (
            r"\runningauthor{TEXT}",
            rf"\runningauthor{{{metadata['authors'][0]['surname']} et al.}}",
        ),
        (
            r"\bibliography{example.bib}",
            rf"\bibliography{{{references_bib_stem}}}",
        ),
    ):
        res = res.replace(old, new)

    return res


def insert_after_tag(
    start: str, to_insert: str | list[str] | tuple[str, ...], tag: str
) -> str:
    """
    Insert text after a given tag
    """
    if isinstance(to_insert, str):
        to_insert = to_insert.splitlines()

    res_l = []
    found_tag = False
    for line in start.splitlines():
        res_l.append(line)
        if tag in line:
            res_l.extend(to_insert)
            found_tag = True

    if not found_tag:
        msg = f"Did not find {tag=} in the text"
        raise AssertionError(msg)

    return "\n".join(res_l)


def insert_author_list_and_affiliations(start: str, metadata: dict[str, Any]) -> str:
    """Insert author list into the text"""
    affiliations = {
        key: (value, i + 1)
        for i, (key, value) in enumerate(metadata["affiliations"].items())
    }

    author_entries = []
    for author in metadata["authors"]:
        author_affiliations = ",".join(
            str(affiliations[key][1]) for key in author["affiliations"]
        )
        if author.get("correspondence_author", False):
            author_text = rf"\Author[{author_affiliations}][{author['email']}]{{{author['given_name']}}}{{{author['surname']}}}"  # noqa: E501
        else:
            author_text = rf"\Author[{author_affiliations}]{{{author['given_name']}}}{{{author['surname']}}}"  # noqa: E501

        author_entries.append(author_text)

    res = insert_after_tag(start, author_entries, "<author-start>")

    affiliations_entries = [
        rf"\affil[{i}]{{{address}}}" for address, i in affiliations.values()
    ]
    res = insert_after_tag(res, affiliations_entries, "<affiliation-start>")
    # breakpoint()

    return res


def apply_replacements(in_text: str, replacements_file: Path) -> str:
    """
    Apply replacements from our replacements file
    """
    replacements = yaml.safe_load(replacements_file.read_text())

    out = copy.deepcopy(in_text)
    for old, new in replacements.items():
        out = out.replace(old, new)

    return out


def remove_stale_latex_artifacts(latex_main: Path) -> None:
    """Remove outputs that can otherwise be reused from an earlier build."""
    for suffix in [".aux", ".bbl", ".bcf", ".blg", ".log", ".out", ".pdf", ".run.xml"]:
        latex_main.with_suffix(suffix).unlink(missing_ok=True)


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


def main(  # noqa: PLR0913
    abstract: Annotated[
        Path,
        typer.Option(
            help=(
                "Path to the abstract file. "
                "It will be used as-is so this file must not contain any headers."
            ),
            dir_okay=False,
            file_okay=True,
        ),
    ],
    # introduction file - must not contain any headers
    # sections - can be more than one, should use headers, used in order given
    # conclusions file - must not contain any headers
    # code and data availability - must not contain any headers
    # author contributions - must not contain any headers
    # competing interests - must not contain any headers
    # acknowledgements - must not contain any headers
    # replacements file - allows us to do last minute replacements
    replacements: Annotated[
        Path,
        typer.Option(
            help=(
                "Path to a yaml file which defines replacements to apply "
                "before compiling the latex. "
                "Allows us to use short-hand without annoying copernicus, "
                "who don't want us to define special commands "
                "in the latex we give them."
            ),
            dir_okay=False,
            file_okay=True,
        ),
    ],
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

    res = insert_author_list_and_affiliations(res, metadata_values)

    # Dummy content for now
    res = insert_after_tag(res, abstract.read_text(), tag="<abstract-start>")
    res = insert_after_tag(
        res, r"As discussed in \citet{dunne2025evolving}", tag="<introduction-start>"
    )
    res = insert_after_tag(
        res,
        r"""\section{Methods}

Methods here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.


Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.

Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.

Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.


Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.

Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.

Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.


Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.

Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.

Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.


Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.

Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.

Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.


Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.

Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
Line here.
""",
        tag="<body-start>",
    )
    res = insert_after_tag(
        res, r"As discussed in \citet{dunne2025evolving}", tag="<conclusions-start>"
    )

    res = res.replace(
        r"\codedataavailability{TEXT}",
        r"\codedataavailability{Code and data availability info}",
    )
    res = res.replace(
        r"\authorcontribution{TEXT}",
        r"\authorcontribution{Author contribution descriptions}",
    )
    res = res.replace(
        r"\competinginterests{TEXT}",
        r"\competinginterests{No competing interests}",
    )

    res = insert_after_tag(
        res, "Acknowledgements go here", tag="<acknowledgements-start>"
    )

    res = apply_replacements(res, replacements)

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
