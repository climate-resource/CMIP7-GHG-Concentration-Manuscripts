"""
Compile a GMD-templated based latex document to PDF
"""

import copy
import json
import re
import shutil
import subprocess
import sys
import textwrap
import tomllib
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import IO, Annotated, Any

import typer
import yaml

from local.historical_ghg_forcing_for_cmip7.value_checks import VALUE_CHECKS
from local.value_checks import ValueCheckCollector, run_value_checks

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


def get_source_file_str(source_file: Path) -> str:
    """Get a string to insert in the compiled latex so we can see the original source"""
    return f"% Source file: {source_file.relative_to(REPO_ROOT)}"


def insert_author_list_and_affiliations(
    start: str, metadata: dict[str, Any], metadata_file: Path
) -> str:
    """Insert author list into the text"""
    affiliations = {
        key: (value, i + 1)
        for i, (key, value) in enumerate(metadata["affiliations"].items())
    }

    author_entries = [get_source_file_str(metadata_file)]
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
        get_source_file_str(metadata_file),
        *(rf"\affil[{i}]{{{address}}}" for address, i in affiliations.values()),
    ]
    res = insert_after_tag(res, affiliations_entries, "<affiliation-start>")

    return res


def insert_file_content_after_tag(
    in_text: str, filepath: Path, tag: str, collector: ValueCheckCollector
) -> str:
    """
    Insert file content after a specific tag
    """
    to_insert = f"{get_source_file_str(filepath)}\n{collector.read_text(filepath)}"

    res = insert_after_tag(in_text, to_insert, tag=tag)

    return res


def apply_replacements(in_text: str, replacements: Mapping[str, str]) -> str:
    """
    Apply replacements from our replacements file
    """
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


def run_pdflatex(latex_main: Path, output_sink: int | IO[Any] | None) -> None:
    """
    Run a single pdflatex pass over `latex_main`

    Parameters
    ----------
    latex_main
        Main latex file to compile

    output_sink
        Where to send pdflatex's output, see [compile_latex][]
    """
    subprocess.run(
        ["pdflatex", "-interaction=nonstopmode", "-halt-on-error", latex_main.name],
        cwd=latex_main.parent,
        check=True,
        stdout=output_sink,
    )


def compile_latex(
    latex_main: Path,
    output_sink: int | IO[Any] | None,
    n_passes_after_bibtex: int = 3,
) -> Path:
    """
    Compile `latex_main` to PDF

    Parameters
    ----------
    latex_main
        Main latex file to compile

    output_sink
        Where to send the output of pdflatex and bibtex.

        They write everything, errors and warnings included,
        to stdout (as well as their log files).
        Anything [subprocess.run][]'s `stdout` accepts can be used, e.g.
        `sys.stderr` (so the output can be dropped with `2>/dev/null`
        without losing the rest of the script's output),
        `subprocess.DEVNULL` (to drop it)
        or `None` (to leave it on stdout).

    n_passes_after_bibtex
        Number of pdflatex passes to do after running bibtex.
        Two isn't enough for this manuscript:
        cross-references still move on the second pass.

    Returns
    -------
    :
        Path to the compiled PDF
    """
    run_pdflatex(latex_main, output_sink=output_sink)

    if aux_file_has_citations(latex_main.with_suffix(".aux")):
        subprocess.run(
            ["bibtex", latex_main.stem],
            cwd=latex_main.parent,
            check=True,
            stdout=output_sink,
        )

    for _ in range(n_passes_after_bibtex):
        run_pdflatex(latex_main, output_sink=output_sink)

    log = latex_main.with_suffix(".log").read_text(encoding="utf-8", errors="replace")
    if "Rerun to get cross-references right" in log:
        print(
            "WARNING: latex still wants another pass "
            "to get the cross-references right. "
            "Consider increasing `n_passes_after_bibtex`."
        )

    built_pdf = latex_main.with_suffix(".pdf")
    if not built_pdf.exists():
        raise FileNotFoundError(built_pdf)

    return built_pdf


@dataclass
class FigureSpec:
    """
    Spec for including a figure in the output latex
    """

    full_path: Path
    """
    Full path to the figure
    """

    tag_to_replace: str
    """
    Tag in the raw latex to replace with the figure's path (relative to the build path)
    """


def split_tag_and_path(value: str, option: str) -> tuple[str, Path]:
    """
    Split a CLI value of the form `tag-to-replace=full-path`

    Parameters
    ----------
    value
        Value to split

    option
        CLI option the value was passed to, for the error message

    Returns
    -------
    :
        The tag to replace and the path
    """
    tag, _, full_path = value.partition("=")
    if not full_path:
        msg = (
            f"Bad value for {option}. "
            "Please format as `tag-to-replace=full-path`. "
            f"Received: {value!r}"
        )
        raise ValueError(msg)

    return tag, Path(full_path)


def pass_figure_specs(raw: list[str]) -> tuple[FigureSpec, ...]:
    """
    Pass CLI figure information into `FigureSpec`s
    """
    res_l = []
    for v in raw:
        tag, full_path = split_tag_and_path(v, option="--figure-file")
        fs = FigureSpec(tag_to_replace=tag, full_path=full_path)
        res_l.append(fs)

    res = tuple(res_l)

    return res


def load_tex_inputs_manifest(manifest_file: Path) -> tuple[list[str], list[str]]:
    """
    Load a manifest of figure and table files

    Parameters
    ----------
    manifest_file
        JSON file with a `figures` and a `tables` entry,
        each of which maps a tag to replace to a full path

    Returns
    -------
    :
        The figures and tables, formatted as for `--figure-file` and `--table-file`
    """
    manifest = json.loads(manifest_file.read_text())

    figure_files = [f"{tag}={path}" for tag, path in manifest["figures"].items()]
    table_files = [f"{tag}={path}" for tag, path in manifest["tables"].items()]

    return figure_files, table_files


def tag_is_used(value: str, text: str) -> bool:
    """
    Check whether the tag in a `tag-to-replace=full-path` value appears in the text
    """
    tag, _, _ = value.partition("=")

    return tag in text


UNREPLACED_FIGURE_TAG = re.compile(r"\\includegraphics(\[[^\]]*\])?\{(<[^}]*>)\}")
"""
A figure which is still a tag, i.e. one we were never given a file for
"""


def check_no_unreplaced_figure_tags(text: str) -> None:
    """
    Check that every figure included in the text has been given a file

    Parameters
    ----------
    text
        Text to check

    Raises
    ------
    AssertionError
        A figure included in `text` (outside a comment) is still a tag
    """
    unreplaced = [
        match.group(2)
        for line in text.splitlines()
        if not line.lstrip().startswith("%")
        for match in UNREPLACED_FIGURE_TAG.finditer(line)
    ]
    if unreplaced:
        msg = (
            "No figure file was given for these tags: "
            f"{', '.join(unreplaced)}. "
            "Pass them with `--figure-file` or `--tex-inputs-manifest`."
        )
        raise AssertionError(msg)


def inline_table_files(
    in_text: str, table_files: list[str], collector: ValueCheckCollector
) -> str:
    """
    Inline table files in place of their tags

    The files' content is put straight into the text,
    rather than being left for latex to pull in with an input command,
    so the output stays a single file (which is what Copernicus wants)
    and our replacements are applied to the tables too.

    Parameters
    ----------
    in_text
        Text in which to inline the tables

    table_files
        Values passed to `--table-file`, as `tag-to-replace=full-path`

    collector
        Collector of value check comments, used to read the files

    Returns
    -------
    :
        `in_text`, with each tag replaced by its file's content
    """
    res = in_text
    for v in table_files:
        tag, full_path = split_tag_and_path(v, option="--table-file")
        if tag not in res:
            msg = f"Did not find {tag=} in the text"
            raise AssertionError(msg)

        res = res.replace(
            tag, f"{get_source_file_str(full_path)}\n{collector.read_text(full_path)}"
        )

    return res


def main(  # noqa: PLR0913, PLR0915
    abstract: Annotated[
        Path,
        typer.Option(
            help=(
                "Path to the abstract file. "
                "It will be used as-is. "
                "The file must not contain any headers, they come from the template."
            ),
            dir_okay=False,
            file_okay=True,
        ),
    ],
    introduction: Annotated[
        Path,
        typer.Option(
            help=(
                "Path to the introduction file. "
                "It will be used as-is. "
                "The file must not contain any headers, they come from the template."
            ),
            dir_okay=False,
            file_okay=True,
        ),
    ],
    section: Annotated[
        list[Path],
        typer.Option(
            help=(
                "Path to use for the body sections. "
                "This can be supplied multiple times. "
                "The files are used in the order they are provided to the CLI. "
                "The files are used as-is and should contain headers as needed. "
            ),
            dir_okay=False,
            file_okay=True,
        ),
    ],
    conclusion: Annotated[
        Path,
        typer.Option(
            help=(
                "Path to the conclusion file. "
                "It will be used as-is. "
                "The file must not contain any headers, they come from the template."
            ),
            dir_okay=False,
            file_okay=True,
        ),
    ],
    code_and_data_availability: Annotated[
        Path,
        typer.Option(
            help=(
                "Path to the code and data availability file. "
                "It will be used as-is. "
                "The file must not contain any headers, they come from the template."
            ),
            dir_okay=False,
            file_okay=True,
        ),
    ],
    author_contribution: Annotated[
        Path,
        typer.Option(
            help=(
                "Path to the author contribution file. "
                "It will be used as-is. "
                "The file must not contain any headers, they come from the template."
            ),
            dir_okay=False,
            file_okay=True,
        ),
    ],
    competing_interests: Annotated[
        Path,
        typer.Option(
            help=(
                "Path to the competing interests file. "
                "It will be used as-is. "
                "The file must not contain any headers, they come from the template."
            ),
            dir_okay=False,
            file_okay=True,
        ),
    ],
    acknowledgements: Annotated[
        Path,
        typer.Option(
            help=(
                "Path to the acknowledgements file. "
                "It will be used as-is. "
                "The file must not contain any headers, they come from the template."
            ),
            dir_okay=False,
            file_okay=True,
        ),
    ],
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
    appendix: Annotated[
        list[Path] | None,
        typer.Option(
            help=(
                "Path to use for the appendices. "
                "This can be supplied multiple times. "
                "The files are used in the order they are provided to the CLI. "
                "As for `--section`, the files should contain headers as needed "
                r"(each `\section` is a new appendix: A, B, ...). "
                r"The `\appendix` and `\noappendix` commands come from this script, "
                "so the files must not contain them."
            ),
            dir_okay=False,
            file_okay=True,
        ),
    ] = None,
    extra: Annotated[
        list[Path] | None,
        typer.Option(
            help=(
                "Paths to copy into the build directory, "
                "but not include in the `main.tex` file. "
                "These are useful "
                r"if you use \input or \include commands in your latex "
                "(so the files need to be in the build directory, "
                "but it is left to latex to add the content). "
                "Replacements are applied to these files "
                "as part of the copying process."
            )
        ),
    ] = None,
    auxiliary: Annotated[
        list[Path] | None,
        typer.Option(
            help=(
                "Paths to copy into the build directory "
                "without any processing e.g. figure files."
            )
        ),
    ] = None,
    figure_file: Annotated[
        list[str] | None,
        typer.Option(
            help=(
                "Figure file to add. "
                "Should be passed as `tag-to-replace-in-latex=path-to-figure`, "
                "e.g. `--figure-file=<n2o-methods-figure>=/path/to/figure.pdf`"
            )
        ),
    ] = None,
    table_file: Annotated[
        list[str] | None,
        typer.Option(
            help=(
                "Table file to inline. "
                "Should be passed as `tag-to-replace-in-latex=path-to-table`, "
                "e.g. `--table-file=<per-gas-table>=/path/to/table.tex`. "
                "The tag is replaced by the file's content "
                "(so it can hold any latex, not only a table)."
            )
        ),
    ] = None,
    tex_inputs_manifest: Annotated[
        Path | None,
        typer.Option(
            help=(
                "JSON file of figure and table files to add. "
                "It should have a `figures` and a `tables` entry, "
                "each of which maps a tag to replace in the latex "
                "to the path of the file, "
                "i.e. the same information as `--figure-file` and `--table-file`. "
                "Unlike those options, entries whose tag is not in the text "
                "are skipped, so the manifest can list more than is used."
            ),
            dir_okay=False,
            file_okay=True,
        ),
    ] = None,
    check_values: Annotated[
        bool,
        typer.Option(
            help=(
                "Check the values behind the `% value-check: {...}` comments "
                "in the latex before compiling. "
                "Fails if any statement no longer holds "
                "or if any comment has no check."
            )
        ),
    ] = True,
) -> None:
    """
    Compile the PDF
    """
    collector = ValueCheckCollector(root=REPO_ROOT)
    if tex_inputs_manifest is not None:
        manifest_figure_files, manifest_table_files = load_tex_inputs_manifest(
            tex_inputs_manifest
        )
    else:
        manifest_figure_files, manifest_table_files = [], []

    figure_specs = pass_figure_specs([*(figure_file or []), *manifest_figure_files])

    with open(metadata, "rb") as fh:
        metadata_values = tomllib.load(fh)

    with open(copernicus_template_dir / clean_copernicus_template_filename) as fh:
        raw = fh.read()

    res = do_basic_replacements(
        raw, metadata_values, references_bib_stem=references_bib_file.stem
    )

    res = insert_author_list_and_affiliations(
        res, metadata_values, metadata_file=metadata
    )

    res = insert_file_content_after_tag(
        res, abstract, tag="<abstract-start>", collector=collector
    )
    res = insert_file_content_after_tag(
        res, introduction, tag="<introduction-start>", collector=collector
    )

    body_text = "\n\n".join(
        f"{get_source_file_str(sf)}\n{collector.read_text(sf)}" for sf in section
    )
    res = insert_after_tag(res, body_text, tag="<body-start>")

    res = insert_file_content_after_tag(
        res, conclusion, tag="<conclusions-start>", collector=collector
    )

    for start_code, source_file in (
        ("codedataavailability", code_and_data_availability),
        ("authorcontribution", author_contribution),
        ("competinginterests", competing_interests),
    ):
        to_replace = rf"\{start_code}{{TEXT}}"
        source_text = collector.read_text(source_file)
        replacement_text = f"{get_source_file_str(source_file)}\n{source_text}"
        replacement = to_replace.replace(
            "TEXT", f"\n{textwrap.indent(replacement_text, prefix=4 * ' ')}\n"
        )
        res = res.replace(to_replace, replacement)

    if appendix:
        appendix_text = "\n\n".join(
            [
                r"\appendix",
                *(
                    f"{get_source_file_str(af)}\n{collector.read_text(af)}"
                    for af in appendix
                ),
                # Otherwise the sections and figures after the appendices
                # keep the appendix numbering
                r"\noappendix",
            ]
        )
        res = insert_after_tag(res, appendix_text, tag="<appendix-start>")

    res = insert_file_content_after_tag(
        res, acknowledgements, tag="<acknowledgements-start>", collector=collector
    )

    # Before any other replacements,
    # so the tables get the same replacements as the rest of the text.
    if table_file is not None:
        res = inline_table_files(res, table_file, collector=collector)

    res = inline_table_files(
        res,
        [v for v in manifest_table_files if tag_is_used(v, res)],
        collector=collector,
    )

    latex_dir = build_dir / "latex"
    latex_dir.mkdir(exist_ok=True, parents=True)

    figure_replacements = {}
    seen = set()
    for fs in figure_specs:
        # Otherwise every figure we are given ends up in the build,
        # whether the text uses it or not.
        if fs.tag_to_replace not in res:
            continue

        dest = latex_dir / fs.full_path.name
        dest.parent.mkdir(exist_ok=True, parents=True)
        if dest in seen:
            msg = f"Multiple figure files will land at {dest!r}"
            raise AssertionError(msg)

        shutil.copy2(fs.full_path, dest)
        seen.add(dest)
        figure_replacements[fs.tag_to_replace] = str(dest.relative_to(latex_dir))

    res = apply_replacements(res, figure_replacements)

    check_no_unreplaced_figure_tags(res)

    replacements_map = yaml.safe_load(replacements.read_text())
    res = apply_replacements(res, replacements_map)

    latex_main = latex_dir / "main.tex"
    with open(latex_main, "w") as fh:
        fh.write(res)

    remove_stale_latex_artifacts(latex_main)

    copy_template_files(template_dir=copernicus_template_dir, build_dir=latex_dir)
    shutil.copy2(references_bib_file, latex_main.parent / references_bib_file.name)

    for extra_file in extra if extra is not None else []:
        raw = collector.read_text(extra_file)
        mapped = apply_replacements(raw, replacements_map)
        (latex_dir / extra_file.name).write_text(mapped)

    for auxiliary_file in auxiliary if auxiliary is not None else []:
        shutil.copy2(auxiliary_file, latex_main.parent / auxiliary_file.name)

    # Last thing before building, so every file has been read
    # (and its value check comments collected)
    if check_values:
        value_check_report = run_value_checks(collector.specs, VALUE_CHECKS)
        print(value_check_report.to_str())
        value_check_report.raise_if_not_ok()

    # stderr, so the latex output can be dropped with `2>/dev/null`
    built_pdf = compile_latex(latex_main, output_sink=sys.stderr)

    output.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(built_pdf, output)


if __name__ == "__main__":
    # Without locals, so a failed value check's message isn't buried
    app = typer.Typer(pretty_exceptions_show_locals=False)
    app.command()(main)
    app()
