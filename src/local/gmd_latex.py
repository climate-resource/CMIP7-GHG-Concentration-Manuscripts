"""
Assembly and compilation of GMD-template based latex manuscripts

Shared by every manuscript we build with Copernicus' template.
[assemble_main_tex][] puts the manuscript's pieces into the template
and checks the result, [compile_latex][] turns it into a PDF.
"""

from __future__ import annotations

import copy
import re
import shutil
import subprocess
import textwrap
import tomllib
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import IO, Any

import yaml

from local.paths import REPO_ROOT
from local.value_checks import ValueCheckCollector


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

    res = insert_after_tag(start, author_entries, "@author-start@")

    affiliations_entries = [
        get_source_file_str(metadata_file),
        *(rf"\affil[{i}]{{{address}}}" for address, i in affiliations.values()),
    ]
    res = insert_after_tag(res, affiliations_entries, "@affiliation-start@")

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


COPERNICUS_TEMPLATE_FILES = (
    "copernicus.bst",
    "copernicus.cls",
    "copernicus.cfg",
    "pdfscreen.sty",
    "pdfscreencop.sty",
)
"""Files from Copernicus' template which latex needs next to `main.tex`"""


def copy_template_files(template_dir: Path, build_dir: Path) -> None:
    """Copy template files into the latex build dir"""
    for fn in COPERNICUS_TEMPLATE_FILES:
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


def pass_figure_specs(
    tags_and_paths: list[tuple[str, Path]],
) -> tuple[FigureSpec, ...]:
    """
    Pass tags and figure paths into `FigureSpec`s
    """
    return tuple(
        FigureSpec(tag_to_replace=tag, full_path=full_path)
        for tag, full_path in tags_and_paths
    )


COMMENT_START = re.compile(r"(?<!\\)%")
"""
The start of a latex comment (a `%`, unless it is escaped)
"""


def split_comment(line: str) -> tuple[str, str]:
    """
    Split a line of latex into the part before its comment and its comment

    Parameters
    ----------
    line
        Line to split

    Returns
    -------
    :
        The part of `line` before its comment
        and its comment (empty if it has none),
        so that joining them gives back `line`
    """
    match = COMMENT_START.search(line)
    if match is None:
        return line, ""

    return line[: match.start()], line[match.start() :]


def tag_in_text(tag: str, text: str) -> bool:
    """
    Check whether a tag appears in the text, outside of comments

    Parameters
    ----------
    tag
        Tag to look for

    text
        Text to look in

    Returns
    -------
    :
        Whether `tag` appears in `text` outside of comments
    """
    return any(tag in split_comment(line)[0] for line in text.splitlines())


def replace_tag(text: str, tag: str, replacement: str) -> str:
    """
    Replace a tag in the text, except where it appears in comments

    Comments can then mention tags
    (e.g. to explain where something comes from)
    without the tag being replaced there.

    Parameters
    ----------
    text
        Text in which to replace the tag

    tag
        Tag to replace

    replacement
        What to replace `tag` with

    Returns
    -------
    :
        `text`, with `tag` replaced by `replacement` outside of comments
    """
    res_l = []
    for line in text.splitlines():
        before_comment, comment = split_comment(line)
        res_l.append(f"{before_comment.replace(tag, replacement)}{comment}")

    return "\n".join(res_l)


UNREPLACED_FIGURE_TAG = re.compile(r"\\includegraphics(\[[^\]]*\])?\{(@[^}@]*@)\}")
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
        for match in UNREPLACED_FIGURE_TAG.finditer(split_comment(line)[0])
    ]
    if unreplaced:
        msg = (
            "No figure file was given for these tags: "
            f"{', '.join(unreplaced)}. "
            "Pass them with `--figure-file` or `--tex-inputs-manifest`."
        )
        raise AssertionError(msg)


UNREPLACED_TAG = re.compile(r"@[a-z0-9][a-z0-9-]*@")
"""
Any tag, i.e. something we should have replaced
"""


def check_no_unreplaced_tags(text: str) -> None:
    """
    Check that every tag in the text has been replaced

    Without this, a tag we were never given anything for
    ends up in the PDF as text.

    Parameters
    ----------
    text
        Text to check

    Raises
    ------
    AssertionError
        `text` (outside comments) still contains tags
    """
    unreplaced = sorted(
        {
            match.group(0)
            for line in text.splitlines()
            for match in UNREPLACED_TAG.finditer(split_comment(line)[0])
        }
    )
    if unreplaced:
        msg = (
            "Nothing was given to replace these tags: "
            f"{', '.join(unreplaced)}. "
            "Pass them with `--figure-file`, `--inline-file` "
            "or `--tex-inputs-manifest`."
        )
        raise AssertionError(msg)


NOT_PROSE = re.compile(
    r"\$[^$]*\$"
    r"|\\(url|texttt|cite\w*|ref|label|href|includegraphics|input)"
    r"(\[[^\]]*\])?\{[^}]*\}"
)
"""
Parts of a line which aren't prose (maths, references, URLs, code, file names)
"""

DISPLAY_MATHS_ENVIRONMENT = re.compile(r"\\(begin|end)\{(equation|align|eqnarray)\*?\}")
"""
The start or end of a display maths environment
"""

FULL_WORD_REFERENCE = re.compile(r"\b(Figures?|Sections?|Equations?)~\(?\\ref\b")
"""
A reference which uses the full word (only right at the start of a sentence)
"""

ABBREVIATED_REFERENCE = re.compile(r"\b(Figs?|Sects?|Eqs?)\.~\(?\\ref\b")
"""
A reference which uses the abbreviation (everywhere except the start of a sentence)
"""

EQUATION_REFERENCE_WITHOUT_PARENTHESES = re.compile(r"\b(Eqs?\.|Equations?)~\\ref\b")
"""
A reference to an equation whose number isn't in parentheses
"""

SPACED_HYPHEN = re.compile(r"\S - \S")
"""
A hyphen used as a dash (a dash is `--`, spaced for a break in a sentence)
"""

HYPHENATED_YEAR_RANGE = re.compile(r"\b(1[0-9]{3}|20[0-9]{2})-(1[0-9]{3}|20[0-9]{2})\b")
"""
A range of years written with a hyphen (a range is written with `--`)
"""

NOT_SENTENCE_ENDS = ("e.g.", "i.e.", "et al.", "cf.", "Fig.", "Figs.", "Eq.", "Eqs.")
"""
Text which ends in a full stop without ending a sentence
"""


def is_sentence_start(before: str) -> bool:
    """
    Check whether the text that comes before something means it starts a sentence

    Parameters
    ----------
    before
        Prose before the thing of interest,
        empty if it starts a paragraph

    Returns
    -------
    :
        Whether the thing of interest starts a sentence
    """
    before = before.rstrip()
    if not before or before.endswith(("{", r"\item")):
        return True

    return before.endswith((".", "?", "!")) and not before.endswith(NOT_SENTENCE_ENDS)


def get_reference_style_problems(code: str, previous: str) -> list[str]:
    """
    Get the ways in which a line's references don't follow the Copernicus style

    Parameters
    ----------
    code
        Line to check, without its comment

    previous
        The previous line of prose, empty if `code` starts a paragraph

    Returns
    -------
    :
        Description of each problem
    """
    problems = []
    for match in FULL_WORD_REFERENCE.finditer(code):
        if not is_sentence_start(code[: match.start()] or previous):
            problems.append(f"{match.group(1)} mid-sentence")

    for match in ABBREVIATED_REFERENCE.finditer(code):
        if is_sentence_start(code[: match.start()] or previous):
            problems.append(f"{match.group(1)}. starts a sentence")

    if EQUATION_REFERENCE_WITHOUT_PARENTHESES.search(code):
        problems.append("equation number without parentheses")

    return problems


def check_copernicus_style(text: str) -> None:
    """
    Check that references and dashes follow the Copernicus style

    The rules are from Copernicus' manuscript preparation guidelines
    (https://publications.copernicus.org/for_authors/manuscript_preparation.html):

    - "Fig.", "Sect." and "Eq." in running text,
      "Figure", "Section" and "Equation" at the start of a sentence
    - equation numbers in parentheses, e.g. "Eq. (1)"
    - ranges with an unspaced en dash (`--`), e.g. "2000--2010"
    - spaced en dashes (` -- `) for breaks in a sentence, not hyphens

    Maths, comments, references, URLs and code are not checked.

    Parameters
    ----------
    text
        Text to check

    Raises
    ------
    AssertionError
        `text` doesn't follow the Copernicus style
    """
    problems = []
    previous = ""
    in_display_maths = False
    for line in text.splitlines():
        code = split_comment(line)[0]
        if not line.strip():
            # A new paragraph
            previous = ""
            continue

        if not code.strip():
            # A comment line, which doesn't end a paragraph
            continue

        environment = DISPLAY_MATHS_ENVIRONMENT.search(code)
        if environment is not None:
            in_display_maths = environment.group(1) == "begin"
            previous = code
            continue

        if in_display_maths:
            continue

        problems.extend(
            f"{problem}: {line.strip()}"
            for problem in get_reference_style_problems(code, previous)
        )

        prose = NOT_PROSE.sub(" ", code)
        if SPACED_HYPHEN.search(prose):
            problems.append(f"hyphen used as a dash: {line.strip()}")

        if HYPHENATED_YEAR_RANGE.search(prose):
            problems.append(f"range of years with a hyphen: {line.strip()}")

        previous = code

    if problems:
        msg = "The text doesn't follow the Copernicus style:\n" + "\n".join(
            f"- {p}" for p in problems
        )
        raise AssertionError(msg)


def inline_files(
    in_text: str, files: list[tuple[str, Path]], collector: ValueCheckCollector
) -> str:
    """
    Inline files in place of their tags

    The files' content is put straight into the text,
    rather than being left for latex to pull in with an input command,
    so the output stays a single file (which is what Copernicus wants)
    and our replacements are applied to the files' content too.

    Parameters
    ----------
    in_text
        Text in which to inline the files

    files
        Tag to replace and the file to replace it with, for each file

    collector
        Collector of value check comments, used to read the files

    Returns
    -------
    :
        `in_text`, with each tag replaced by its file's content
    """
    res = in_text
    for tag, full_path in files:
        if not tag_in_text(tag, res):
            msg = f"Did not find {tag=} in the text (outside of comments)"
            raise AssertionError(msg)

        res = replace_tag(
            res,
            tag,
            f"{get_source_file_str(full_path)}\n{collector.read_text(full_path)}",
        )

    return res


def assemble_main_tex(  # noqa: PLR0913, PLR0915
    *,
    abstract: Path,
    introduction: Path,
    sections: Sequence[Path],
    conclusion: Path,
    code_and_data_availability: Path,
    author_contribution: Path,
    competing_interests: Path,
    acknowledgements: Path,
    replacements: Path,
    metadata: Path,
    references_bib_file: Path,
    copernicus_template_dir: Path,
    clean_copernicus_template_filename: str,
    build_dir: Path,
    appendices: Sequence[Path] = (),
    extra: Sequence[Path] = (),
    auxiliary: Sequence[Path] = (),
    figure_files: Sequence[tuple[str, Path]] = (),
    files_to_inline: Sequence[tuple[str, Path]] = (),
    files_to_inline_if_tagged: Sequence[tuple[str, Path]] = (),
) -> tuple[Path, ValueCheckCollector]:
    r"""
    Assemble the manuscript's `main.tex` in the build directory

    The pieces are put into Copernicus' template,
    the figures and the files to inline go in place of their tags,
    the result is checked (no tags left over, Copernicus' style)
    and then our replacements are applied.
    Everything latex needs to compile `main.tex`
    (the template's class and style files, the bibliography, the figures)
    is copied next to it.

    Parameters
    ----------
    abstract
        Abstract (no headers, they come from the template)

    introduction
        Introduction (no headers, they come from the template)

    sections
        Body sections, in order (with their own headers)

    conclusion
        Conclusion (no headers, they come from the template)

    code_and_data_availability
        Code and data availability statement

    author_contribution
        Author contribution statement

    competing_interests
        Competing interests statement

    acknowledgements
        Acknowledgements

    replacements
        YAML file of replacements to apply to the assembled latex

        Allows us to use short-hand without annoying Copernicus,
        who don't want us to define special commands in the latex we give them.

    metadata
        The manuscript's metadata file (title, authors, affiliations)

    references_bib_file
        Bibtex references file

    copernicus_template_dir
        Directory in which Copernicus' latex template was extracted

    clean_copernicus_template_filename
        Name of the clean template file (in `copernicus_template_dir`)

    build_dir
        Directory in which to build (`main.tex` goes in its `latex` sub-directory)

    appendices
        Appendices, in order (each `\section` is a new appendix: A, B, ...)

        The `\appendix` and `\noappendix` commands are added here,
        so the files must not contain them.

    extra
        Files to copy into the build directory, with replacements applied,
        but not include in `main.tex` (for latex's `\input` or `\include`)

    auxiliary
        Files to copy into the build directory without any processing

    figure_files
        Tag to replace and the figure to replace it with, for each figure

        Figures whose tag is not in the text are skipped.

    files_to_inline
        Tag to replace and the file to inline in its place, for each file

        Every tag must be in the text.

    files_to_inline_if_tagged
        As `files_to_inline`, except files whose tag is not in the text are skipped

    Returns
    -------
    :
        Path to `main.tex` and the collector of the value check comments
        in everything which went into it

    Raises
    ------
    AssertionError
        Two figures would land at the same path in the build directory
    """
    collector = ValueCheckCollector(root=REPO_ROOT)
    figure_specs = pass_figure_specs(list(figure_files))

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
        res, abstract, tag="@abstract-start@", collector=collector
    )
    res = insert_file_content_after_tag(
        res, introduction, tag="@introduction-start@", collector=collector
    )

    body_text = "\n\n".join(
        f"{get_source_file_str(sf)}\n{collector.read_text(sf)}" for sf in sections
    )
    res = insert_after_tag(res, body_text, tag="@body-start@")

    res = insert_file_content_after_tag(
        res, conclusion, tag="@conclusions-start@", collector=collector
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

    if appendices:
        appendix_text = "\n\n".join(
            [
                r"\appendix",
                *(
                    f"{get_source_file_str(af)}\n{collector.read_text(af)}"
                    for af in appendices
                ),
                # Otherwise the sections and figures after the appendices
                # keep the appendix numbering
                r"\noappendix",
            ]
        )
        res = insert_after_tag(res, appendix_text, tag="@appendix-start@")

    res = insert_file_content_after_tag(
        res, acknowledgements, tag="@acknowledgements-start@", collector=collector
    )

    # Before any other replacements,
    # so the inlined files get the same replacements as the rest of the text
    # (including the figure replacements, as they can include figures).
    res = inline_files(res, list(files_to_inline), collector=collector)
    res = inline_files(
        res,
        [
            (tag, full_path)
            for tag, full_path in files_to_inline_if_tagged
            if tag_in_text(tag, res)
        ],
        collector=collector,
    )

    latex_dir = build_dir / "latex"
    latex_dir.mkdir(exist_ok=True, parents=True)

    seen = set()
    for fs in figure_specs:
        # Otherwise every figure we are given ends up in the build,
        # whether the text uses it or not.
        if not tag_in_text(fs.tag_to_replace, res):
            continue

        dest = latex_dir / fs.full_path.name
        dest.parent.mkdir(exist_ok=True, parents=True)
        if dest in seen:
            msg = f"Multiple figure files will land at {dest!r}"
            raise AssertionError(msg)

        shutil.copy2(fs.full_path, dest)
        seen.add(dest)
        res = replace_tag(res, fs.tag_to_replace, str(dest.relative_to(latex_dir)))

    # The figure check first, as its message is more specific
    check_no_unreplaced_figure_tags(res)
    check_no_unreplaced_tags(res)
    check_copernicus_style(res)

    replacements_map = yaml.safe_load(replacements.read_text())
    res = apply_replacements(res, replacements_map)

    latex_main = latex_dir / "main.tex"
    with open(latex_main, "w") as fh:
        fh.write(res)

    remove_stale_latex_artifacts(latex_main)

    copy_template_files(template_dir=copernicus_template_dir, build_dir=latex_dir)
    shutil.copy2(references_bib_file, latex_main.parent / references_bib_file.name)

    for extra_file in extra:
        raw = collector.read_text(extra_file)
        mapped = apply_replacements(raw, replacements_map)
        (latex_dir / extra_file.name).write_text(mapped)

    for auxiliary_file in auxiliary:
        shutil.copy2(auxiliary_file, latex_main.parent / auxiliary_file.name)

    return latex_main, collector
