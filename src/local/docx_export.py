"""
Export of a compiled GMD-template based latex manuscript to Word

The aim is a copy which is good enough for collaborators to review
(e.g. after uploading it to Google Docs), not a polished conversion.
The title, authors and abstract are dropped.

The work is done by pandoc.
Before handing over, we rewrite the bits of the latex pandoc doesn't understand
(Copernicus' macros, PDF figures)
and use the numbers from latex's `.aux` file for headings, references and captions,
so that e.g. "Fig. 3" means the same thing in the Word copy as in the PDF.
This means that the PDF must be compiled before exporting.
"""

from __future__ import annotations

import re
import shutil
import subprocess
import sys
from pathlib import Path

import pypdfium2

from local.gmd_latex import split_comment

COPERNICUS_SECTION_MACROS = {
    "codeavailability": "Code availability",
    "dataavailability": "Data availability",
    "codedataavailability": "Code and data availability",
    "sampleavailability": "Sample availability",
    "videosupplement": "Video supplement",
    "authorcontribution": "Author contributions",
    "competinginterests": "Competing interests",
    "disclaimer": "Disclaimer",
    "copyrightstatement": "Copyright statement",
}
"""Copernicus macros which take their section's text as an argument, and its heading"""

DROPPED_COMMANDS = ("appendix", "noappendix", "clearpage")
"""Commands (without arguments) which mean nothing in Word"""

SECTIONING_COMMAND = re.compile(
    r"\\(section|subsection|subsubsection)(\*?)(?![a-zA-Z])"
)
"""A sectioning command, its level and whether it is starred"""

NEWLABEL = re.compile(r"^\\newlabel\{([^}]*)\}\{\{([^}]*)\}")
"""A label and its number in a `.aux` file"""

TOC_ENTRY_START = "\\@writefile{toc}{\\contentsline {"
"""The start of a table of contents entry in a `.aux` file"""

FLOAT_ENVIRONMENT = re.compile(
    r"\\begin\{(figure|table)(\*?)\}(?P<body>.*?)\\end\{\1\2\}", flags=re.DOTALL
)
"""A figure or table environment"""

INCLUDEGRAPHICS = re.compile(
    r"\\includegraphics(?P<options>\[[^\]]*\])?\{(?P<path>[^}]*)\}"
)
"""An `\\includegraphics` command"""


def find_group_end(text: str, start: int) -> int:
    """
    Find the end of a brace group

    Parameters
    ----------
    text
        Text to search

    start
        Index of the group's opening brace

    Returns
    -------
    :
        Index just after the group's closing brace
    """
    if text[start] != "{":
        msg = (
            f"Expected an opening brace at {start}, found {text[start : start + 20]!r}"
        )
        raise ValueError(msg)

    depth = 0
    i = start
    while i < len(text):
        char = text[i]
        if char == "\\":
            # Skip escaped characters, e.g. `\{`
            i += 2
            continue

        if char == "{":
            depth += 1
        elif char == "}":
            depth -= 1
            if depth == 0:
                return i + 1

        i += 1

    msg = f"Unbalanced braces in group starting with {text[start : start + 50]!r}"
    raise ValueError(msg)


def strip_comments(text: str) -> str:
    """
    Strip comments from latex

    Lines which are only a comment are removed entirely,
    so they don't split paragraphs.

    Parameters
    ----------
    text
        Latex to strip

    Returns
    -------
    :
        `text` without comments
    """
    out = []
    for line in text.splitlines():
        before, comment = split_comment(line)
        if comment and not before.strip():
            continue

        out.append(before)

    return "\n".join(out)


def parse_aux(aux_file: Path) -> tuple[dict[str, str], list[tuple[str, str]]]:
    """
    Parse the numbers latex gave to labels and headings

    Parameters
    ----------
    aux_file
        `.aux` file written by latex

    Returns
    -------
    :
        The number of each label
        and the level and numbered title of each heading, in order
        (e.g. `("subsection", "3.1 Motivation")`
        or `("section", "Appendix C: Emissions")`)
    """
    labels = {}
    headings = []
    for line in aux_file.read_text().splitlines():
        if match := NEWLABEL.match(line):
            key, number = match.groups()
            if not key.endswith("@cref"):
                labels[key] = number

        elif line.startswith(TOC_ENTRY_START):
            level_end = line.index("}", len(TOC_ENTRY_START))
            level = line[len(TOC_ENTRY_START) : level_end]
            title_start = line.index("{", level_end)
            title = line[title_start + 1 : find_group_end(line, title_start) - 1]
            if title.startswith("\\numberline {"):
                number_start = len("\\numberline ")
                number_end = find_group_end(title, number_start)
                number = title[number_start + 1 : number_end - 1]
                title = f"{number} {title[number_end:]}"

            headings.append((level, title))

    return labels, headings


def replace_brace_macro(text: str, macro: str, prefix: str = "") -> str:
    r"""
    Replace a macro which takes one argument with its argument

    For example, with `prefix="\section*{Competing interests}"`,
    `\competinginterests{None}`
    becomes `\section*{Competing interests}` followed by `None`.

    Parameters
    ----------
    text
        Latex in which to replace

    macro
        Name of the macro (without the backslash)

    prefix
        Text to put in front of the argument

    Returns
    -------
    :
        `text` with the macro replaced
    """
    pattern = re.compile(rf"\\{macro}(?![a-zA-Z])\s*\{{")
    while match := pattern.search(text):
        group_start = match.end() - 1
        group_end = find_group_end(text, group_start)
        body = text[group_start + 1 : group_end - 1]
        text = f"{text[: match.start()]}{prefix}\n{body}\n{text[group_end:]}"

    return text


def rewrite_copernicus_macros(text: str) -> str:
    """
    Rewrite the Copernicus macros which pandoc doesn't understand

    Parameters
    ----------
    text
        Latex to rewrite

    Returns
    -------
    :
        `text`, with Copernicus' macros replaced by standard latex
    """
    for macro, default_heading in (
        ("introduction", "Introduction"),
        ("conclusions", "Conclusions"),
    ):
        text = re.sub(
            rf"\\{macro}(?![a-zA-Z])(\[(?P<heading>[^\]]*)\])?",
            lambda m, h=default_heading: f"\\section{{{m.group('heading') or h}}}",
            text,
        )

    for macro, heading in COPERNICUS_SECTION_MACROS.items():
        text = replace_brace_macro(text, macro, prefix=f"\\section*{{{heading}}}")

    # Table footnotes
    text = replace_brace_macro(text, "belowtable")
    # Table rules
    text = re.sub(r"\\(top|middle|bottom)hline(?![a-zA-Z])", r"\\hline", text)
    # pandoc drops the captions of starred floats
    text = re.sub(r"\\(begin|end)\{(figure|table)\*\}", r"\\\1{\2}", text)

    text = text.replace(
        r"\begin{acknowledgements}", r"\section*{Acknowledgements}"
    ).replace(r"\end{acknowledgements}", "")

    for command in DROPPED_COMMANDS:
        text = re.sub(rf"\\{command}(?![a-zA-Z])", "", text)

    # We pass the bibliography to pandoc ourselves
    text = re.sub(r"\\bibliography(style)?\{[^}]*\}", "", text)
    # natbib configuration
    text = re.sub(r"\\def\\NAT@\w+\{[^}]*\}", "", text)

    return text


def number_headings(text: str, headings: list[tuple[str, str]]) -> str:
    """
    Give each (unstarred) heading its numbered title from latex

    Parameters
    ----------
    text
        Latex in which to number the headings

    headings
        Level and numbered title of each heading, in order (from [parse_aux][])

    Returns
    -------
    :
        `text` with numbered headings

    Raises
    ------
    ValueError
        The headings in `text` don't match `headings`
        (e.g. because the `.aux` file is out of date)
    """
    out = []
    pos = 0
    n_found = 0
    for match in SECTIONING_COMMAND.finditer(text):
        if match.start() < pos:
            # Inside a heading we already replaced
            continue

        level, starred = match.groups()
        if starred:
            continue

        group_start = match.end()
        while text[group_start].isspace():
            group_start += 1

        if text[group_start] == "[":
            # Short title for the table of contents
            group_start = text.index("]", group_start) + 1

        group_end = find_group_end(text, group_start)

        if n_found >= len(headings):
            msg = (
                f"More headings in the latex than in the .aux file ({len(headings)}). "
                "Is the .aux file out of date?"
            )
            raise ValueError(msg)

        aux_level, title = headings[n_found]
        if aux_level != level:
            msg = (
                f"Heading {n_found + 1} is a {level} in the latex "
                f"but a {aux_level} in the .aux file ({title!r}). "
                "Is the .aux file out of date?"
            )
            raise ValueError(msg)

        out.append(text[pos : match.start()])
        out.append(f"\\{level}{{{title}}}")
        pos = group_end
        n_found += 1

    out.append(text[pos:])

    if n_found != len(headings):
        msg = (
            f"Found {n_found} headings in the latex, "
            f"but {len(headings)} in the .aux file. "
            "Is the .aux file out of date?"
        )
        raise ValueError(msg)

    return "".join(out)


def replace_refs(text: str, labels: dict[str, str]) -> str:
    """
    Replace references with the numbers latex gave them

    Parameters
    ----------
    text
        Latex in which to replace references

    labels
        Number of each label (from [parse_aux][])

    Returns
    -------
    :
        `text` with references replaced.
        As in latex, references to unknown labels become "??"
        (with a warning).
    """

    def get_number(key: str) -> str:
        try:
            return labels[key]
        except KeyError:
            print(
                f"WARNING No number for label {key!r} in the .aux file, using '??'",
                file=sys.stderr,
            )
            return "??"

    text = re.sub(r"\\eqref\{([^}]*)\}", lambda m: f"({get_number(m.group(1))})", text)
    text = re.sub(r"\\ref\{([^}]*)\}", lambda m: get_number(m.group(1)), text)

    return text


def number_captions(text: str, labels: dict[str, str]) -> str:
    """
    Put "Figure N: " or "Table N: " in front of each caption

    Floats without a label are left as they are.

    Parameters
    ----------
    text
        Latex in which to number the captions

    labels
        Number of each label (from [parse_aux][])

    Returns
    -------
    :
        `text` with numbered captions
    """

    def number_caption(match: re.Match[str]) -> str:
        env = match.group(0)
        label = re.search(r"\\label\{([^}]*)\}", env)
        caption = re.search(r"\\caption(\[[^\]]*\])?\{", env)
        if label is None or caption is None or label.group(1) not in labels:
            return env

        prefix = f"{match.group(1).capitalize()} {labels[label.group(1)]}: "
        return f"{env[: caption.end()]}{prefix}{env[caption.end() :]}"

    return FLOAT_ENVIRONMENT.sub(number_caption, text)


def convert_figures(text: str, latex_dir: Path, out_dir: Path, dpi: int) -> str:
    """
    Convert PDF figures to PNG, which Word and Google Docs can show

    Parameters
    ----------
    text
        Latex in which to convert the figures

    latex_dir
        Directory relative to which the figures' paths are given

    out_dir
        Directory in which to write the PNGs

    dpi
        Resolution of the PNGs

    Returns
    -------
    :
        `text`, with each PDF figure replaced by its PNG
    """
    out_dir.mkdir(parents=True, exist_ok=True)

    def convert(match: re.Match[str]) -> str:
        path = latex_dir / match.group("path")
        if path.suffix.lower() != ".pdf":
            return match.group(0)

        png = out_dir / f"{path.stem}.png"
        pdf = pypdfium2.PdfDocument(path)
        try:
            pdf[0].render(scale=dpi / 72).to_pil().save(png)
        finally:
            pdf.close()

        return f"\\includegraphics{match.group('options') or ''}{{{png}}}"

    return INCLUDEGRAPHICS.sub(convert, text)


def prepare_latex_for_pandoc(
    latex: str, aux_file: Path, latex_dir: Path, figures_dir: Path, dpi: int
) -> str:
    """
    Rewrite latex so pandoc can turn it into a reviewable Word document

    Parameters
    ----------
    latex
        Latex of the manuscript

    aux_file
        `.aux` file written by latex when compiling `latex`

    latex_dir
        Directory relative to which the figures' paths are given

    figures_dir
        Directory in which to write the converted figures

    dpi
        Resolution of the converted figures

    Returns
    -------
    :
        Latex ready for pandoc
    """
    labels, headings = parse_aux(aux_file)

    res = strip_comments(latex)
    res = rewrite_copernicus_macros(res)
    res = number_headings(res, headings)
    res = number_captions(res, labels)
    res = replace_refs(res, labels)
    # After the headings, as the titles from the .aux file contain these too
    res = re.sub(r"\\protect(?![a-zA-Z])\s*", "", res)
    res = re.sub(r"\\nobreakdash\s*-", "-", res)
    res = convert_figures(res, latex_dir=latex_dir, out_dir=figures_dir, dpi=dpi)
    # citeproc puts the references at the end, give them a heading
    res = res.replace(r"\end{document}", "\\section*{References}\n\\end{document}")

    return res


def export_to_docx(  # noqa: PLR0913
    latex_main: Path,
    output: Path,
    work_dir: Path,
    aux_file: Path | None = None,
    bibliography: Path | None = None,
    dpi: int = 200,
) -> None:
    """
    Export a compiled latex manuscript to Word

    Parameters
    ----------
    latex_main
        Main latex file of the manuscript

    output
        Path in which to write the Word document

    work_dir
        Directory in which to do the conversion

    aux_file
        `.aux` file written by latex when compiling `latex_main`.
        If not supplied, `latex_main` with the suffix `.aux`.

    bibliography
        Bibtex file of the manuscript's references.
        If not supplied, `references.bib` next to `latex_main`.

    dpi
        Resolution of the figures in the Word document
    """
    pandoc = shutil.which("pandoc")
    if pandoc is None:
        msg = (
            "pandoc is needed to export to Word, see https://pandoc.org/installing.html"
        )
        raise FileNotFoundError(msg)

    latex_dir = latex_main.parent
    if aux_file is None:
        aux_file = latex_main.with_suffix(".aux")

    if bibliography is None:
        bibliography = latex_dir / "references.bib"

    if not aux_file.exists():
        msg = (
            f"{aux_file} doesn't exist. "
            "Compile the PDF first, the export uses latex's numbering."
        )
        raise FileNotFoundError(msg)

    if aux_file.stat().st_mtime < latex_main.stat().st_mtime:
        print(
            f"WARNING {aux_file} is older than {latex_main}, "
            "numbers may be out of date. Compile the PDF again to be sure.",
            file=sys.stderr,
        )

    work_dir.mkdir(parents=True, exist_ok=True)
    pandoc_input = work_dir / "main-for-docx.tex"
    pandoc_input.write_text(
        prepare_latex_for_pandoc(
            latex_main.read_text(),
            aux_file=aux_file,
            latex_dir=latex_dir,
            figures_dir=work_dir / "figures",
            dpi=dpi,
        )
    )

    output.parent.mkdir(parents=True, exist_ok=True)
    res = subprocess.run(  # noqa: S603 # inputs are ours
        [
            pandoc,
            str(pandoc_input),
            "--from=latex",
            "--to=docx",
            "--citeproc",
            f"--bibliography={bibliography}",
            f"--resource-path={work_dir}:{latex_dir}",
            f"--output={output}",
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    if res.stderr:
        print(res.stderr, file=sys.stderr)

    if res.returncode != 0:
        msg = f"pandoc failed (exit code {res.returncode}), see its output above"
        raise RuntimeError(msg)
