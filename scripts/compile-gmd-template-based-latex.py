"""
Compile a GMD-templated based latex document to PDF

The work is done by [local.gmd_latex][], this is just its command-line interface.
"""

import json
import shutil
import sys
from pathlib import Path
from typing import Annotated

import typer

from local.cmip_ghg_generation import DEFAULT_BUNDLE_DIR
from local.gmd_latex import assemble_main_tex, compile_latex
from local.historical_ghg_forcing_for_cmip7.value_checks import get_value_checks
from local.value_checks import run_value_checks


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


def load_tex_inputs_manifest(
    manifest_file: Path,
) -> tuple[list[tuple[str, Path]], list[tuple[str, Path]]]:
    """
    Load a manifest of figure files and files to inline

    Parameters
    ----------
    manifest_file
        JSON file with a `figures` and an `inline` entry,
        each of which maps a tag to replace to a full path

    Returns
    -------
    :
        The figures and the files to inline, as (tag, path) pairs
    """
    manifest = json.loads(manifest_file.read_text())

    figure_files = [(tag, Path(path)) for tag, path in manifest["figures"].items()]
    inline_files = [(tag, Path(path)) for tag, path in manifest["inline"].items()]

    return figure_files, inline_files


def main(  # noqa: PLR0913
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
                "e.g. `--figure-file=@n2o-methods-figure@=/path/to/figure.pdf`"
            )
        ),
    ] = None,
    inline_file: Annotated[
        list[str] | None,
        typer.Option(
            help=(
                "File whose text to put in place of a tag, "
                "e.g. a generated table. "
                "Should be passed as `tag-to-replace=path-to-file`, "
                "e.g. `--inline-file=@per-gas-table@=/path/to/table.tex`. "
                "Tags in comments are left as they are."
            )
        ),
    ] = None,
    tex_inputs_manifest: Annotated[
        Path | None,
        typer.Option(
            help=(
                "JSON file of figure files and files to inline. "
                "It should have a `figures` and an `inline` entry, "
                "each of which maps a tag to replace in the latex "
                "to the path of the file, "
                "i.e. the same information as "
                "`--figure-file` and `--inline-file`. "
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
    bundle_dir: Annotated[
        Path,
        typer.Option(
            help=(
                "Directory which holds the original run's bundle. "
                "Only used to check values (see `--check-values`)."
            ),
            dir_okay=True,
            file_okay=False,
        ),
    ] = DEFAULT_BUNDLE_DIR,
) -> None:
    """
    Compile the PDF
    """
    # The manifest gives (tag, path) pairs straight away,
    # only the CLI values need splitting
    if tex_inputs_manifest is not None:
        manifest_figure_files, manifest_inline_files = load_tex_inputs_manifest(
            tex_inputs_manifest
        )
    else:
        manifest_figure_files, manifest_inline_files = [], []

    cli_figure_files = [
        split_tag_and_path(v, option="--figure-file") for v in figure_file or []
    ]
    cli_inline_files = [
        split_tag_and_path(v, option="--inline-file") for v in inline_file or []
    ]

    latex_main, collector = assemble_main_tex(
        abstract=abstract,
        introduction=introduction,
        sections=section,
        conclusion=conclusion,
        code_and_data_availability=code_and_data_availability,
        author_contribution=author_contribution,
        competing_interests=competing_interests,
        acknowledgements=acknowledgements,
        replacements=replacements,
        metadata=metadata,
        references_bib_file=references_bib_file,
        copernicus_template_dir=copernicus_template_dir,
        clean_copernicus_template_filename=clean_copernicus_template_filename,
        build_dir=build_dir,
        appendices=appendix or (),
        extra=extra or (),
        auxiliary=auxiliary or (),
        figure_files=[*cli_figure_files, *manifest_figure_files],
        files_to_inline=cli_inline_files,
        files_to_inline_if_tagged=manifest_inline_files,
    )

    # Last thing before building, so every file has been read
    # (and its value check comments collected)
    if check_values:
        value_check_report = run_value_checks(
            collector.specs, get_value_checks(bundle_dir=bundle_dir)
        )
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
