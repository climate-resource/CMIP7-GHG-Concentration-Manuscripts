"""
Export a compiled GMD-template based latex manuscript to Word

The work is done by [local.docx_export][], this is just its command-line interface.
The resulting `.docx` can be uploaded to Google Drive and opened with Google Docs.
"""

from pathlib import Path
from typing import Annotated

import typer

from local.docx_export import export_to_docx


def main(  # noqa: PLR0913
    latex_main: Annotated[
        Path,
        typer.Option(
            help=(
                "Main latex file of the manuscript. "
                "It must have been compiled already "
                "(the export uses the numbers in the `.aux` file)."
            ),
            dir_okay=False,
            file_okay=True,
            exists=True,
        ),
    ],
    output: Annotated[
        Path,
        typer.Option(
            help="Path in which to write the Word document",
            dir_okay=False,
            file_okay=True,
        ),
    ],
    work_dir: Annotated[
        Path,
        typer.Option(
            help="Directory in which to do the conversion",
            dir_okay=True,
            file_okay=False,
        ),
    ],
    aux_file: Annotated[
        Path | None,
        typer.Option(
            help=(
                "`.aux` file written when compiling `--latex-main`. "
                "Defaults to `--latex-main` with the suffix `.aux`."
            ),
            dir_okay=False,
            file_okay=True,
        ),
    ] = None,
    bibliography: Annotated[
        Path | None,
        typer.Option(
            help=(
                "Bibtex file of the manuscript's references. "
                "Defaults to `references.bib` next to `--latex-main`."
            ),
            dir_okay=False,
            file_okay=True,
        ),
    ] = None,
    dpi: Annotated[int, typer.Option(help="Resolution of the figures")] = 200,
) -> None:
    """
    Export to Word
    """
    export_to_docx(
        latex_main=latex_main,
        output=output,
        work_dir=work_dir,
        aux_file=aux_file,
        bibliography=bibliography,
        dpi=dpi,
    )


if __name__ == "__main__":
    app = typer.Typer(pretty_exceptions_show_locals=False)
    app.command()(main)
    app()
