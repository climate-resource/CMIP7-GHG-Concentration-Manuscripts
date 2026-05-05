"""
Write jupyter-book config for a specific build

We build multiple outputs here,
so we want to be more custom than is supported by having a single config
and we want to avoid trying to keep multiple configs in sync.
"""

import shutil
from pathlib import Path
from typing import Annotated

import typer
import yaml

REPO_ROOT = Path(__file__).parents[1]


def main(
    base: Annotated[
        Path,
        typer.Option(
            help="Path to the base config",
            dir_okay=False,
            file_okay=True,
        ),
    ],
    output: Annotated[
        Path,
        typer.Option(
            help="Path in which to write the output config",
            dir_okay=False,
            file_okay=True,
        ),
    ],
    title: Annotated[str, typer.Option(help="Title of the output")],
    description: Annotated[str, typer.Option(help="Description of the output")],
    references_bib_file: Annotated[
        Path, typer.Option(help="Bibtex references file to use")
    ],
    source_file: Annotated[
        Path, typer.Option(help="Path to the source file we will compile")
    ],
    source_file_raw: Annotated[
        Path,
        typer.Option(help="Path to the raw source file (which exists in the git repo)"),
    ],
    references_bib_file_out: Annotated[
        Path | None,
        typer.Option(
            help=(
                "Where to write the bibtex file. "
                "If not supplied, we write it in the same directory as `--source-file`"
            )
        ),
    ] = None,
) -> None:
    if references_bib_file_out is None:
        references_bib_file_out = output.parent / references_bib_file.name

    base_config = yaml.safe_load(base.read_text(encoding="utf-8"))
    latex_targetname = f"{source_file.stem}.tex"

    res = {
        **base_config,
        "title": title,
        "description": description,
        "bibtex_bibfiles": [
            str(
                references_bib_file_out.absolute().relative_to(
                    source_file.parent.absolute()
                )
            )
        ],
    }
    res["repository"]["path_to_book"] = str(
        source_file_raw.absolute().relative_to(REPO_ROOT)
    )
    # Avoid Sphinx's default fallback of "Project name not set", which
    # generates projectnamenotset.tex for LaTeX builds.
    res.setdefault("sphinx", {}).setdefault("config", {})["project"] = title
    res.setdefault("latex", {}).setdefault("latex_documents", {})["targetname"] = (
        latex_targetname
    )
    output.parent.mkdir(exist_ok=True, parents=True)
    output.write_text(yaml.dump(res, sort_keys=True), encoding="utf-8")

    shutil.copy2(references_bib_file, references_bib_file_out)


if __name__ == "__main__":
    typer.run(main)
