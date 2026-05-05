"""
Write jupyter-book table of contents (TOC) file for a specific build

This assumes that the book is being built based on a single file.

We build multiple outputs here,
so we want to be more custom than is supported by having a single TOC
and we want to avoid trying to keep multiple configs in sync.
"""

from pathlib import Path
from typing import Annotated

import typer
import yaml

REPO_ROOT = Path(__file__).parents[1]


def main(
    output: Annotated[
        Path,
        typer.Option(
            help="Path in which to write the output TOC file",
            dir_okay=False,
            file_okay=True,
        ),
    ],
    source_file: Annotated[
        Path, typer.Option(help="Path to the source file we will compile")
    ],
) -> None:
    # Also need this intro file for some reason...
    with open(source_file.parent / "intro.md", "w") as fh:
        fh.write("# Introduction")

    res = {
        "format": "jb-book",
        "root": "intro",
        "options": {"maxdepth": 3},
        "chapters": [{"file": source_file.stem}],
    }

    output.parent.mkdir(exist_ok=True, parents=True)
    output.write_text(yaml.dump(res, sort_keys=True), encoding="utf-8")


if __name__ == "__main__":
    typer.run(main)
