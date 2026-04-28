"""
Compile jupyter-book latex to PDF
"""

import copy
import shutil
import subprocess
from pathlib import Path
from typing import Annotated

import typer

REPO_ROOT = Path(__file__).parents[1]


def main(
    source: Annotated[
        Path,
        typer.Option(
            help="Path to the source file",
            dir_okay=False,
            file_okay=True,
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
    references_bib_file: Annotated[
        Path, typer.Option(help="Bibtex references file to use")
    ],
) -> None:
    with open(source, encoding="utf-8") as fh:
        tex_jupyter_book_l = [v.strip() for v in fh.readlines()]

    tex_mod_l = copy.deepcopy(tex_jupyter_book_l)
    # Assumes use of biber in preamble
    for line in [
        "",
        "% Bibliography managed by biber.",
        "% Command injected in `create_latex_pdf.py`",
        r"\printbibliography[heading=bibintoc,title={References}]",
        "",
    ]:
        tex_mod_l.insert(-1, line)

    tex_mod = "\n".join(tex_mod_l)

    source_tmp_stem = f"{source.stem}-prepped"
    source_tmp = source.parent / f"{source_tmp_stem}.tex"
    with open(source_tmp, "w", encoding="utf-8") as fh:
        fh.write(tex_mod)

    subprocess.run(
        f"xelatex {source_tmp.name}", shell=True, cwd=source_tmp.parent, check=True
    )

    shutil.copy(references_bib_file, source_tmp.parent / references_bib_file.name)
    subprocess.run(
        f"biber {source_tmp.name.replace('.tex', '.bcf')}",
        shell=True,
        cwd=source_tmp.parent,
        check=True,
    )

    subprocess.run(
        f"xelatex {source_tmp.name}", shell=True, cwd=source_tmp.parent, check=True
    )
    subprocess.run(
        f"xelatex {source_tmp.name}", shell=True, cwd=source_tmp.parent, check=True
    )

    built_pdf = source_tmp.parent / f"{source_tmp_stem}.pdf"
    if not built_pdf.exists():
        raise FileNotFoundError(built_pdf)

    output.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(built_pdf, output)


if __name__ == "__main__":
    typer.run(main)
