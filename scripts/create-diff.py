"""
Create diff between versions
"""

import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Annotated

import typer
from git import Repo


def run(cmd, cwd=None):
    """Run a shell command with error checking."""
    print(f"INFO Running: {cmd}")
    result = subprocess.run(cmd, shell=True, cwd=cwd, check=False)
    if result.returncode != 0:
        sys.exit(f"ERROR Command failed: {cmd}")


def main(
    notebook_stem: Annotated[
        str,
        typer.Argument(
            help="""Stem of the notebook for which to create the diff.

e.g. `user-guide-historical`.
The file `notebooks/{notebook_stem}.py` must exist."""
        ),
    ],
    out_title: Annotated[
        str,
        typer.Option(help="Title to use when creating the PDFs"),
    ],
    compare_against: Annotated[
        str,
        typer.Option(
            help="Git object to compare against. Can be a commit ID or branch name."
        ),
    ],
    diff_dir: Annotated[
        Path,
        typer.Option(help="Directory in which to store diffing components"),
    ] = Path("book-diff"),
):
    """Create the diff PDF"""
    if not diff_dir.exists():
        diff_dir.mkdir(parents=True)

    REPO_ROOT = Path(__file__).parents[1]
    repo = Repo(REPO_ROOT)
    # Before doing anything more, check we can go to our `compare_commit` and back
    # i.e. there aren't any untracked clashes we need to commit.
    current_branch = repo.active_branch

    # repo.git.checkout(compare_against)
    # print(f"Checked out {compare_against}")
    # shutil.rmtree("book/_build/latex")
    # run(f'uv run python scripts/create_latex_pdf.py {notebook_stem} "{out_title}"')
    # run("git checkout -- book")
    # if (diff_dir / compare_against).exists():
    #     shutil.rmtree(diff_dir / compare_against)
    # shutil.copytree("book/_build/latex", diff_dir / compare_against)
    #
    # repo.git.checkout(current_branch)
    # print("Checked out the current branch")
    # shutil.rmtree("book/_build/latex")
    # run(f'uv run python scripts/create_latex_pdf.py {notebook_stem} "{out_title}"')
    # run("git checkout -- book")
    # if (diff_dir / current_branch.name).exists():
    #     shutil.rmtree(diff_dir / current_branch.name)
    # shutil.copytree("book/_build/latex", diff_dir / current_branch.name)

    pwd = os.getcwd()
    os.chdir(diff_dir)

    run(
        " ".join(
            [
                # unix-specific?
                "latexdiff",
                f"{compare_against}/projectnamenotset.tex",
                f"{current_branch.name}/projectnamenotset.tex",
                "--flatten",
                # unix-specific?
                "> diff.tex",
            ]
        )
    )

    if (Path(compare_against) / "projectnamenotset.bbl").exists():
        run(
            " ".join(
                [
                    # unix-specific?
                    "latexdiff",
                    f"{compare_against}/projectnamenotset.bbl",
                    f"{current_branch.name}/projectnamenotset.bbl",
                    # unix-specific?
                    "> diff.bbl",
                ]
            )
        )
    else:
        print(
            f"WARNING: no .bbl for {compare_against}, "
            f"just using the one from {current_branch.name}"
        )
        shutil.copyfile(f"{current_branch.name}/projectnamenotset.bbl", "diff.bbl")

    for file in Path(compare_against).glob("*.sty"):
        shutil.copy(file, ".")

    for file in Path(current_branch.name).glob("*.sty"):
        shutil.copy(file, ".")

    for file in Path(compare_against).glob("*.cls"):
        shutil.copy(file, ".")

    for file in Path(current_branch.name).glob("*.cls"):
        shutil.copy(file, ".")

    for file in Path(compare_against).glob("*.png"):
        shutil.copy(file, ".")

    for file in Path(current_branch.name).glob("*.png"):
        shutil.copy(file, ".")

    subprocess.run("xelatex -interaction=nonstopmode diff.tex", shell=True, check=False)
    subprocess.run("xelatex -interaction=nonstopmode diff.tex", shell=True, check=False)

    out_file = f"{notebook_stem}-diff_{current_branch.name}_{compare_against}.pdf"
    shutil.copy("diff.pdf", out_file)
    print(f"INFO: diff available in {diff_dir / out_file}")

    os.chdir(pwd)


if __name__ == "__main__":
    typer.run(main)
