"""
Create the PDF version of the scenario user guide

Copy-paste number 2: when we get to 3, we should probably stop
copy-pasting and extract the common patterns into something useful,
but we can do that as part of
https://github.com/climate-resource/CMIP7-GHG-Concentration-Manuscripts/issues/4
"""

import re
from pathlib import Path

import jupytext
import nbconvert
import nbformat
import papermill
from loguru import logger
from nbconvert.exporters import PDFExporter
from nbconvert.preprocessors import TagRemovePreprocessor
from nbconvert.preprocessors.base import Preprocessor


class LatexSubSuperScriptPreprocessor(Preprocessor):
    """
    Pre-processor which converts markdown sub and superscript tags into latex equivalent
    """

    def preprocess_cell(
        self,
        cell: nbformat.NotebookNode,
        resources: nbconvert.exporters.ResourcesDict,
        cell_index: int,
    ) -> tuple[nbformat.NotebookNode, nbconvert.exporters.ResourcesDict]:
        """
        Pre-process a cell before passing it to the exporter

        Parameters
        ----------
        cell
            Cell to pre-process

        resources
            Resources available to help with pre-processing

        cell_index
            Cell index

        Returns
        -------
        :
            Processed cell and resources
        """
        if cell.cell_type == "markdown" and (
            "<sub>" in cell.source or "<sup>" in cell.source
        ):
            cell.source = re.sub(
                r"<sub>(.*?)</sub>", r"\\textsubscript{\1}", cell.source
            )
            cell.source = re.sub(
                r"<sup>(.*?)</sup>", r"\\textsuperscript{\1}", cell.source
            )

        return cell, resources


def main():
    """
    Create the PDF of the user guide
    """
    force_rerun = False
    force_rerun = True

    REPO_ROOT = Path(__file__).parents[1]

    IN_FILE = REPO_ROOT / "notebooks" / "user-guide-scenarios.py"
    OUT_FILE = REPO_ROOT / "output-pdfs" / "user-guide-scenarios.pdf"
    UNEXECUTED_NOTEBOOKS_DIR = REPO_ROOT / "notebooks-unexecuted"
    EXECUTED_NOTEBOOKS_DIR = REPO_ROOT / "notebooks-executed"

    executed_notebook = EXECUTED_NOTEBOOKS_DIR / IN_FILE.name.replace(".py", ".ipynb")
    if not executed_notebook.exists() or force_rerun:
        # Execute notebook from .py

        notebook_jupytext = jupytext.read(IN_FILE)
        unexecuted_notebook = UNEXECUTED_NOTEBOOKS_DIR / IN_FILE.name.replace(
            ".py", ".ipynb"
        )
        if unexecuted_notebook.exists():
            unexecuted_notebook.unlink()

        unexecuted_notebook.parent.mkdir(parents=True, exist_ok=True)
        logger.info(
            f"Converting {IN_FILE.relative_to(REPO_ROOT)} "
            f"to {unexecuted_notebook.relative_to(REPO_ROOT)} "
            "using jupytext"
        )
        jupytext.write(notebook_jupytext, unexecuted_notebook, fmt="ipynb")

        if executed_notebook.exists():
            executed_notebook.unlink()

        notebook_parameters = {}
        logger.info(
            f"Executing {unexecuted_notebook.relative_to(REPO_ROOT)} "
            f"to {executed_notebook.relative_to(REPO_ROOT)} "
            f"using papermill with {notebook_parameters=}"
        )
        papermill.execute_notebook(
            unexecuted_notebook, executed_notebook, notebook_parameters
        )

    exporter = PDFExporter(
        # template_file="article",
        template_file="report",
        # verbose=True,
    )
    # Inherits from LatexExporter (https://github.com/jupyter/nbconvert/blob/main/nbconvert/exporters/latex.py#L18)
    # so we could also do funky latex stuff if we wanted
    exporter.register_preprocessor(
        TagRemovePreprocessor(
            remove_cell_tags=("remove_cell",),
            remove_all_outputs_tags=("remove_output",),
            remove_input_tags=("remove_input",),
        ),
        enabled=True,
    )

    exporter.register_preprocessor(
        LatexSubSuperScriptPreprocessor(),
        enabled=True,
    )

    logger.info(
        f"Converting {executed_notebook.relative_to(REPO_ROOT)} to PDF using nbconvert"
    )
    output = exporter.from_filename(executed_notebook)

    OUT_FILE.parent.mkdir(exist_ok=True, parents=True)
    logger.info(f"Writing result to {OUT_FILE.relative_to(REPO_ROOT)}")
    with open(OUT_FILE, "wb") as f:
        f.write(output[0])


if __name__ == "__main__":
    main()
