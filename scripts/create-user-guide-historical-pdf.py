"""
Create the PDF version of the historical user guide
"""

from pathlib import Path

import jupytext
import papermill
from loguru import logger
from nbconvert.exporters import PDFExporter
from nbconvert.preprocessors import TagRemovePreprocessor


def main():
    force_rerun = False
    force_rerun = True

    REPO_ROOT = Path(__file__).parents[1]

    IN_FILE = REPO_ROOT / "notebooks" / "user-guide-historical.py"
    OUT_FILE = REPO_ROOT / "output-pdfs" / "user-guide-historical.pdf"
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
