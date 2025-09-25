"""
Create the PDF version of the historical user guide
"""

import datetime as dt
import os.path
from pathlib import Path

import jupytext
import papermill
from nbconvert.exporters import PDFExporter
from nbconvert.preprocessors import TagRemovePreprocessor


def main():
    force_rerun = False
    force_rerun = True

    output_title = "Historical dataset - data description and user guide"

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
        jupytext.write(notebook_jupytext, unexecuted_notebook, fmt="ipynb")

        if executed_notebook.exists():
            executed_notebook.unlink()

        papermill.execute_notebook(unexecuted_notebook, executed_notebook)

    exporter = PDFExporter()
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

    modified_date = dt.datetime.fromtimestamp(
        os.path.getmtime(executed_notebook), tz=dt.timezone.utc
    ).strftime("%B %-d, %Y")
    with open(executed_notebook) as fh:
        output = exporter.from_file(
            fh,
            resources={
                "metadata": {
                    "name": output_title,
                    "path": executed_notebook,
                    "modified_date": modified_date,
                }
            },
        )

    # Write to output html file
    OUT_FILE.parent.mkdir(exist_ok=True, parents=True)
    with open(OUT_FILE, "wb") as f:
        f.write(output[0])


if __name__ == "__main__":
    main()
