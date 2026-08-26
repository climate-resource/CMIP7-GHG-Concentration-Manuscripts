"""
Generation of the N2O methods figure
"""

from __future__ import annotations

from pathlib import Path


def generate_n2o_methods_figure(outfile: Path) -> Path:
    """
    Generate the N2O methods figure
    """
    outfile.parent.mkdir(exist_ok=True, parents=True)
    with open(outfile, "w") as fh:
        fh.write("hi")
