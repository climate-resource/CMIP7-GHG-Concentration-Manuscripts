"""
Helper function to convert jupytext .py into .md files
"""

#!/usr/bin/env python3
import subprocess
from pathlib import Path

# Source and target directories
source_dir = Path("notebooks")
output_dir = Path("book")

# Make sure output directory exists
output_dir.mkdir(parents=True, exist_ok=True)

# Loop over all .py files in the source directory
for py_file in source_dir.glob("*.py"):
    base = py_file.stem  # filename without extension
    output_file = output_dir / f"{base}.md"

    # Run the jupytext command
    subprocess.run(
        [
            "uv",
            "run",
            "jupytext",
            "--to",
            "myst",
            str(py_file),
            "--output",
            str(output_file),
        ],
        check=True,
    )

    print(f"Converted {py_file} -> {output_file}")
