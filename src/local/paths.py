"""
Path helpers
"""

from pathlib import Path

REPO_ROOT = Path(__file__).parents[2]

DATA_RAW_DIR = REPO_ROOT / "data" / "raw"
"""Directory in which we keep raw data

Almost everything in here is downloaded or copied in,
hence almost all of it is ignored by git.
The exceptions are explicitly un-ignored in `.gitignore`.
"""

FIGURES_DIR = REPO_ROOT / "figures"
"""Directory in which we write generated figures"""
