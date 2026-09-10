- PRIMAP regression

PRIMAP_REGRESSION_DATA_FILE = Path("manuscript-outputs") / "ch4_primap-regression-data.nc"
"""Where the re-run notebook saves the PRIMAP regression data we want

Relative to the bundle's root directory,
because that is the notebook's working directory.
"""

PRIMAP_REGRESSION_YEARS_FILE = Path("manuscript-outputs") / "ch4_primap-regression-years.json"
"""Where the re-run notebook saves the PRIMAP regression years information we want

Relative to the bundle's root directory,
because that is the notebook's working directory.
"""

    save_cell_primap_data = f"""
# Added for the CMIP7 GHG manuscript.
# The original run never saved these pieces out,
# but we need them for our plotting
from pathlib import Path

years_to_fill_with_regression
primap_regression_data_file = Path("{PRIMAP_REGRESSION_DATA_FILE.as_posix()}")
primap_regression_data_file.parent.mkdir(exist_ok=True, parents=True)
primap_regression_data.to_netcdf(primap_regression_data_file)
primap_regression_data_file
"""

    save_cell_primap_years = f"""
years_to_fill_with_regression
years_to_fill_with_regression_file = Path("{PRIMAP_REGRESSION_YEARS_FILE.as_posix()}")
years_to_fill_with_regression_file.parent.mkdir(exist_ok=True, parents=True)
with open(years_to_fill_with_regression_file, "w") as fh:
    json.dump([int(v) for v in years_to_fill_with_regression], fh)

years_to_fill_with_regression_file
"""

- need to colour/dash/marker the bits of the lat. grad. PC extension
    - N2O: obs. network vs constant extension
    - CH4: obs. network vs. regression driven vs. ice core optimised vs. constant extension

-
