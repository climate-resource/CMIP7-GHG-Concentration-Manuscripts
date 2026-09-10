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


Redo the figures

I want to make a new version of the current ch4_methods and n2o_methods figure.
I want you to do this from scratch, but you are welcome to copy paste the existing code as needed.
The hard part about these figures seems to be that there are a few constraints:

1. the flying carpet has to be square
1. the maps have a fixed aspect ratio
1. the colour bars are a pain

All the other panels are flexible in terms of their aspect ratio I believe.
Hopefully this gives us the flexibility to basically place the fixed panels first, then fit everything else in around them.

There's a couple of things which I think would work best, but am ultimately flexible about

1. the observational network values panel goes in the top left
1. the figure roughly goes in order of the data processing steps
1. the flying carpet is in the bottom row (bottom right might be based to handle aspect ratio stuff, bottom left is also fine)
1. putting all the maps on the same row will probably work best, so we only pay their fixed aspect ratio price once
1. you might have to make the figure very big in order to get the number of columns you need and essentially trick matplotlib's internal margin book keeping. That is fine, we will let latex scale the figure back down in the paper, so the actual size doesn't matter that much, just the ratios between things
