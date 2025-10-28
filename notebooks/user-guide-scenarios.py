# ---
# jupyter:
#   authors:
#   - name: Zebedee Nicholls
#   - name: Florence Bockting
#   - name: Mika Plf{\"u}ger
#   jupytext:
#     notebook_metadata_filter: title,authors
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.17.3
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
#   title: 'CMIP Greenhouse Gas (GHG) Concentration Scenario Dataset:
#
#     Data Description and User Guide'
# ---

# %% [markdown] editable=true slideshow={"slide_type": ""}
# # Overview
#
# Here we provide a short description of the draft scenario dataset
# and a guide for users.
# This is intended to provide a short introduction for users of the data.
# The full details of the dataset's construction
# and evaluation against other data sources
# will be provided in the full manuscript which is being prepared.

# %% [markdown] editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
# ## Imports

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
import calendar
from functools import partial
from typing import Any

import cftime
import matplotlib
import matplotlib.pyplot as plt
import nc_time_axis  # noqa: F401
import numpy as np
import numpy.typing as npt
import pandas as pd
import pandas_indexing as pix
import pandas_openscm
import seaborn as sns
import tqdm.auto

from local.data_loading import (
    fetch_and_load_ghg_dataset,
    fetch_and_load_ghg_dataset_scenarios,
    get_ghg_dataset_local_files,
)
from local.esgf.db_helpers import create_all_tables, get_sqlite_engine
from local.esgf.search.search_query import KnownIndexNode
from local.paths import REPO_ROOT

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
pandas_openscm.register_pandas_accessors()

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
local_data_root_dir = REPO_ROOT / "data" / "raw" / "esgf"
local_data_root_dir.mkdir(exist_ok=True, parents=True)
sqlite_file = REPO_ROOT / "download-test-scenario-database.db"
# # Obviously we wouldn't delete the database every time
# # in production, but while experimenting it's handy
# # to always start with a clean slate.
# if sqlite_file.exists():
#     sqlite_file.unlink()

engine = get_sqlite_engine(sqlite_file)
create_all_tables(engine)


# %% [markdown] editable=true slideshow={"slide_type": ""}
# # Dataset construction
#
# The dataset is constructed following a similar methodology
# to {cite:t}`meinshausen_shared_2020`.
# The full method will be described in a forthcoming paper.
# In brief, the method is:
#
# 1. **retrieve concentrations** of GHG gases
#
#     1. for gases covered under the Montreal Protocol (TODO REF) and whose
#        concentration evolution is already specified in WMO 2022 (TODO, ref),
#        we simply use the WMO 2022 concentrations
#
#     1. for all other gases, we
#
#         1. start with harmonised, complete emissions sets for each scenario
#            of interest
#
#         1. run MAGICC (TODO REF) to translate these emissions into
#             global-mean concentrations
#
#            1. MAGICC is run in the same configuration which was used in AR6
#              (as described/evaluated in Cross-Chapter Box 7.1[^1]).
#              This represents our best estimate, in line with the last IPCC report,
#              of the concentrations that result from the emissions.
#              Note, given that CMIP7 models have not yet been run, this will,
#              almost by definition, produce different concentrations
#              than those outputed by the ESMs run in emissions-driven mode.
#              We nonetheless look forward to reading lots of papers saying that
#              "CMIP7 input concentrations are biased up/down/whatever".
#         1. harmonise the global-mean concentrations to the historical
#              concentrations using gradient-aware harmonisation[^2] to ensure a smooth
#              transition in both the absolute values and the gradient, improving on the
#              CMIP6 data which had an abrupt jump in the gradient (most notable for
#              methane, further details below)
#
# 3. **seasonality and latitudinal gradients**: We use the same statistical
#    models for seasonality and latitudinal gradients as were used over the
#    historical period. These are then applied to the future, based on the same
#    input drivers.
#     - one exception is CO<sub>2</sub>, which uses GPP as the input driver
#       rather than a regression on a combined global-mean surface temperature and
#       CO<sub>2</sub> concentration metric (as is used in the product that covers
#       history)
#     - arguably we should harmonise the inputs to these models too.
#       We don't do this step as such a harmonisation would only provide
#       second-order corrections.
# 5. **combine the global-means, seasonality and latitudinal gradients** to
#    produce our complete, gridded projections
# 6. **calculate lower resolution products** from the gridded product
# 7. **release** all the data
#
# [^1]: https://www.ipcc.ch/report/ar6/wg1/chapter/chapter-7/
# [^2]: https://github.com/climate-resource/gradient-aware-harmonisation

# %% [markdown] editable=true slideshow={"slide_type": ""}
# # Finding and accessing the data

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ## ESGF
#
# The Earth System Grid Federation (ESGF, TODO REF) provides
# access to a range of climate data.
#
# The scenario data of interest here,
# which is a draft dataset
# can be found under "MIP era" `CMIP6Plus` (for draft datasets),
# "institution ID" `CR`
# and "source version" `0.1.0`
# (also under "source IDs" of the form `CR-*-0-1-0`,
# although this is less useful as the ESGF search API
# does not appear to support glob/star expressions
# for specific facet searches).
#
# It is possible to filter searches on ESGF
# via the user interface.
# Searches can often be encoded in URLs[^3] too
# (although these URLs sometimes move,
# so we make no guarantee that this link will always be live).
#
# These datasets are a draft only.
# The final datasets will maintain the same form; however,
# the numbers are currently being finalised
# and the final names have been changed since the draft datasets
# were published,
# so please take care to treat the values shown here as drafts only.
#
# To download the data, we recommend accessing it directly via the
# ESGF user interfaces via links as discussed above.
# Alternately, there are tools dedicated to accessing ESGF data,
# with two prominent examples being:
#
# 1. esgpull: https://esgf.github.io/esgf-download
# 2. intake-esgf: https://intake-esgf.readthedocs.io
#
# Please see these tools' docs for usage instructions.
#
# [^3]: An example URL: [https://esgf-node.ornl.gov/search?project=input4MIPs&activeFacets=%7B%22source_version%22%3A%<br>
# 220.1.0%22%2C%22institution_id%22%3A%22CR%22%2C%22mip_era%22%3A%22CMIP6Plus%22%7D](
# https://esgf-node.ornl.gov/search?project=input4MIPs&activeFacets=%7B%22source_version%22%3A%220.1.0%22%2C%22institution_id%22%3A%22CR%22%2C%22mip_era%22%3A%22CMIP6Plus%22%7D)

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ## Zenodo
#
# While it aims to be, the ESGF is technically not a permanent archive
# and does not issue DOIs.
# In order to provide more reliable, citable access to the data,
# we will also provide the final scenario datasets on Zenodo
# (TODO REF, although we have not done this step for the draft datasets).
# When ready, we will update this guide to use the final scenario data
# and include the zenodo link to the source code and input data used to process it.

# %% [markdown] editable=true slideshow={"slide_type": ""}
# # Data description

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ## Format
#
# The data is provided in netCDF format [TODO citation].
# This self-describing format allows the data
# to be placed in the same file as metadata
# (in the so-called "file header").
# The scenario datasets are only 78 years long,
# so each scenario-specific dataset is provided in a single file
# (unlike the historical dataset,
# which is split into multiple files).

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ## Scenario information
#
# Determining the scenario to which each dataset applies
# is not trivial
# (there was [discussion](https://github.com/PCMDI/input4MIPs_CVs/discussions/64)
# about how to make this more trivial,
# but ultimately backwards-compatibility was prioritised).
#
# Each dataset is given a unique source ID.
# This source ID appears both in the global attributes of the file
# as well as in the filename (as the 5th underscore separated element).
# The concept of a source ID is a bit of a perculiar one
# to CMIP forcings data.
# In practice, it is simply a unique identifier for a collection of datasets
# (and it's best not to read more than that into it).
#
# As a result of the way that scenario data is handled in CMIP/input4MIPs,
# the scenario name appears as part of the source ID,
# rather than as a standalone attribute/identifier.
# This means that its value must be extracted from the other parts of the source ID.
# In general, this is not a trivial problem
# (users of forcings data more generally are directed to
# [input4mips-cvs.readthedocs.io/en/latest/extracting-scenario-from-source-id](https://input4mips-cvs.readthedocs.io/en/latest/extracting-scenario-from-source-id/)).
# For the greenhouse gas concentrations,
# the scenario identifier is simply
# the second hyphen-separated element of the source ID.
# For example, for the source ID `CR-ml-0-1-0`,
# the scenario identifier is `ml`.
# A Python function for doing this extraction is below.


# %% editable=true slideshow={"slide_type": ""}
def extract_scenario_id(source_id: str) -> str:
    """
    Extract scenario ID from a GHG concentration source ID

    Parameters
    ----------
    source_id
        Source ID from which to extract the scenario ID

    Returns
    -------
    :
        Scenario ID
    """
    return source_id.split("-")[1]


print(f"{extract_scenario_id('CR-ml-0-1-0')=}")
print(f"{extract_scenario_id('CR-l-0-1-0')=}")

# %% [markdown] editable=true slideshow={"slide_type": ""}
# These scenario IDs can then be used to find details of the complete scenario.
# These details will be provided both via the CMIP CVs
# (see https://github.com/WCRP-CMIP/CMIP7-CVs)
# and the final ScenarioMIP paper
# (revisions of
# [this paper](https://doi.org/10.5194/egusphere-2024-3765) are expected soon).
# As above, note that the scenario IDs
# have changed since publication of the draft dataset. Table \ref{tab:scenario_ids}
# provides an overview
# of the changes. Scenario IDs of the final datasets can be confirmed
# [here](https://github.com/WCRP-CMIP/CMIP7-CVs/discussions/1#discussioncomment-14585785).
#
# **Table 4.1:** Mapping of draft scenario IDs to final scenario IDs.
#
# \begin{table}[ht]
# \centering
# \caption{Mapping of draft scenario IDs to final scenario IDs}
# \label{tab:scenario_ids}
# \begin{tabular}{cc}
# \hline
# Draft dataset & Final dataset \\
# \hline
# \textbf{vllo} & \textbf{vl} \\
# \textbf{vlho} & \textbf{ln} \\
# l    & l  \\
# ml   & ml \\
# m    & m  \\
# hl   & hl \\
# h    & h  \\
# \hline
# \end{tabular}
# \end{table}
#

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ## Grids and frequencies provided
#
# We provide five combinations of grids and time sampling
# (also referred to as *frequency*,
# although this is a bit of a misuse as the units of frequency are per time,
# which doesn't match the convention for these metadata values).
# The grid and frequency information for each file can be found in its netCDF header
# under the attributes `grid_label` (for grid) and `frequency` (for time sampling).
# The `grid_label` and `frequency` also appear in each file's name,
# which allows files to be filtered without needing to load them first.
# The five combinations of grid and time sampling are:
#
# 1. global-, annual-mean (`grid_label="gm"`, `frequency="yr"`)
# 1. global-, monthly-mean (`grid_label="gm"`, `frequency="mon"`)
# 1. hemispheric-, annual-mean (`grid_label="gr1z"`, `frequency="yr"`)
# 1. hemispheric-, monthly-mean (`grid_label="gr1z"`, `frequency="mon"`)
# 1. 15-degree latitudinal, monthly-mean (`grid_label="gnz"`, `frequency="mon"`)

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ## Species provided
#
# We provide concentrations for 43 greenhouse gas concentrations and species,
# as well as three equivalent species.
# The species are:
#
# <!-- Note: generated using `scripts/generate-ghg-listing.py` --->
# - major greenhouse gases (3)
#     - CH<sub>4</sub>, CO<sub>2</sub>, N<sub>2</sub>O
# - ozone-depleting substances (17)
#     - CFCs (5)
#         - CFC-11, CFC-113, CFC-114, CFC-115, CFC-12
#     - HCFCs (3)
#         - HCFC-141b, HCFC-142b, HCFC-22
#     - Halons (3)
#         - Halon 1211, Halon 1301, Halon 2402
#     - other ozone-depleting substances (6)
#         - CCl<sub>4</sub>, CH<sub>2</sub>Cl<sub>2</sub>, CH<sub>3</sub>Br,
#           CH<sub>3</sub>CCl<sub>3</sub>, CH<sub>3</sub>Cl, CHCl<sub>3</sub>
# - ozone fluorinated compounds (23)
#     - HFCs (11)
#         - HFC-125, HFC-134a, HFC-143a, HFC-152a, HFC-227ea, HFC-23, HFC-236fa,
#           HFC-245fa, HFC-32, HFC-365mfc, HFC-4310mee
#     - PFCs (9)
#         - C<sub>2</sub>F<sub>6</sub>, C<sub>3</sub>F<sub>8</sub>,
#           C<sub>4</sub>F<sub>10</sub>, C<sub>5</sub>F<sub>12</sub>,
#           C<sub>6</sub>F<sub>14</sub>, C<sub>7</sub>F<sub>16</sub>,
#           C<sub>8</sub>F<sub>18</sub>, CC<sub>4</sub>F<sub>8</sub>,
#           CF<sub>4</sub>
#     - other (3)
#         - NF<sub>3</sub>, SF<sub>6</sub>, SO<sub>2</sub>F<sub>2</sub>
#
# ### Equivalent species
#
# For most models, you will not use all 43 species.
# As a result, we provide equivalent species too.
# There are two options if you don't want to use all 43 species.
#
# #### Option 1
#
# Use CO<sub>2</sub>, CH<sub>4</sub>, N<sub>2</sub>O and CFC-12 directly.
# Use CFC-11 equivalent to capture the radiative effect of all other species.
#
# #### Option 2
#
# Use CO<sub>2</sub>, CH<sub>4</sub> and N<sub>2</sub>O directly.
# Use CFC-12 equivalent
# to capture the radiative effect of all ozone depleting substances (ODSs)
# and HFC-134a equivalent
# to capture the radiative effect of all other fluorinated gases.

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ## Uncertainty
#
# At present, we provide no analysis of the uncertainty associated with these datasets.
#
# On top of the uncertainties in the historical data
# (which are, in radiative forcing terms, small),
# the scenario datasets are subject to uncertainties
# from the modelling process required to produce them.
# This means that, in radiative forcing terms, the uncertainty in these concentrations
# is much larger than the historical data.
# Nonetheless, it is very likely to be small
# compared to other uncertainties in the climate system,
# but this statement is not presently based on any robust analysis
# (rather it is based on expert judgement and past experience).

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ## Differences compared to CMIP6
#
# There are two major differences from CMIP6.
# The first is that the scenarios are different.
# By definition, this changes the concentrations.
# (There is also no 1:1 mapping between CMIP6 and CMIP7 scenarios,
# so the overall spacing and ensemble of scenarios is different too).
# The second is the transition from the history, observation-based period
# to the scenario, model-based projections.
# In CMIP6, this transition simply occured over a single year.
# In places, this led to notable changes in gradient at this transition point
# (further analysis of this is provided below).
# In CMIP7, we instead use a more sophisticated harmonisation algorithm
# (using the [gradient-aware-harmonisation](https://github.com/climate-resource/gradient-aware-harmonisation)
# package developed as part of the ESA project).
# This leads to smoother, more realistic transitions
# between the history, observation-based period
# and the scenario, model-based projections.

# %% [markdown] editable=true slideshow={"slide_type": ""}
# # User guide

# %% [markdown] editable=true slideshow={"slide_type": ""}
# Having downloaded the data, using it is relatively straightforward
# (scenario identification issue discussed above notwithstanding).

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ## Annual-, global-mean data
#
# We start with the annual-, global-mean data.
# Like all our datasets, this is composed of one file per scenario,
# given seven files in total.
#
# For yearly data, the time labels in the filename are years
# (for months, the month is included e.g. you will see `202201-210012`
# rather than `2022-2100` in the filename,
# the files also have different values for the `frequency` attribute).
# Global-mean data is identified by the 'grid label' `gm`,
# which appears in the filename.
# Below we show the filenames for the CO<sub>2</sub> output.
#
# **Note: in the draft datasets, the time axis starts in 2023.
# This will be updated to a 2022 start for the final datasets,
# in line with the rest of the scenario datasets.**

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
# Ensure data is downloaded
query_kwargs_co2_yearly_global = dict(
    ghg="co2",
    time_sampling="yr",
    grid="gm",
    cmip_era="CMIP6Plus",
    source_version="0.1.0",
    institution_id="CR",
    target_mip="ScenarioMIP",
    engine=engine,
)
fetch_and_load = partial(
    fetch_and_load_ghg_dataset_scenarios,
    local_data_root_dir=local_data_root_dir,
    # index_node=KnownIndexNode.DKRZ,
    # cmip_era="CMIP6",
    # source_id="UoM-CMIP-1-2-0",
    index_node=KnownIndexNode.ORNL,
)

_ = fetch_and_load(**query_kwargs_co2_yearly_global)

# Get file paths
co2_yearly_global_fps = get_ghg_dataset_local_files(**query_kwargs_co2_yearly_global)

# %% editable=true slideshow={"slide_type": ""} tags=["remove_input"]
for fp in sorted(co2_yearly_global_fps)[::-1]:
    print(f"- {fp.name}")

# %% [markdown] editable=true slideshow={"slide_type": ""}
# Output for other gases are named identically,
# with `co2` being replaced by the other gas name.
# For example, for methane the filenames are:

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
# Ensure data is downloaded
query_kwargs_ch4_yearly_global = {
    **query_kwargs_co2_yearly_global,
    "ghg": "ch4",
}
_ = fetch_and_load(**query_kwargs_ch4_yearly_global)

# %% editable=true slideshow={"slide_type": ""} tags=["remove_input"]
# Get file paths
ch4_yearly_global_fps = get_ghg_dataset_local_files(**query_kwargs_ch4_yearly_global)
for fp in sorted(ch4_yearly_global_fps)[::-1]:
    print(f"- {fp.name}")

# %% [markdown] editable=true slideshow={"slide_type": ""}
# As described above, the data is netCDF files.
# This means that metadata can be trivially inspected
# using a tool like `ncdump`.
# There is a lot of metadata included in these files.
# In general, you should not need to parse this metadata directly.
# However, if you have specific questions,
# please feel free to contact the emails given in the `contact` attribute.

# %% editable=true slideshow={"slide_type": ""} tags=["remove_input"]
# !ncdump -h {co2_yearly_global_fps[0]} | fold -w 80 -s

# %% [markdown] editable=true slideshow={"slide_type": ""}
# Using a tool like [xarray](https://github.com/pydata/xarray),
# loading and working with the data is trivial.

# %% editable=true slideshow={"slide_type": ""}
import xarray as xr

time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)
ds_example_co2_yearly_global = xr.open_dataset(
    co2_yearly_global_fps[-1], decode_times=time_coder
)

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
# Force values to compute to avoid dask getting involved
ds_example_co2_yearly_global = ds_example_co2_yearly_global.compute()

# %% editable=true slideshow={"slide_type": ""}
ds_example_co2_yearly_global

# %% editable=true slideshow={"slide_type": ""} tags=["remove_input"]
ds_example_co2_yearly_global["co2"].plot.scatter(alpha=0.4)
plt.show()

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ## Space- and time-average nature of the data
#
# All of our data represents the mean over each cell.
# This is indicated by the `cell_methods` attribute
# of all of our output variables.

# %% editable=true slideshow={"slide_type": ""}
ds_example_co2_yearly_global["co2"].attrs["cell_methods"]

# %% [markdown] editable=true slideshow={"slide_type": ""}
# This mean is both in space and time.
# The time bounds covered by each step
# are specified by the `time_bnds` variable
# (when there is spatial information,
# equivalent `lat_bnds` and `lon_bnds`
# information is also included).
# This variable specifies the start (inclusive)
# and end (exclusive) of the time period
# covered by each data point.

# %% editable=true slideshow={"slide_type": ""}
ds_example_co2_yearly_global["time_bnds"]

# %% [markdown] editable=true slideshow={"slide_type": ""}
# As a result of the time average that the data represents,
# it is inappropriate to plot this data
# using a line plot
# (the mean of the lines joining the points
# is not the same as the data given in the files).
# Instead, the data should be plotted (and used)
# as a scatter or a step plot, as shown below.
# (The same logic applies to any spatial plots
# which could be created from our datasets
# that include spatial dimensions).

# %% editable=true slideshow={"slide_type": ""} tags=["remove_input"]
ds_plt = ds_example_co2_yearly_global.isel(time=slice(None, 5))

fig, ax = plt.subplots(figsize=(8, 4))
ds_plt["co2"].plot.scatter(ax=ax)

for bounds, val in zip(ds_plt["time_bnds"].values, ds_plt["co2"].values):
    ax.plot(bounds, [val, val], color="tab:blue", linewidth=1.0, alpha=0.7)

xticks = [cftime.DatetimeGregorian(y, 1, 1) for y in range(2023, 2029)]
ax.set_xticks(xticks)
ax.set_xlim(xticks[0], xticks[-1])
ax.grid()

plt.show()

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ## Monthly-, global-mean data
#
# If you want to have information at a finer level
# of temporal detail, we also provide monthly files.
#
# For monthly data, the time labels in the filename are months.
# Below we show the filenames for the CO<sub>2</sub> output.

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
# Ensure data is downloaded
query_kwargs_co2_monthly_global = {
    **query_kwargs_co2_yearly_global,
    "time_sampling": "mon",
}
_ = fetch_and_load(**query_kwargs_co2_monthly_global)

# %% editable=true slideshow={"slide_type": ""} tags=["remove_input"]
# Get file paths
co2_monthly_global_fps = get_ghg_dataset_local_files(**query_kwargs_co2_monthly_global)
for fp in sorted(co2_monthly_global_fps)[::-1]:
    print(f"- {fp.name}")

# %% [markdown] editable=true slideshow={"slide_type": ""}
# Again, the data can be trivially loaded with [xarray](https://github.com/pydata/xarray).

# %% editable=true slideshow={"slide_type": ""}
ds_example_co2_monthly_global = xr.open_mfdataset(
    co2_monthly_global_fps[-1], decode_times=time_coder
)

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
# Force values to compute to avoid dask getting involved
ds_example_co2_monthly_global = ds_example_co2_monthly_global.compute()

# %% editable=true slideshow={"slide_type": ""}
ds_example_co2_monthly_global

# %% [markdown] editable=true slideshow={"slide_type": ""}
# For this data, the time bounds show that each point
# is the average a month, not a year.

# %% editable=true slideshow={"slide_type": ""}
ds_example_co2_monthly_global["time_bnds"]

# %% [markdown] editable=true slideshow={"slide_type": ""}
# As above, as a result of the time average that the data represents,
# it is inappropriate to plot this data using a line plot.
# Scatter or step plots should be used instead.

# %% editable=true slideshow={"slide_type": ""} tags=["remove_input"]
ds_plt = ds_example_co2_monthly_global.isel(time=slice(None, 5 * 12))

fig, ax = plt.subplots(figsize=(8, 4))
ds_plt["co2"].plot.scatter(ax=ax)

for bounds, val in zip(ds_plt["time_bnds"].values, ds_plt["co2"].values):
    ax.plot(bounds, [val, val], color="tab:blue", linewidth=1.0, alpha=0.7)

xticks = [cftime.DatetimeGregorian(y, 1, 1) for y in range(2023, 2029)]
ax.set_xticks(xticks)
ax.set_xlim(xticks[0], xticks[-1])
ax.grid()

plt.show()

# %% [markdown] editable=true slideshow={"slide_type": ""}
# The monthly data includes seasonality.
# Plotting the monthly and yearly data
# on the same axes makes particularly clear
# why a line plot is inappropriate.

# %% editable=true slideshow={"slide_type": ""} tags=["remove_input"]
fig, ax = plt.subplots(figsize=(8, 4))

for ds_plt, label, colour in (
    (
        ds_example_co2_monthly_global.isel(time=slice(None, 5 * 12)),
        "monthly",
        "tab:blue",
    ),
    (ds_example_co2_yearly_global.isel(time=slice(None, 5)), "yearly", "tab:orange"),
):
    ds_plt["co2"].plot.scatter(ax=ax, label=label, color=colour, s=10)

    for bounds, val in zip(ds_plt["time_bnds"].values, ds_plt["co2"].values):
        ax.plot(bounds, [val, val], color=colour, linewidth=1.0, alpha=0.7)

ax.legend()

xticks = [cftime.DatetimeGregorian(y, 1, 1) for y in range(2023, 2029)]
ax.set_xticks(xticks)
ax.set_xlim(xticks[0], xticks[-1])
ax.grid()

plt.show()

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ## Monthly-, latitudinally-resolved data
#
# We also provide data with spatial,
# specifically latitudinal, resolution.
# This data comes on a 15-degree latitudinal grid
# (see below for details of the grid and latitudinal bounds).
# These files are identified by the grid label `gnz`.
# We only provide these files with monthly resolution.
#
# For completeness, we note that we also provide hemispheric means.
# These are not shown here,
# but are identified by the grid label `gr1z`.
#
# Below we show the filenames for the latitudinally-resolved data
# for CO<sub>2</sub>

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
# Ensure data is downloaded
query_kwargs_co2_monthly_lat = {
    **query_kwargs_co2_yearly_global,
    "time_sampling": "mon",
    "grid": "gnz",
}
_ = fetch_and_load(**query_kwargs_co2_monthly_lat)

# %% editable=true slideshow={"slide_type": ""} tags=["remove_input"]
# Get file paths
co2_monthly_lat_fps = get_ghg_dataset_local_files(**query_kwargs_co2_monthly_lat)
for fp in sorted(co2_monthly_lat_fps)[::-1]:
    print(f"- {fp.name}")

# %% [markdown] editable=true slideshow={"slide_type": ""}
# Again, the data can be trivially loaded with [xarray](https://github.com/pydata/xarray).

# %% editable=true slideshow={"slide_type": ""}
ds_example_co2_monthly_lat = xr.open_mfdataset(
    co2_monthly_lat_fps[-1],
    decode_times=time_coder,
    data_vars=None,
    compat="no_conflicts",
)

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
# Force values to compute to avoid dask getting involved
ds_example_co2_monthly_lat = ds_example_co2_monthly_lat.compute()

# %% editable=true slideshow={"slide_type": ""}
ds_example_co2_monthly_lat

# %% [markdown] editable=true slideshow={"slide_type": ""}
# For this data, the latitudinal bounds show the area
# over which each point is the average.

# %% editable=true slideshow={"slide_type": ""}
ds_example_co2_monthly_lat["lat_bnds"]

# %% [markdown] editable=true slideshow={"slide_type": ""}
# As above, but this time for the spatial axis,
# it is inappropriate to plot this data using a line plot.
# Scatter or step plots should be used instead.

# %% editable=true slideshow={"slide_type": ""} tags=["remove_input"]
ds_plt = ds_example_co2_monthly_lat.isel(time=slice(None, 12))


def get_label_for_month(inds: xr.Dataset) -> str:
    """
    Get the label for a given month of data
    """
    year = int(inds["time"].dt.year)
    month_name = calendar.month_name[int(inds["time"].dt.month)]

    return f"{year} - {month_name}"


mosaic_flat = [get_label_for_month(ds_plt.sel(time=time)) for time in ds_plt["time"]]

mosaic = [mosaic_flat[3 * i : 3 * (i + 1)] for i in range(len(mosaic_flat) // 3)]

fig, axes_d = plt.subplot_mosaic(mosaic, figsize=(8, 9), sharey=True, sharex=True)

for time in ds_plt["time"]:
    ds_plt_time = ds_plt.sel(time=time)
    label = get_label_for_month(ds_plt_time)

    axes_d[label].scatter(
        x=ds_plt_time["co2"].values,
        y=ds_plt_time["lat"].values,
        s=10,
        label=label,
    )

    for bounds, val in zip(ds_plt_time["lat_bnds"].values, ds_plt_time["co2"].values):
        axes_d[label].plot(
            [val, val], bounds, color="tab:blue", linewidth=1.0, alpha=0.7
        )

    yticks = np.arange(-90, 91, 15.0)
    axes_d[label].set_yticks(yticks)
    axes_d[label].set_ylim(yticks[0], yticks[-1])
    # axes_d[label].set_ylabel("Latitude (degrees north)")

    # axes_d[label].set_xlabel("co2 [ppm]")
    axes_d[label].grid()
    axes_d[label].set_title(label, fontsize="small")

for month in [1, 4, 7, 10]:
    axes_d[f"2023 - {calendar.month_name[month]}"].set_ylabel(
        "Latitude (degrees north)"
    )

for month in range(10, 13):
    axes_d[f"2023 - {calendar.month_name[month]}"].set_xlabel("co2 [ppm]")
# ax.legend(loc="center left", bbox_to_anchor=(1.05, 0.5))

plt.tight_layout()
plt.show()

# %% [markdown] editable=true slideshow={"slide_type": ""}
# We can compare the global-mean data
# to the data at each latitude.
# The strength of the latitudinal gradient varies also by gas (not shown).

# %% editable=true slideshow={"slide_type": ""} tags=["remove_input"]
fig, ax = plt.subplots(figsize=(8, 4))

time_slice = slice(None, 5 * 12)

ds_plt = ds_example_co2_monthly_global.isel(time=time_slice)
ds_plt["co2"].plot.scatter(
    ax=ax, label="global-mean", color="tab:blue", s=30, zorder=10.0
)

for bounds, val in zip(ds_plt["time_bnds"].values, ds_plt["co2"].values):
    ax.plot(bounds, [val, val], color="tab:blue", linewidth=1.0, alpha=0.7)

ds_all_lats = ds_example_co2_monthly_lat.isel(time=time_slice)

for i, lat in enumerate(sorted(ds_all_lats["lat"])[::-1]):
    ds_plt = ds_all_lats.sel(lat=lat)
    colour = matplotlib.colormaps["magma"](i / len(ds_all_lats["lat"]))

    ds_plt["co2"].plot.scatter(
        ax=ax, label=f"{float(lat)}", color=colour, marker="x", s=10
    )

    for bounds, val in zip(ds_plt["time_bnds"].values, ds_plt["co2"].values):
        ax.plot(bounds, [val, val], color=colour, linewidth=1.0, alpha=0.7)

ax.legend(loc="center left", bbox_to_anchor=(1.05, 0.5))

xticks = [cftime.DatetimeGregorian(y, 1, 1) for y in range(2023, 2029)]
ax.set_xticks(xticks)
ax.set_xlim(xticks[0], xticks[-1])
ax.grid()
ax.set_title("")

plt.show()

# %% [markdown] editable=true slideshow={"slide_type": ""}
# The data can also be plotted in a so-called "magic carpet"
# to see the variation in space and time simultaneously.

# %% editable=true slideshow={"slide_type": ""} tags=["remove_input"]
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(projection="3d")

tmp = ds_example_co2_monthly_lat["co2"].isel(time=range(10 * 12)).copy()
tmp = tmp.assign_coords(time=tmp["time"].dt.year + tmp["time"].dt.month / 12)
# Interpolate so the plot shows the step nature
tmp = tmp.interp(
    coords=dict(
        time=np.linspace(
            tmp["time"].values[0], tmp["time"].values[-1], tmp["time"].size * 10
        )
    ),
    method="nearest",
).interp(
    coords=dict(
        lat=np.linspace(
            tmp["lat"].values[0], tmp["lat"].values[-1], tmp["lat"].size * 10
        )
    ),
    method="nearest",
)

tmp.plot.surface(
    x="time",
    y="lat",
    ax=ax,
    cmap="magma_r",
    levels=30,
    # alpha=0.7,
)

ax.view_init(15, -135, 0)

plt.tight_layout()
plt.show()

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ## Transition from history
#
# Each dataset transitions smoothly from the historical data.

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ### Annual-, global-mean data

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
query_kwargs_co2_yearly_global_history = {
    **query_kwargs_co2_yearly_global,
    "cmip_era": "CMIP7",
    "source_id": "CR-CMIP-1-0-0",
    "source_version": "1.0.0",
    "target_mip": "CMIP",
}

fetch_and_load_history = partial(
    fetch_and_load_ghg_dataset,
    local_data_root_dir=local_data_root_dir,
    # index_node=KnownIndexNode.DKRZ,
    index_node=KnownIndexNode.ORNL,
)

ds_history_co2_yearly_global = fetch_and_load_history(
    **query_kwargs_co2_yearly_global_history
).compute()

# %% editable=true slideshow={"slide_type": ""} tags=["remove_input"]
fig, ax = plt.subplots()


def time_sel_func(
    da: xr.DataArray, start_year: int, end_year: int
) -> npt.NDArray[bool]:
    """
    Select times of interest
    """
    return (da.dt.year >= start_year) & (da.dt.year <= end_year)


ds_history_co2_yearly_global["co2"].sel(
    time=time_sel_func(
        ds_history_co2_yearly_global["time"], start_year=2010, end_year=2030
    )
).plot.scatter(
    ax=ax,
    edgecolors="none",
    label=f"{ds_history_co2_yearly_global.attrs['source_id']} (history)",
)
ds_example_co2_yearly_global["co2"].sel(
    time=time_sel_func(
        ds_example_co2_yearly_global["time"], start_year=2010, end_year=2030
    )
).plot.scatter(
    ax=ax,
    edgecolors="none",
    label=f"{ds_example_co2_yearly_global.attrs['source_id']} (scenario)",
)

ax.set_xlim(cftime.DatetimeGregorian(2010, 1, 1), cftime.DatetimeGregorian(2030, 1, 1))
ax.legend()

plt.show()

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ### Monthly-, global-mean data
#
# Note that the transition from history to scenarios
# is clearly wrong in the draft dataset.
# This will be fixed before the final dataset is published.

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
query_kwargs_co2_monthly_global_history = {
    **query_kwargs_co2_yearly_global_history,
    "time_sampling": "mon",
}

ds_history_co2_monthly_global = fetch_and_load_history(
    **query_kwargs_co2_monthly_global_history
).compute()

# %% editable=true slideshow={"slide_type": ""} tags=["remove_input"]
fig, ax = plt.subplots()


def time_sel_func(
    da: xr.DataArray, start_year: int, end_year: int
) -> npt.NDArray[bool]:
    """
    Select times of interest
    """
    return (da.dt.year >= start_year) & (da.dt.year <= end_year)


ds_history_co2_monthly_global["co2"].sel(
    time=time_sel_func(
        ds_history_co2_monthly_global["time"], start_year=2010, end_year=2030
    )
).plot.scatter(
    ax=ax,
    edgecolors="none",
    label=f"{ds_history_co2_yearly_global.attrs['source_id']} (history)",
)
ds_example_co2_monthly_global["co2"].sel(
    time=time_sel_func(
        ds_example_co2_monthly_global["time"], start_year=2010, end_year=2030
    )
).plot.scatter(
    ax=ax,
    edgecolors="none",
    label=f"{ds_example_co2_yearly_global.attrs['source_id']} (scenario)",
)

ax.set_xlim(cftime.DatetimeGregorian(2010, 1, 1), cftime.DatetimeGregorian(2030, 1, 1))
ax.legend()

plt.show()

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ### Monthly-, 15-degree latitudinally binned data
#
# Note that the transition from history to scenarios
# is clearly wrong in the draft dataset.
# This will be fixed before the final dataset is published.

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
# Get file paths
query_kwargs_co2_monthly_lat_history = {
    **query_kwargs_co2_monthly_global_history,
    "grid": "gnz",
}

ds_history_co2_monthly_lat_history = fetch_and_load_history(
    **query_kwargs_co2_monthly_lat_history
).compute()


# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
def to_data_frame(
    ds: xr.Dataset,
    unstack_col: str,
    assign_metadata: dict[str, Any] | None,
    ds_var: str | None = None,
) -> pd.DataFrame:
    """
    Convert an [xr.Dataset][xarray.Dataset] to [pd.DataFrame][pandas.DataFrame]

    Definitely not a general function,
    just a helper for the specific kinds of conversions we want to do here.
    """
    if ds_var is None:
        ds_var = ds.attrs["variable_id"]

    res = (
        ds[ds_var]
        .to_dataframe()[ds_var]
        .unstack(unstack_col)
        .pix.assign(unit=ds[ds_var].attrs["units"])
    )
    if assign_metadata is not None:
        res = res.pix.assign(**assign_metadata)

    res = res.openscm.eiim()

    return res


# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
pdf_l = [
    to_data_frame(
        ds_example_co2_monthly_lat,
        unstack_col="time",
        assign_metadata={"experiment": "scenario"},
    ),
    to_data_frame(
        ds_history_co2_monthly_lat_history,
        unstack_col="time",
        assign_metadata={"experiment": "historical"},
    ),
]

pdf = pd.concat([v.reorder_levels(pdf_l[0].index.names) for v in pdf_l]).sort_index(
    axis="columns"
)
pdf.columns = [v.year + (v.month * 2 - 1) / 24 for v in pdf.columns]
pdf

# %% editable=true slideshow={"slide_type": ""} tags=["remove_input"]
start_year = 2015
end_year = 2030

sns.relplot(
    data=pdf.loc[:, start_year:end_year].openscm.to_long_data(),
    x="time",
    y="value",
    col="lat",
    col_wrap=3,
    col_order=sorted(pdf.index.get_level_values("lat").unique())[::-1],
    kind="scatter",
    hue="experiment",
)

plt.show()

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ## Full scenario set
#
# Having seen the transition for a single scenario,
# we now show the full scenario set
# including the transition from history.

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ### Annual-, global-mean data

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
gases_to_show = ["co2", "ch4", "n2o", "cfc12eq", "hfc134aeq"]
ds_gases_full_yearly_d = {}
for gas in gases_to_show:
    ds_gases_full_yearly_d[gas] = {}
    for key, target_mip, source_version, institution_id, cmip_era in (
        ("history", "CMIP", "1.0.0", "CR", "CMIP7"),
        ("scenarios", "ScenarioMIP", "0.1.0", "CR", "CMIP6Plus"),
    ):
        query_kwargs = {
            "ghg": gas,
            "time_sampling": "yr",
            "grid": "gm",
            "target_mip": target_mip,
            "source_version": source_version,
            "institution_id": institution_id,
            "cmip_era": cmip_era,
            "engine": engine,
        }
        ds = fetch_and_load(**query_kwargs)

        # Unify time axis days to simplify
        ds["time"] = [
            cftime.DatetimeProlepticGregorian(v.year, v.month, 15)
            for v in ds["time"].values
        ]

        # compute to avoid dask weirdness
        ds_gases_full_yearly_d[gas][key] = ds.compute()

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
palette = {
    "history": "k",
    "history-cmip6": "tab:grey",
    "historical": "k",
    "historical-cmip6": "tab:grey",
    "vl": "#24a4ff",
    "vllo": "#24a4ff",
    "ln": "#4a0daf",
    "vlho": "#4a0daf",
    "l": "#00cc69",
    "ml": "#f5ac00",
    "m": "#ffa9dc",
    "h": "#700000",
    "hl": "#8f003b",
    "ssp119": "#00a9cf",
    "ssp126": "#003466",
    "ssp245": "#f69320",
    "ssp370": "#df0000",
    "ssp434": "#2274ae",
    "ssp460": "#b0724e",
    "ssp534-over": "#92397a",
    "ssp585": "#980002",
}

plt.rcParams["axes.xmargin"] = 0

hue_order = [
    "h",
    "hl",
    "m",
    "ml",
    "l",
    # "ln",
    "vlho",
    # "vl",
    "vllo",
    "historical",
]

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
start_year = 1950

pdf_l = [
    v.sort_index(axis="columns")
    for ghg in tqdm.auto.tqdm(ds_gases_full_yearly_d)
    for v in [
        to_data_frame(
            ds=ds_gases_full_yearly_d[ghg]["scenarios"]
            .sel(
                time=ds_gases_full_yearly_d[ghg]["scenarios"]["time"].dt.year
                >= start_year
            )
            .groupby("time.year")
            .mean(),
            unstack_col="year",
            assign_metadata={"ghg": ghg},
        ).openscm.update_index_levels_from_other(
            {"experiment": ("scenario", lambda x: x)}
        ),
        to_data_frame(
            ds=ds_gases_full_yearly_d[ghg]["history"]
            .sel(
                time=ds_gases_full_yearly_d[ghg]["history"]["time"].dt.year
                >= start_year
            )
            .groupby("time.year")
            .mean(),
            unstack_col="year",
            assign_metadata={"experiment": "historical", "ghg": ghg},
        ),
    ]
]
pdf = pd.concat(
    [v.reorder_levels(pdf_l[0].index.names) for v in pdf_l], axis="rows"
).sort_index(axis="columns")

# pdf

# %% editable=true slideshow={"slide_type": ""} tags=["remove_input"]
sns.relplot(
    data=pdf.openscm.to_long_data(),
    x="time",
    y="value",
    hue="experiment",
    palette=palette,
    hue_order=hue_order,
    kind="scatter",
    col="ghg",
    col_wrap=3,
    facet_kws=dict(sharey=False),
)

plt.show()

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ### Monthly-, global-mean data

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
ds_gases_full_monthly_d = {}
for gas in gases_to_show:
    ds_gases_full_monthly_d[gas] = {}
    for key, target_mip, source_version, institution_id, cmip_era in (
        ("history", "CMIP", "1.0.0", "CR", "CMIP7"),
        ("scenarios", "ScenarioMIP", "0.1.0", "CR", "CMIP6Plus"),
    ):
        query_kwargs = {
            "ghg": gas,
            "time_sampling": "mon",
            "grid": "gm",
            "target_mip": target_mip,
            "source_version": source_version,
            "institution_id": institution_id,
            "cmip_era": cmip_era,
            "engine": engine,
        }
        ds = fetch_and_load(**query_kwargs)

        # Unify time axis days to simplify
        ds["time"] = [
            cftime.DatetimeProlepticGregorian(v.year, v.month, 15)
            for v in ds["time"].values
        ]

        # compute to avoid dask weirdness
        ds_gases_full_monthly_d[gas][key] = ds.compute()

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
start_year = 2015
end_year = 2030

pdf_l = [
    v.sort_index(axis="columns")
    for ghg in tqdm.auto.tqdm(ds_gases_full_monthly_d)
    for v in [
        to_data_frame(
            ds=ds_gases_full_monthly_d[ghg]["scenarios"].sel(
                time=(
                    (
                        ds_gases_full_monthly_d[ghg]["scenarios"]["time"].dt.year
                        >= start_year
                    )
                    & (
                        ds_gases_full_monthly_d[ghg]["scenarios"]["time"].dt.year
                        <= end_year
                    )
                )
            ),
            unstack_col="time",
            assign_metadata={"ghg": ghg},
        ).openscm.update_index_levels_from_other(
            {"experiment": ("scenario", lambda x: x)}
        ),
        to_data_frame(
            ds=ds_gases_full_monthly_d[ghg]["history"].sel(
                time=(
                    (
                        ds_gases_full_monthly_d[ghg]["history"]["time"].dt.year
                        >= start_year
                    )
                    & (
                        ds_gases_full_monthly_d[ghg]["history"]["time"].dt.year
                        <= end_year
                    )
                )
            ),
            unstack_col="time",
            assign_metadata={"experiment": "historical", "ghg": ghg},
        ),
    ]
]
pdf = pd.concat(
    [v.reorder_levels(pdf_l[0].index.names) for v in pdf_l], axis="rows"
).sort_index(axis="columns")
pdf.columns = [v.year + (v.month * 2 - 1) / 24 for v in pdf.columns]

# pdf

# %% editable=true slideshow={"slide_type": ""} tags=["remove_input"]
sns.relplot(
    data=pdf.openscm.to_long_data(),
    x="time",
    y="value",
    hue="experiment",
    palette=palette,
    hue_order=hue_order,
    kind="scatter",
    col="ghg",
    col_wrap=3,
    edgecolors="none",
    facet_kws=dict(sharey=False),
)

plt.show()

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ### Monthly-, 15-degree latitudinal data

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
ds_gases_full_monthly_lat_d = {}
for gas in gases_to_show:
    ds_gases_full_monthly_lat_d[gas] = {}
    for key, target_mip, source_version, institution_id, cmip_era in (
        ("history", "CMIP", "1.0.0", "CR", "CMIP7"),
        ("scenarios", "ScenarioMIP", "0.1.0", "CR", "CMIP6Plus"),
    ):
        query_kwargs = {
            "ghg": gas,
            "time_sampling": "mon",
            "grid": "gnz",
            "target_mip": target_mip,
            "source_version": source_version,
            "institution_id": institution_id,
            "cmip_era": cmip_era,
            "engine": engine,
        }
        ds = fetch_and_load(**query_kwargs)

        # Unify time axis days to simplify
        ds["time"] = [
            cftime.DatetimeProlepticGregorian(v.year, v.month, 15)
            for v in ds["time"].values
        ]

        # compute to avoid dask weirdness
        ds_gases_full_monthly_lat_d[gas][key] = ds.compute()

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
start_year = 2015
end_year = 2030

pdf_l = [
    v.sort_index(axis="columns")
    for ghg in tqdm.auto.tqdm(ds_gases_full_monthly_lat_d)
    for v in [
        to_data_frame(
            ds=ds_gases_full_monthly_lat_d[ghg]["scenarios"].sel(
                time=(
                    (
                        ds_gases_full_monthly_lat_d[ghg]["scenarios"]["time"].dt.year
                        >= start_year
                    )
                    & (
                        ds_gases_full_monthly_lat_d[ghg]["scenarios"]["time"].dt.year
                        <= end_year
                    )
                )
            ),
            unstack_col="time",
            assign_metadata={"ghg": ghg},
        ).openscm.update_index_levels_from_other(
            {"experiment": ("scenario", lambda x: x)}
        ),
        to_data_frame(
            ds=ds_gases_full_monthly_lat_d[ghg]["history"].sel(
                time=(
                    (
                        ds_gases_full_monthly_lat_d[ghg]["history"]["time"].dt.year
                        >= start_year
                    )
                    & (
                        ds_gases_full_monthly_lat_d[ghg]["history"]["time"].dt.year
                        <= end_year
                    )
                )
            ),
            unstack_col="time",
            assign_metadata={"experiment": "historical", "ghg": ghg},
        ),
    ]
]
pdf = pd.concat(
    [v.reorder_levels(pdf_l[0].index.names) for v in pdf_l], axis="rows"
).sort_index(axis="columns")
pdf.columns = [v.year + (v.month * 2 - 1) / 24 for v in pdf.columns]

# pdf

# %% editable=true slideshow={"slide_type": ""} tags=["remove_input"]
sns.relplot(
    data=pdf.loc[pix.isin(lat=[-82.5, 7.5, 82.5])].openscm.to_long_data(),
    x="time",
    y="value",
    hue="experiment",
    palette=palette,
    hue_order=hue_order,
    kind="scatter",
    row="ghg",
    col="lat",
    edgecolors="none",
    facet_kws=dict(sharey="row"),
)

plt.show()

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ## Differences from CMIP6
#
# ### File formats and naming
#
# The file formats are generally close to CMIP6.
# There are three key changes:
#
# 1. the global-mean and hemispheric-mean data are split into separate files.
#    In CMIP6, this data was in the same file (with a grid label of `GMNHSH`).
#    We have split this for two reasons:
#    a) `GMNHSH` is not a grid label recognised in the CMIP CVs and
#    b) having global-mean and hemispheric-mean data in the same file
#       required us to introduce a 'sector' coordinate,
#       which was confusing and does not follow the CF-conventions.
# 1. the files are split into different time components.
#    The CMIP6 data had the scenarios and their extensions in a single file.
#    The CMIP7 extensions are not defined yet,
#    so the scenarios (up to 2100) will be in one file,
#    with the extensions being in a separate file
#    (and under separate source IDs).
# 1. the variable names were simplified.
#    They are now simply the names of the gases,
#    for example we now use "co2" rather than "mole_fraction_of_carbon_dioxide".
#    A full mapping is provided below.
#
# #### Variable name mapping
#
# ```python
# CMIP6_TO_CMIP7_VARIABLE_MAP = {
#     # name in CMIP6: name in CMIP7
#     "mole_fraction_of_carbon_dioxide_in_air": "co2",
#     "mole_fraction_of_methane_in_air": "ch4",
#     "mole_fraction_of_nitrous_oxide_in_air": "n2o",
#     "mole_fraction_of_c2f6_in_air": "c2f6",
#     "mole_fraction_of_c3f8_in_air": "c3f8",
#     "mole_fraction_of_c4f10_in_air": "c4f10",
#     "mole_fraction_of_c5f12_in_air": "c5f12",
#     "mole_fraction_of_c6f14_in_air": "c6f14",
#     "mole_fraction_of_c7f16_in_air": "c7f16",
#     "mole_fraction_of_c8f18_in_air": "c8f18",
#     "mole_fraction_of_c_c4f8_in_air": "cc4f8",
#     "mole_fraction_of_carbon_tetrachloride_in_air": "ccl4",
#     "mole_fraction_of_cf4_in_air": "cf4",
#     "mole_fraction_of_cfc11_in_air": "cfc11",
#     "mole_fraction_of_cfc113_in_air": "cfc113",
#     "mole_fraction_of_cfc114_in_air": "cfc114",
#     "mole_fraction_of_cfc115_in_air": "cfc115",
#     "mole_fraction_of_cfc12_in_air": "cfc12",
#     "mole_fraction_of_ch2cl2_in_air": "ch2cl2",
#     "mole_fraction_of_methyl_bromide_in_air": "ch3br",
#     "mole_fraction_of_ch3ccl3_in_air": "ch3ccl3",
#     "mole_fraction_of_methyl_chloride_in_air": "ch3cl",
#     "mole_fraction_of_chcl3_in_air": "chcl3",
#     "mole_fraction_of_halon1211_in_air": "halon1211",
#     "mole_fraction_of_halon1301_in_air": "halon1301",
#     "mole_fraction_of_halon2402_in_air": "halon2402",
#     "mole_fraction_of_hcfc141b_in_air": "hcfc141b",
#     "mole_fraction_of_hcfc142b_in_air": "hcfc142b",
#     "mole_fraction_of_hcfc22_in_air": "hcfc22",
#     "mole_fraction_of_hfc125_in_air": "hfc125",
#     "mole_fraction_of_hfc134a_in_air": "hfc134a",
#     "mole_fraction_of_hfc143a_in_air": "hfc143a",
#     "mole_fraction_of_hfc152a_in_air": "hfc152a",
#     "mole_fraction_of_hfc227ea_in_air": "hfc227ea",
#     "mole_fraction_of_hfc23_in_air": "hfc23",
#     "mole_fraction_of_hfc236fa_in_air": "hfc236fa",
#     "mole_fraction_of_hfc245fa_in_air": "hfc245fa",
#     "mole_fraction_of_hfc32_in_air": "hfc32",
#     "mole_fraction_of_hfc365mfc_in_air": "hfc365mfc",
#     "mole_fraction_of_hfc4310mee_in_air": "hfc4310mee",
#     "mole_fraction_of_nf3_in_air": "nf3",
#     "mole_fraction_of_sf6_in_air": "sf6",
#     "mole_fraction_of_so2f2_in_air": "so2f2",
#     "mole_fraction_of_cfc11eq_in_air": "cfc11eq",
#     "mole_fraction_of_cfc12eq_in_air": "cfc12eq",
#     "mole_fraction_of_hfc134aeq_in_air": "hfc134aeq",
# }
# ```

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ### Data comparisons
#
# Comparing the data from CMIP6 and CMIP7 shows changes in two key areas:
#
# 1. the scenarios are simply different
# 2. the transition from historical to scenarios is more carefully harmonised

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ### Annual-, global-mean data

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
scenario_group_map = {
    "h": "high",
    "hl": "high",
    "ssp585": "high",
    "ssp534-over": "high",
    "ssp370": "high",
    "ssp460": "high",
    "m": "continuing-trends",
    "ml": "continuing-trends",
    "ssp434": "continuing-trends",
    "ssp245": "continuing-trends",
    "l": "low",
    # "vl": "low",
    # "ln": "low",
    "vllo": "low",
    "vlho": "low",
    "ssp126": "low",
    "ssp119": "low",
    "historical": "historical",
}

hue_order_incl_cmip6 = [
    "h",
    "hl",
    "m",
    "ml",
    "l",
    # "ln",
    "vlho",
    # "vl",
    "vllo",
    "historical",
    "ssp585",
    "ssp534-over",
    "ssp370",
    "ssp460",
    "ssp434",
    "ssp245",
    "ssp126",
    "ssp119",
]

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
ds_gases_full_yearly_multi_phase_d = {}
for gas in gases_to_show:
    ds_gases_full_yearly_multi_phase_d[gas] = {}
    for key, target_mip, source_version, institution_id, cmip_era in (
        ("history-cmip7", "CMIP", "1.0.0", "CR", "CMIP7"),
        ("scenarios-cmip7", "ScenarioMIP", "0.1.0", "CR", "CMIP6Plus"),
        ("history-cmip6", "CMIP", "1.2.0", "UoM", "CMIP6"),
        ("scenarios-cmip6", "ScenarioMIP", "1.2.1", "UoM", "CMIP6"),
    ):
        query_kwargs = {
            "ghg": gas,
            "time_sampling": "yr",
            "grid": "gm",
            "target_mip": target_mip,
            "source_version": source_version,
            "institution_id": institution_id,
            "cmip_era": cmip_era,
            "engine": engine,
        }
        ds = fetch_and_load(**query_kwargs)

        # Unify time axis days to simplify
        ds["time"] = [
            cftime.DatetimeProlepticGregorian(v.year, v.month, 15)
            for v in ds["time"].values
        ]

        # compute to avoid dask weirdness
        ds_gases_full_yearly_multi_phase_d[gas][key] = ds.compute()

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
start_year = 2000
end_year = 2100
pdf_l = [
    v.sort_index(axis="columns")
    for ghg in tqdm.auto.tqdm(ds_gases_full_yearly_multi_phase_d)
    for v in [
        to_data_frame(
            ds=ds_gases_full_yearly_multi_phase_d[ghg]["scenarios-cmip7"]
            .sel(
                time=(
                    (
                        ds_gases_full_yearly_multi_phase_d[ghg]["scenarios-cmip7"][
                            "time"
                        ].dt.year
                        >= start_year
                    )
                    & (
                        ds_gases_full_yearly_multi_phase_d[ghg]["scenarios-cmip7"][
                            "time"
                        ].dt.year
                        <= end_year
                    )
                )
            )
            .groupby("time.year")
            .mean(),
            unstack_col="year",
            assign_metadata={"ghg": ghg, "cmip_era": "CMIP7"},
        ).openscm.update_index_levels_from_other(
            {"experiment": ("scenario", lambda x: x)}
        ),
        to_data_frame(
            ds=ds_gases_full_yearly_multi_phase_d[ghg]["history-cmip7"]
            .sel(
                time=(
                    (
                        ds_gases_full_yearly_multi_phase_d[ghg]["history-cmip7"][
                            "time"
                        ].dt.year
                        >= start_year
                    )
                    & (
                        ds_gases_full_yearly_multi_phase_d[ghg]["history-cmip7"][
                            "time"
                        ].dt.year
                        <= end_year
                    )
                )
            )
            .groupby("time.year")
            .mean(),
            unstack_col="year",
            assign_metadata={
                "experiment": "historical",
                "ghg": ghg,
                "cmip_era": "CMIP7",
            },
        ),
        to_data_frame(
            ds=ds_gases_full_yearly_multi_phase_d[ghg]["scenarios-cmip6"]
            .sel(
                time=(
                    (
                        ds_gases_full_yearly_multi_phase_d[ghg]["scenarios-cmip6"][
                            "time"
                        ].dt.year
                        >= start_year
                    )
                    & (
                        ds_gases_full_yearly_multi_phase_d[ghg]["scenarios-cmip6"][
                            "time"
                        ].dt.year
                        <= end_year
                    )
                )
            )
            .groupby("time.year")
            .mean(),
            unstack_col="year",
            assign_metadata={"ghg": ghg, "cmip_era": "CMIP6"},
            ds_var=ghg,
        ).openscm.update_index_levels_from_other(
            {"experiment": ("scenario", lambda x: x)}
        ),
        to_data_frame(
            ds=ds_gases_full_yearly_multi_phase_d[ghg]["history-cmip6"]
            .sel(
                time=(
                    (
                        ds_gases_full_yearly_multi_phase_d[ghg]["history-cmip6"][
                            "time"
                        ].dt.year
                        >= start_year
                    )
                    & (
                        ds_gases_full_yearly_multi_phase_d[ghg]["history-cmip6"][
                            "time"
                        ].dt.year
                        <= end_year
                    )
                )
            )
            .groupby("time.year")
            .mean(),
            unstack_col="year",
            assign_metadata={
                "experiment": "historical",
                "ghg": ghg,
                "cmip_era": "CMIP6",
            },
            ds_var=ghg,
        ),
    ]
]
pdf = pd.concat(
    [v.reorder_levels(pdf_l[0].index.names) for v in pdf_l], axis="rows"
).sort_index(axis="columns")

pdf = pdf.openscm.update_index_levels_from_other(
    {"scenario_group": ("experiment", lambda x: scenario_group_map[x])}
)
pdf_grouped = {sg: sgdf for sg, sgdf in pdf.groupby("scenario_group")}
tmp_l = []
for sg, sgdf in pdf_grouped.items():
    if sg == "historical":
        continue

    tmp_l.append(sgdf)
    tmp_l.append(
        pdf_grouped["historical"].openscm.set_index_levels({"scenario_group": sg})
    )

pdf = pd.concat([v.reorder_levels(tmp_l[0].index.names) for v in tmp_l])
pdf

# %% editable=true slideshow={"slide_type": ""} tags=["remove_input"]
start_year = 2000
end_year = 2100

sns.relplot(
    data=pdf.loc[
        :, (pdf.columns >= start_year) & (pdf.columns <= end_year + 1)
    ].openscm.to_long_data(),
    x="time",
    y="value",
    hue="experiment",
    palette=palette,
    style="cmip_era",
    # markers={"CMIP6": ".", "CMIP7": "o"},
    # edgecolors="none",
    markers={"CMIP6": "+", "CMIP7": 4},
    hue_order=hue_order_incl_cmip6,
    kind="scatter",
    row="ghg",
    col="scenario_group",
    col_order=["low", "continuing-trends", "high"],
    facet_kws=dict(sharey=False),
    s=100,
    alpha=0.8,
)

plt.show()

# %% [markdown] editable=true slideshow={"slide_type": ""}
# If we zoom in, the harmonisation issues in the CMIP6
# data are more clear, most obviously for CH<sub>4</sub>
# where the transition from history to scenarios
# also includes a notable increase in the rate of change
# that is independent of scenario
# (i.e. is not driven by a change in emissions).

# %% editable=true slideshow={"slide_type": ""} tags=["remove_input"]
start_year = 2005
end_year = 2030

sns.relplot(
    data=pdf.loc[
        :, (pdf.columns >= start_year) & (pdf.columns <= end_year + 1)
    ].openscm.to_long_data(),
    x="time",
    y="value",
    hue="experiment",
    palette=palette,
    style="cmip_era",
    # markers={"CMIP6": ".", "CMIP7": "o"},
    # edgecolors="none",
    markers={"CMIP6": "+", "CMIP7": 4},
    hue_order=hue_order_incl_cmip6,
    kind="scatter",
    row="ghg",
    col="scenario_group",
    col_order=["low", "continuing-trends", "high"],
    facet_kws=dict(sharey="row"),
    s=100,
    alpha=0.8,
)

plt.show()

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ### Monthly-, global-mean data

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
ds_gases_full_monthly_multi_phase_d = {}
for gas in gases_to_show:
    ds_gases_full_monthly_multi_phase_d[gas] = {}
    for key, target_mip, source_version, institution_id, cmip_era in (
        ("history-cmip7", "CMIP", "1.0.0", "CR", "CMIP7"),
        ("scenarios-cmip7", "ScenarioMIP", "0.1.0", "CR", "CMIP6Plus"),
        ("history-cmip6", "CMIP", "1.2.0", "UoM", "CMIP6"),
        ("scenarios-cmip6", "ScenarioMIP", "1.2.1", "UoM", "CMIP6"),
    ):
        query_kwargs = {
            "ghg": gas,
            "time_sampling": "mon",
            "grid": "gm",
            "target_mip": target_mip,
            "source_version": source_version,
            "institution_id": institution_id,
            "cmip_era": cmip_era,
            "engine": engine,
        }
        ds = fetch_and_load(**query_kwargs)

        # Unify time axis days to simplify
        ds["time"] = [
            cftime.DatetimeProlepticGregorian(v.year, v.month, 15)
            for v in ds["time"].values
        ]

        # compute to avoid dask weirdness
        ds_gases_full_monthly_multi_phase_d[gas][key] = ds.compute()

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
start_year = 2000
end_year = 2100
pdf_l = [
    v.sort_index(axis="columns")
    for ghg in tqdm.auto.tqdm(ds_gases_full_monthly_multi_phase_d)
    for v in [
        to_data_frame(
            ds=ds_gases_full_monthly_multi_phase_d[ghg]["scenarios-cmip7"].sel(
                time=(
                    (
                        ds_gases_full_monthly_multi_phase_d[ghg]["scenarios-cmip7"][
                            "time"
                        ].dt.year
                        >= start_year
                    )
                    & (
                        ds_gases_full_monthly_multi_phase_d[ghg]["scenarios-cmip7"][
                            "time"
                        ].dt.year
                        <= end_year
                    )
                )
            ),
            unstack_col="time",
            assign_metadata={"ghg": ghg, "cmip_era": "CMIP7"},
        ).openscm.update_index_levels_from_other(
            {"experiment": ("scenario", lambda x: x)}
        ),
        to_data_frame(
            ds=ds_gases_full_monthly_multi_phase_d[ghg]["history-cmip7"].sel(
                time=(
                    (
                        ds_gases_full_monthly_multi_phase_d[ghg]["history-cmip7"][
                            "time"
                        ].dt.year
                        >= start_year
                    )
                    & (
                        ds_gases_full_monthly_multi_phase_d[ghg]["history-cmip7"][
                            "time"
                        ].dt.year
                        <= end_year
                    )
                )
            ),
            unstack_col="time",
            assign_metadata={
                "experiment": "historical",
                "ghg": ghg,
                "cmip_era": "CMIP7",
            },
        ),
        to_data_frame(
            ds=ds_gases_full_monthly_multi_phase_d[ghg]["scenarios-cmip6"].sel(
                time=(
                    (
                        ds_gases_full_monthly_multi_phase_d[ghg]["scenarios-cmip6"][
                            "time"
                        ].dt.year
                        >= start_year
                    )
                    & (
                        ds_gases_full_monthly_multi_phase_d[ghg]["scenarios-cmip6"][
                            "time"
                        ].dt.year
                        <= end_year
                    )
                )
            ),
            unstack_col="time",
            assign_metadata={"ghg": ghg, "cmip_era": "CMIP6"},
            ds_var=ghg,
        ).openscm.update_index_levels_from_other(
            {"experiment": ("scenario", lambda x: x)}
        ),
        to_data_frame(
            ds=ds_gases_full_monthly_multi_phase_d[ghg]["history-cmip6"].sel(
                time=(
                    (
                        ds_gases_full_monthly_multi_phase_d[ghg]["history-cmip6"][
                            "time"
                        ].dt.year
                        >= start_year
                    )
                    & (
                        ds_gases_full_monthly_multi_phase_d[ghg]["history-cmip6"][
                            "time"
                        ].dt.year
                        <= end_year
                    )
                )
            ),
            unstack_col="time",
            assign_metadata={
                "experiment": "historical",
                "ghg": ghg,
                "cmip_era": "CMIP6",
            },
            ds_var=ghg,
        ),
    ]
]
pdf = pd.concat(
    [v.reorder_levels(pdf_l[0].index.names) for v in pdf_l], axis="rows"
).sort_index(axis="columns")

pdf = pdf.openscm.update_index_levels_from_other(
    {"scenario_group": ("experiment", lambda x: scenario_group_map[x])}
)
pdf_grouped = {sg: sgdf for sg, sgdf in pdf.groupby("scenario_group")}
tmp_l = []
for sg, sgdf in pdf_grouped.items():
    if sg == "historical":
        continue

    tmp_l.append(sgdf)
    tmp_l.append(
        pdf_grouped["historical"].openscm.set_index_levels({"scenario_group": sg})
    )

pdf = pd.concat([v.reorder_levels(tmp_l[0].index.names) for v in tmp_l])
pdf.columns = [v.year + (v.month * 2 - 1) / 24 for v in pdf.columns]
pdf

# %% editable=true slideshow={"slide_type": ""} tags=["remove_input"]
start_year = 2010
end_year = 2030

sns.relplot(
    data=pdf.loc[
        :, (pdf.columns >= start_year) & (pdf.columns <= end_year + 1)
    ].openscm.to_long_data(),
    x="time",
    y="value",
    hue="experiment",
    palette=palette,
    style="cmip_era",
    # markers={"CMIP6": ".", "CMIP7": "o"},
    # edgecolors="none",
    markers={"CMIP6": "+", "CMIP7": 4},
    hue_order=hue_order_incl_cmip6,
    kind="scatter",
    row="ghg",
    col="scenario_group",
    col_order=["low", "continuing-trends", "high"],
    facet_kws=dict(sharey="row"),
    s=30,
    alpha=0.8,
)

plt.show()

# %% [markdown]
# # Bibliography
# ```{bibliography}
# :style: plain
# :filter: {"user-guide-scenarios"} & docnames
# ```
