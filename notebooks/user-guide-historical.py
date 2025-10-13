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
#   title: 'CMIP Greenhouse Gas (GHG) Concentration Historical Dataset:
#
#     Data Description and User Guide'
# ---

# %% [markdown] editable=true slideshow={"slide_type": ""}
# # Overview
#
# Here we provide a short description of the historical dataset
# and a guide for users.
# This is intended to provide a short introduction for users of the data.
# The full details of the dataset's construction
# and evaluation against other data sources
# will be provided in the full manuscript which is being prepared.

# %% [markdown] editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
# When ready, we will point to this manuscript here.
# [TODO cross-link]

# %% [markdown] editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
# ## Imports

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
import calendar
from functools import partial

import cftime
import matplotlib
import matplotlib.pyplot as plt
import nc_time_axis  # noqa: F401
import numpy as np
from local.data_loading import fetch_and_load_ghg_dataset, get_ghg_dataset_local_files
from local.esgf.db_helpers import create_all_tables, get_sqlite_engine
from local.esgf.search.search_query import KnownIndexNode
from local.paths import REPO_ROOT

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
local_data_root_dir = REPO_ROOT / "data" / "raw" / "esgf"
local_data_root_dir.mkdir(exist_ok=True, parents=True)
sqlite_file = REPO_ROOT / "download-test-database.db"
# # Obviously we wouldn't delete the database every time
# # in production, but while experimenting it's handy
# # to always start with a clean slate.
# if sqlite_file.exists():
#     sqlite_file.unlink()

engine = get_sqlite_engine(sqlite_file)
create_all_tables(engine)

# %% [markdown]
# # Finding and accessing the data

# %% [markdown]
# ## ESGF
#
# The Earth System Grid Federation (ESGF) provides access to a range of climate data.
#
# The historical data of interest here,
# which is the data to be used
# for historical and piControl simulations within CMIP [TODO ref Dunne paper],
# can be found under the "source ID", `CR-CMIP-1-0-0`.
# The concept of a source ID is a bit of a perculiar one
# to CMIP forcings data.
# In practice, it is simply a unique identifier for a collection of datasets
# (and it's best not to read more than that into it).
# It is possible to filter searches on ESGF
# via the user interface.
# Searches can often be encoded in URLs too
# (although these URLs sometimes move,
# so we make no guarantee that this link will always be live)
# e.g. [https://esgf-node.ornl.gov/search?project=input4MIPs&activeFacets=%7B%22source_id%22%3A%22CR-CMIP-1-0-0%22%7D]().
#
# To download the data, we recommend accessing it directly via the ESGF user interfaces
# via links like the one above.
# Alternately, there are tools dedicated to accessing ESGF data,
# with two prominent examples being:
#
# 1. esgpull: https://esgf.github.io/esgf-download
# 2. intake-esgf: https://intake-esgf.readthedocs.io
#
# Please see these tools' docs for usage instructions.

# %% [markdown]
# ## Zenodo
#
# While it aims to be, the ESGF is technically not a permanent archive
# and does not issue DOIs.
# In order to provide more reliable, citable access to the data,
# we also provide it on Zenodo.
# The data, as well as all the source code and input data used to process it,
# can be found at https://doi.org/10.5281/zenodo.14892947.

# %% [markdown] editable=true slideshow={"slide_type": ""}
# # Data description

# %% [markdown]
# ## Format
#
# The data is provided in netCDF format [TODO citation].
# This self-describing format allows the data
# to be placed in the same file as metadata
# (in the so-called "file header").
# To facilitate simpler use of the data,
# each dataset is split across multiple files.
# The advantage of this is that users do not need to load all years of data
# if they are only interested in data for a certain range,
# which can significantly improve data loading times.
# To get the complete dataset,
# the files can simply be concatenated in time.

# %% [markdown]
# ## Grids and frequencies provided
#
# We provide five combinations of grids and time sampling
# (also referred to as frequency,
# although this is a bit of a misuse as the units of frequency are per time,
# which doesn't match the convention for these metadata values).
# The grid and frequency information for each file can be found in its netCDF header
# under the attributes `grid_label` (for grid) and `frequency` (for time sampling).
# The `grid_label` and `frequency` also appear in each file's name,
# which allows files to be filtered without needing to load them first.
#
# The five combinations of grid and time sampling are:
#
# 1. global-, annual-mean (`grid_label="gm"`, `frequency="yr"`)
# 1. global-, monthly-mean (`grid_label="gm"`, `frequency="mon"`)
# 1. hemispheric-, annual-mean (`grid_label="gr1z"`, `frequency="yr"`)
# 1. hemispheric-, monthly-mean (`grid_label="gr1z"`, `frequency="mon"`)
# 1. 15-degree latitudinal, monthly-mean (`grid_label="gnz"`, `frequency="mon"`)

# %% [markdown]
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

# %% [markdown]
# ## Uncertainty
#
# At present, we provide no analysis of the uncertainty associated with these datasets.
# In radiative forcing terms, the uncertainty in these concentrations
# is very likely to be small compared to other uncertainties in the climate system,
# but this statement is not based on any robust analysis
# (rather it is based on expert judgement).
# It is also worth noting that the uncertainty increases as we go further back in time,
# particularly as we shift from using surface flasks to relying on ice cores instead.

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ## Differences compared to CMIP6
#
# At present, the changes from CMIP6 are minor,
# with the maximum difference in effective radiative forcing terms
# being 0.05 W / m<sup>2</sup>
# (and generally much smaller than this, particularly after 1850).
# For more details, see the plots in the user guide below
# and the forthcoming manuscript.

# %% [markdown] editable=true slideshow={"slide_type": ""}
# # User guide

# %% [markdown] editable=true slideshow={"slide_type": ""}
# Having downloaded the data, using it is quite straightforward.

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
# Ensure data is downloaded
query_kwargs_co2_yearly_global = dict(
    ghg="co2",
    time_sampling="yr",
    grid="gm",
    cmip_era="CMIP7",
    source_id="CR-CMIP-1-0-0",
    engine=engine,
)
fetch_and_load = partial(
    fetch_and_load_ghg_dataset,
    local_data_root_dir=local_data_root_dir,
    # index_node=KnownIndexNode.DKRZ,
    # cmip_era="CMIP6",
    # source_id="UoM-CMIP-1-2-0",
    index_node=KnownIndexNode.ORNL,
)
_ = fetch_and_load(**query_kwargs_co2_yearly_global)

# Get file paths
co2_yearly_global_fps = get_ghg_dataset_local_files(**query_kwargs_co2_yearly_global)

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ## Annual-, global-mean data
#
# We start with the annual-, global-mean data.
# Like all our datasets, this is composed of three files,
# each covering a different time period:
#
# 1. year 1 to year 999
# 2. year 1000 to year 1749
# 3. year 1750 to year 2022
#
# For yearly data, the time labels in the filename are years
# (for months, the month is included e.g. you will see `000101-09912`
# rather than `0001-0999` in the filename,
# the files also have different values for the `frequency` attribute).
# Global-mean data is identified by the 'grid label' `gm`,
# which appears in the filename.
# Below we show the filenames for the CO<sub>2</sub> output.

# %% editable=true slideshow={"slide_type": ""} tags=["remove_input"]
for fp in co2_yearly_global_fps:
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
for fp in ch4_yearly_global_fps:
    print(f"- {fp.name}")

# %% [markdown] editable=true slideshow={"slide_type": ""}
# As described above, the data is netCDF files.
# This means that metadata can be trivially inspected
# using a tool like `ncdump`.
# As you can see, there is a lot of metadata included in these files.
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
ds_co2_yearly_global = xr.open_mfdataset(co2_yearly_global_fps, decode_times=time_coder)

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
# Force values to compute to avoid dask getting involved
ds_co2_yearly_global = ds_co2_yearly_global.compute()

# %% editable=true slideshow={"slide_type": ""}
ds_co2_yearly_global

# %% editable=true slideshow={"slide_type": ""} tags=["remove_input"]
ds_co2_yearly_global["co2"].plot.scatter(alpha=0.4)
plt.show()

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ## Space- and time-average nature of the data
#
# All of our data represents the mean over each cell.
# This is indicated by the `cell_methods` attribute
# of all of our output variables.

# %% editable=true slideshow={"slide_type": ""}
ds_co2_yearly_global["co2"].attrs["cell_methods"]

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
ds_co2_yearly_global["time_bnds"]

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
ds_plt = ds_co2_yearly_global.isel(time=slice(-5, None))

fig, ax = plt.subplots(figsize=(8, 4))
ds_plt["co2"].plot.scatter(ax=ax)

for bounds, val in zip(ds_plt["time_bnds"].values, ds_plt["co2"].values):
    ax.plot(bounds, [val, val], color="tab:blue", linewidth=1.0, alpha=0.7)

xticks = [cftime.DatetimeGregorian(y, 1, 1) for y in range(2018, 2024)]
ax.set_xticks(xticks)
ax.set_xlim(xticks[0], xticks[-1])
ax.grid()

plt.show()

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ## Monthly-, global-mean data
#
# If you want to have information at a finer level
# of temporal detail, we also provide monthly files.
# Like the global datasets, these come in three files.
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
for fp in co2_monthly_global_fps:
    print(f"- {fp.name}")

# %% [markdown] editable=true slideshow={"slide_type": ""}
# Again, the data can be trivially loaded with [xarray](https://github.com/pydata/xarray).

# %% editable=true slideshow={"slide_type": ""}
ds_co2_monthly_global = xr.open_mfdataset(
    co2_monthly_global_fps, decode_times=time_coder
)

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
# Force values to compute to avoid dask getting involved
ds_co2_monthly_global = ds_co2_monthly_global.compute()

# %% editable=true slideshow={"slide_type": ""}
ds_co2_monthly_global

# %% [markdown] editable=true slideshow={"slide_type": ""}
# For this data, the time bounds show that each point
# is the average a month, not a year.

# %% editable=true slideshow={"slide_type": ""}
ds_co2_monthly_global["time_bnds"]

# %% [markdown]
# As above, as a result of the time average that the data represents,
# it is inappropriate to plot this data using a line plot.
# Scatter or step plots should be used instead.

# %% editable=true slideshow={"slide_type": ""} tags=["remove_input"]
ds_plt = ds_co2_monthly_global.isel(time=slice(-5 * 12, None))

fig, ax = plt.subplots(figsize=(8, 4))
ds_plt["co2"].plot.scatter(ax=ax)

for bounds, val in zip(ds_plt["time_bnds"].values, ds_plt["co2"].values):
    ax.plot(bounds, [val, val], color="tab:blue", linewidth=1.0, alpha=0.7)

xticks = [cftime.DatetimeGregorian(y, 1, 1) for y in range(2018, 2024)]
ax.set_xticks(xticks)
ax.set_xlim(xticks[0], xticks[-1])
ax.grid()

plt.show()

# %% [markdown]
# The monthly data includes seasonality.
# Plotting the monthly and yearly data
# on the same axes makes particularly clear
# why a line plot is inappropriate.

# %% editable=true slideshow={"slide_type": ""} tags=["remove_input"]
fig, ax = plt.subplots(figsize=(8, 4))

for ds_plt, label, colour in (
    (ds_co2_monthly_global.isel(time=slice(-5 * 12, None)), "monthly", "tab:blue"),
    (ds_co2_yearly_global.isel(time=slice(-5, None)), "yearly", "tab:orange"),
):
    ds_plt["co2"].plot.scatter(ax=ax, label=label, color=colour, s=10)

    for bounds, val in zip(ds_plt["time_bnds"].values, ds_plt["co2"].values):
        ax.plot(bounds, [val, val], color=colour, linewidth=1.0, alpha=0.7)

ax.legend()

xticks = [cftime.DatetimeGregorian(y, 1, 1) for y in range(2018, 2024)]
ax.set_xticks(xticks)
ax.set_xlim(xticks[0], xticks[-1])
ax.grid()

plt.show()

# %% [markdown]
# ## Monthly-, latitudinally-resolved data
#
# We also provide data with spatial,
# specifically latituindal, resolution.
# This data comes on a 15-degree latituindal grid
# (see below for details of the grid and latituindal bounds).
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
for fp in co2_monthly_lat_fps:
    print(f"- {fp.name}")

# %% [markdown] editable=true slideshow={"slide_type": ""}
# Again, the data can be trivially loaded with [xarray](https://github.com/pydata/xarray).

# %% editable=true slideshow={"slide_type": ""}
ds_co2_monthly_lat = xr.open_mfdataset(
    co2_monthly_lat_fps, decode_times=time_coder, data_vars=None, compat="no_conflicts"
)

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
# Force values to compute to avoid dask getting involved
ds_co2_monthly_lat = ds_co2_monthly_lat.compute()

# %% editable=true slideshow={"slide_type": ""}
ds_co2_monthly_lat

# %% [markdown] editable=true slideshow={"slide_type": ""}
# For this data, the latitudinal bounds show the area
# over which each point is the average.

# %% editable=true slideshow={"slide_type": ""}
ds_co2_monthly_lat["lat_bnds"]

# %% [markdown]
# As above, but this time for the spatial axis,
# it is inappropriate to plot this data using a line plot.
# Scatter or step plots should be used instead.

# %% editable=true slideshow={"slide_type": ""} tags=["remove_input"]
ds_plt = ds_co2_monthly_lat.isel(time=slice(-12, None))


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
    axes_d[f"2022 - {calendar.month_name[month]}"].set_ylabel(
        "Latitude (degrees north)"
    )

for month in range(10, 13):
    axes_d[f"2022 - {calendar.month_name[month]}"].set_xlabel("co2 [ppm]")
# ax.legend(loc="center left", bbox_to_anchor=(1.05, 0.5))

plt.tight_layout()
plt.show()

# %% [markdown] editable=true slideshow={"slide_type": ""}
# We can compare the global-mean data
# to the data at each latitude.
# The strength of the latitudinal gradient varies also by gas (not shown).

# %% editable=true slideshow={"slide_type": ""} tags=["remove_input"]
fig, ax = plt.subplots(figsize=(8, 4))

time_slice = slice(-5 * 12, None)

ds_plt = ds_co2_monthly_global.isel(time=time_slice)
ds_plt["co2"].plot.scatter(
    ax=ax, label="global-mean", color="tab:blue", s=30, zorder=10.0
)

for bounds, val in zip(ds_plt["time_bnds"].values, ds_plt["co2"].values):
    ax.plot(bounds, [val, val], color="tab:blue", linewidth=1.0, alpha=0.7)

ds_all_lats = ds_co2_monthly_lat.isel(time=time_slice)

for i, lat in enumerate(sorted(ds_all_lats["lat"])[::-1]):
    ds_plt = ds_all_lats.sel(lat=lat)
    colour = matplotlib.colormaps["magma"](i / len(ds_all_lats["lat"]))

    ds_plt["co2"].plot.scatter(
        ax=ax, label=f"{float(lat)}", color=colour, marker="x", s=10
    )

    for bounds, val in zip(ds_plt["time_bnds"].values, ds_plt["co2"].values):
        ax.plot(bounds, [val, val], color=colour, linewidth=1.0, alpha=0.7)

ax.legend(loc="center left", bbox_to_anchor=(1.05, 0.5))

xticks = [cftime.DatetimeGregorian(y, 1, 1) for y in range(2018, 2024)]
ax.set_xticks(xticks)
ax.set_xlim(xticks[0], xticks[-1])
ax.grid()
ax.set_title("")

plt.show()

# %% [markdown]
# The data can also be plotted in a so-called "magic carpet"
# to see the variation in space and time simultaneously.

# %% editable=true slideshow={"slide_type": ""} tags=["remove_input"]
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(projection="3d")

tmp = ds_co2_monthly_lat["co2"].isel(time=range(-10 * 12, 0)).copy()
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
# ## Differences from CMIP6
#
# ### File formats and naming
#
# The file formats are generally close to CMIP6.
# There are three key changes:
#
# 1. we have split the global-mean and hemispheric-mean data into separate files.
#    In CMIP6, this data was in the same file (with a grid label of `GMNHSH`).
#    We have split this for two reasons:
#    a) `GMNHSH` is not a grid label recognised in the CMIP CVs and
#    b) having global-mean and hemispheric-mean data in the same file
#       required us to introduce a 'sector' coordinate,
#       which was confusing and does not follow the CF-conventions.
# 1. we have split the files into different time components.
#    One file goes from year 1 to year 999 (inclusive).
#    The next file goes from year 1000 to year 1749 (inclusive).
#    The last file goes from year 1750 to year 2022 (inclusive).
#    This simplifies handling and allows groups to avoid loading data
#    they are not interested in (for CMIP, this generally means data pre-1750).
# 1. we have simplified the names of all the variables.
#    They are now simply the names of the gases,
#    for example we now use "co2" rather than "mole_fraction_of_carbon_dioxide".
#    A full mapping is provided below.
#
# There is one more minor change.
# The data now starts in year one, rather than year zero.
# We do this because year zero doesn't exist in most calendars
# (and we want to avoid users of the data having to hack around this
# when using standard data analysis tools).
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
# Comparing the data from CMIP6 and CMIP7 shows minor changes
# (although doing this comparison requires a bit of care
# because of the changes in file formats).

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
gases_to_show = ["co2", "ch4", "n2o", "cfc12eq", "hfc134aeq"]
ds_gases_full_d = {}
for gas in gases_to_show:
    ds_gases_full_d[gas] = {}
    for source_id, cmip_era in (
        ("CR-CMIP-1-0-0", "CMIP7"),
        ("UoM-CMIP-1-2-0", "CMIP6"),
    ):
        query_kwargs = {
            "ghg": gas,
            "time_sampling": "yr",
            "grid": "gm",
            "source_id": source_id,
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
        ds_gases_full_d[gas][cmip_era] = ds.compute()

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
from typing import Callable

import numpy.typing as npt


def sel_times(
    ds_d: dict[str, dict[str, xr.Dataset]],
    sel_func: Callable[[xr.DataArray], npt.NDArray[bool]],
) -> dict[str, dict[str, xr.Dataset]]:
    """
    Select times from our dictionary of [xr.Dataset][]'s
    """
    res = {
        gas: {
            key: value.sel(time=sel_func(value["time"])) for key, value in tmp.items()
        }
        for gas, tmp in ds_d.items()
    }

    return res


def plot_overview_and_deltas(
    ds_d: dict[str, dict[str, xr.Dataset]],
    axes_d: dict[str, matplotlib.axes.Axes],
):
    """
    Plot overviews of timeseries and deltas between CMIP7 and CMIP6
    """
    for ax_name, ax in axes_d.items():
        if ax_name.endswith("_delta"):
            continue

        gas = ax_name

        for cmip_era, ds in ds_d[gas].items():
            label = f"{cmip_era} ({ds.attrs['source_id']})"
            ds[gas].plot.scatter(
                ax=axes_d[gas], label=label, alpha=0.7, edgecolors="none"
            )

        ax.legend()
        ax.set_title(gas)
        ax.xaxis.set_tick_params(labelbottom=True)

        ax_delta = axes_d[f"{gas}_delta"]

        da_cmip6 = ds_d[gas]["CMIP6"][gas]
        da_cmip7 = ds_d[gas]["CMIP7"][gas]
        overlapping_times = np.intersect1d(da_cmip6["time"], da_cmip7["time"])
        delta = da_cmip7.sel(time=overlapping_times) - da_cmip6.sel(
            time=overlapping_times
        )
        ax_delta.set_title("CMIP7 - CMIP6", fontsize="small")
        delta.plot.scatter(
            ax=ax_delta,
            color="tab:grey",
            edgecolors="none",
            s=10,
        )
        ax_delta.axhline(0.0, color="k", linestyle="--")

        ax_delta.xaxis.set_tick_params(labelbottom=True)


plt_mosaic = [
    ["co2", "ch4", "n2o"],
    ["co2", "ch4", "n2o"],
    ["co2_delta", "ch4_delta", "n2o_delta"],
    ["cfc12eq", "hfc134aeq", ""],
    ["cfc12eq", "hfc134aeq", ""],
    ["cfc12eq_delta", "hfc134aeq_delta", ""],
]
get_default_delta_mosaic = partial(
    plt.subplot_mosaic,
    mosaic=plt_mosaic,
    figsize=(12, 8),
    sharex=True,
)


def remove_empty_axes(
    axes_d: dict[str, matplotlib.axes.Axes],
) -> dict[str, matplotlib.axes.Axes]:
    """Remove empty axes"""
    res = {}
    for k, v in axes_d.items():
        if k:
            res[k] = v

        else:
            v.remove()

    return res


# %% [markdown] editable=true slideshow={"slide_type": ""}
# #### Atmospheric concentrations: Year 1 - 2022

# %% editable=true slideshow={"slide_type": ""} tags=["remove_input"]
fig, axes_d = get_default_delta_mosaic()
axes_d = remove_empty_axes(axes_d)

plot_overview_and_deltas(
    ds_gases_full_d,
    axes_d,
)

plt.tight_layout()
plt.show()

# %% [markdown] editable=true slideshow={"slide_type": ""}
# #### Atmospheric concentrations: Year 1750 - 2022

# %% editable=true slideshow={"slide_type": ""} tags=["remove_input"]
fig, axes_d = get_default_delta_mosaic()
axes_d = remove_empty_axes(axes_d)

min_year = 1750
plot_overview_and_deltas(
    sel_times(ds_gases_full_d, lambda x: x.dt.year >= min_year),
    axes_d,
)

plt.tight_layout()
plt.show()

# %% [markdown] editable=true slideshow={"slide_type": ""}
# #### Atmospheric concentrations: Year 1957 - 2022
#
# 1957 is the start of the Scripps ground-based record.
# Before this, data is based on ice cores alone.

# %% editable=true slideshow={"slide_type": ""} tags=["remove_input"]
fig, axes_d = get_default_delta_mosaic()
axes_d = remove_empty_axes(axes_d)

min_year = 1957
plot_overview_and_deltas(
    sel_times(ds_gases_full_d, lambda x: x.dt.year >= min_year),
    axes_d,
)

plt.tight_layout()
plt.show()

# %% [markdown] editable=true slideshow={"slide_type": ""}
# #### Approximate radiative effect: Year 1 - 2022
#
# As seen above, in atmospheric concentration terms
# the differences are small.
# However, this can be put on a common scale
# by comparing the differences in radiative effect terms.
# This gives an approximation of the size of the difference
# that would be seen by an Earth System Model's (ESM's) radiation code.
# This uses basic linear approximations,
# assuming that the radiative effect of each gas
# is simply its atmospheric concentration multiplied by a constant.
# This isn't the same as effective radiative forcing (ERF).
# For that comparison, see the later sections focussed on ERF.

# %% [markdown] editable=true slideshow={"slide_type": ""}
# Values below come from Table 7.SM.7 of
# IPCC AR7 WG1 Ch. 7 Supplementary Material
# (https://www.ipcc.ch/report/ar6/wg1/downloads/report/IPCC_AR6_WGI_Chapter07_SM.pdf).

# %% editable=true slideshow={"slide_type": ""}
from openscm_units import unit_registry

Q = unit_registry.Quantity

RADIATIVE_EFFICIENCIES = {
    "co2": Q(1.33e-5, "W / m^2 / ppb"),
    "ch4": Q(3.88e-4, "W / m^2 / ppb"),
    "n2o": Q(3.2e-3, "W / m^2 / ppb"),
    "cfc12eq": Q(0.358, "W / m^2 / ppb"),
    "hfc134aeq": Q(0.167, "W / m^2 / ppb"),
}

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
ds_gases_full_radiative_effect_d = {}
target_units = "W / m^2"
for gas, gas_ds in ds_gases_full_d.items():
    ds_gases_full_radiative_effect_d[gas] = {}
    for mip_era, ds in gas_ds.items():
        tmp = ds.copy()

        tmp[gas][:] = (
            (Q(tmp[gas].values, tmp[gas].attrs["units"]) * RADIATIVE_EFFICIENCIES[gas])
            .to(target_units)
            .m
        )
        tmp[gas].attrs["units"] = target_units
        tmp[gas].attrs["long_name"] = "approx. radiative effect"

        ds_gases_full_radiative_effect_d[gas][mip_era] = tmp

# %% editable=true slideshow={"slide_type": ""} tags=["remove_input"]
fig, axes_d = get_default_delta_mosaic()
axes_d = remove_empty_axes(axes_d)

plot_overview_and_deltas(
    ds_gases_full_radiative_effect_d,
    axes_d,
)

for name, ax in axes_d.items():
    if name.endswith("_delta"):
        continue

    ax.set_ylim([0, 6.0])

plt.tight_layout()
plt.show()

# %% [markdown]
# #### Approximate radiative effect: Year 1750 - 2022
#
# This is the period relevant for historical simulations in CMIP.

# %% editable=true slideshow={"slide_type": ""} tags=["remove_input"]
fig, axes_d = get_default_delta_mosaic()
axes_d = remove_empty_axes(axes_d)

min_year = 1750
plot_overview_and_deltas(
    sel_times(ds_gases_full_radiative_effect_d, lambda x: x.dt.year >= min_year),
    axes_d,
)

for name, ax in axes_d.items():
    if name.endswith("_delta"):
        continue

    ax.set_ylim([0, 6.0])

plt.tight_layout()
plt.show()

# %% [markdown] editable=true slideshow={"slide_type": ""}
# #### Approximate effective radiative forcing: Year 1750 - 2022
#
# The above isn't effective radiative forcing.
# For that, you have to normalise the data to some reference year.
# There are a few different choices for this reference year.
# In IPCC reports, it is 1750 so that is what we show here.
# It should be noted that some ESMs may make other choices,
# but these would not have a great effect on the interpretation
# of the difference between the CMIP6 and CMIP7 datasets.
#
# Note that this approximation is linear,
# which is a particularly strong approximation for CO<sub>2</sub>
# because of its logarithmic forcing nature.
# We show this approximation here nonetheless
# because it provides an order of magnitude estimate
# for the change from CMIP6 in ERF terms.
# The forthcoming manuscripts will explore the subtleties
# of this quantification in more detail.

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
ds_gases_full_erf_d = {}
reference_year = 1750
for gas, gas_ds in ds_gases_full_radiative_effect_d.items():
    ds_gases_full_erf_d[gas] = {}
    for mip_era, ds in gas_ds.items():
        tmp = ds.copy()

        tmp[gas][:] = (
            tmp[gas][:]
            - tmp.sel(time=ds["time"].dt.year == reference_year)[gas][:].values
        )
        tmp[gas].attrs["long_name"] = "approx. ERF"
        # tmp.attrs = ds.attrs

        ds_gases_full_erf_d[gas][mip_era] = tmp

# %% editable=true slideshow={"slide_type": ""} tags=["remove_input"]
fig, axes_d = get_default_delta_mosaic()
axes_d = remove_empty_axes(axes_d)

min_year = 1750
plot_overview_and_deltas(
    sel_times(ds_gases_full_erf_d, lambda x: x.dt.year >= min_year),
    axes_d,
)

for name, ax in axes_d.items():
    if name.endswith("_delta"):
        continue

    ax.set_ylim([0, 2.0])

plt.tight_layout()
plt.show()

# %% [markdown]
# In summary, in ERF terms, the differences from CMIP6 are very small.
# For all gases, they are less than around 0.025 W / m<sup>2</sup>.
# Compared to the estimated total greenhouse gas forcing and uncertainty in IPCC AR6
# (see Section 7.3.5.2 of AR6 WG1 Chapter 7,
# https://www.ipcc.ch/report/ar6/wg1/chapter/chapter-7/),
# estimated to be 3.84 W / m<sup>2</sup>
# (very likely range of 3.46 to 4.22 W / m<sup>2</sup>),
# such differences are particularly small.

# %% [markdown]
# #### Atmospheric concentrations including seasonality: Year 2000 - 2022
#
# The final comparisons we show are atmospheric concentrations including seasonality.
# Given that most greenhouse gases
# are well-mixed with lifetimes much greater than a year,
# these differences are unlikely to be of huge interest to ESMs.
# However, for other applications, such seasonality differences may matter more.

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
ds_gases_full_monthly_d = {}
for gas in gases_to_show:
    ds_gases_full_monthly_d[gas] = {}
    for source_id, cmip_era in (
        ("CR-CMIP-1-0-0", "CMIP7"),
        ("UoM-CMIP-1-2-0", "CMIP6"),
    ):
        query_kwargs = {
            "ghg": gas,
            "time_sampling": "mon",
            "grid": "gm",
            "source_id": source_id,
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
        ds_gases_full_monthly_d[gas][cmip_era] = ds.compute()

# %% editable=true slideshow={"slide_type": ""} tags=["remove_input"]
fig, axes_d = get_default_delta_mosaic()
axes_d = remove_empty_axes(axes_d)

min_year = 2000
max_year = 2022
plot_overview_and_deltas(
    sel_times(
        ds_gases_full_monthly_d,
        lambda x: (x.dt.year >= min_year) & (x.dt.year <= max_year),
    ),
    axes_d,
)

plt.tight_layout()
plt.show()

# %% [markdown] editable=true slideshow={"slide_type": ""}
# Like the annual-means,
# the atmospheric concentrations including seasonality
# are reasonably consistent between CMIP6 and CMIP7.
# There are some areas of change.
# Full details of these changes will be provided
# in the forthcoming manuscripts.
