# ---
# jupyter:
#   authors:
#   - name: Anna Lanteri
#   - name: Zebedee Nicholls
#   - name: Florence Bockting
#   - name: Mika Pflüger
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
#   title: 'Adding EO data to the generation of CMIP Greenhouse Gas (GHG) Concentration
#     Historical Dataset:
#
#     Data Description and User Guide'
# ---

# %% [markdown]
# ```{role} raw-latex(raw)
# :format: latex
# ```

# %% [markdown] editable=true slideshow={"slide_type": ""}
# # Overview
#
# Here we provide a short description of the historical dataset of concentrations for CO{raw-latex}`\textsubscript{2}` and CH{raw-latex}`\textsubscript{4}` updated with satellite data
# and a guide for users.
# This is intended to provide a short introduction for users of the data:
# its construction, key features, metadata
# and relationship to CMIP6 and CMIP7 forcing data.
#
# The dataset is an extension on the forcings for CMIP7 datasets for CO{raw-latex}`\textsubscript{2}` and CH{raw-latex}`\textsubscript{4}`, achieved by pre-processing and including satellite data to the GHG concentration generation pipeline. In this document, we will describe both the original dataset and the updated verssion, to have a complete and self-contained description. This means this document contains duplicated information from the user guide: 'CMIP7 Greenhouse Gas (GHG)
# Concentration Forcing Historical
# Dataset'.

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
from myst_nb import glue

from local.data_loading import fetch_and_load_ghg_dataset
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

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
# TODO: update this to point at the final, published EO-update dataset
# once it exists (e.g. once it's on ESGF and can be fetched the same way
# as the original dataset above). For now, this points directly at a local
# dev-run output bundle from CMIP-GHG-Concentration-Generation, so there is
# no ESGF fetch step for this data.
EO_UPDATE_DATA_ROOT = (
    REPO_ROOT.parent
    / "CMIP-GHG-Concentration-Generation"
    / "output-bundles"
    / "dev-test-run"
    / "data"
    / "processed"
    / "esgf-ready"
    / "input4MIPs"
    / "CMIP6Plus"
    / "CMIP"
    / "CR"
    / "CR-CMIP-testing"
    / "atmos"
)
EO_UPDATE_DATA_VERSION = "v20260907"
# The two gases use different statistical fits for the satellite-extension step
EO_UPDATE_FIT_SUFFIX = {
    "co2": "SAT_LINEAR_SEASONAL_LAT_STD_WEIGHT_FIT",
    "ch4": "SAT_NONLINEAR_LAT_STD_WEIGHT_FIT",
}


def get_eo_update_local_files(ghg, time_sampling, grid):
    """Get the local files for the EO-updated (satellite-extended) dataset"""
    d = EO_UPDATE_DATA_ROOT / time_sampling / ghg / grid / EO_UPDATE_DATA_VERSION
    return sorted(d.glob(f"*_{EO_UPDATE_FIT_SUFFIX[ghg]}.nc"))


# %% [markdown] jp-MarkdownHeadingCollapsed=true
# # Dataset construction
#
# The original CMIP7 GHG concentration dataset (produced without using satellite data) is constructed following the methodology of
# {raw-latex}`\textcite{meinshausen_historical_2017}`.
# The methods are described in full in that paper
# and will be clarified and described again
# in the forthcoming manuscript describing this dataset's construction.
#
# In brief, the dataset for each greenhouse gas is constructed via the following steps
# (for full details and code, see GitHub repository
# [CMIP-GHG-Concentration-Generation](https://github.com/climate-resource/CMIP-GHG-Concentration-Generation)
# )[^gh-code]:
#
# 1. collect as many ground-based observations as possible from ground-based networks such as the NOAA
#    {raw-latex}`\parencite{lan_atmospheric_co2_2025,lan_atmospheric_ch4_2025}`
#    and AGAGE
#    {raw-latex}`\parencite{prinn_history_2000,prinn2018history,rigby2008renewed,rigby2017role}`
#    networks
#    (the full set of input sources are documented in the `references*`
#    global attributes of the output files
#    and will be discussed in more detail in a forthcoming paper)
#     - these are only available over the last few decades at most
#       (less for some greenhouse gases)
#     - these are spatially sparse because sampling stations
#       are discrete points and there are not an infinite number of stations
#       (at most, usually around 30, often far fewer)
# 4. bin the ground-based observations in space and time
#    {raw-latex}`\parencite[15-degree latitudinal bins, 60-degree longitudinal bins, monthly time bins,
#    following][]{meinshausen_historical_2017}`,
#    averaging over input stations and observations that fall in the same cell
# 5. interpolate the binned data in space using a standard 2D linear interpolation
#    as in {raw-latex}`\textcite{meinshausen_historical_2017}`,
#    to derive a dataset with spatial coverage
# 6. use the interpolated, ground-based data
#    to derive a statistical model for seasonal variation and latitudinal gradients
#    specific to each greenhouse gas
#      - the exact form of the statistical model varies by gas
#        (we use empirical orthogonal function methods,
#        linear regression and basic scaling arguments)
#        but is generally driven by either concentrations of the gas itself,
#        global-mean temperature or purely statistical regressions/extensions
# 7. use the models, plus ice core or other proxy records,
#    to extend global-mean concentrations, seasonality and latitudinal gradients
#    over the full time period of the dataset (i.e. back to year 1)
#      - where ice cores or proxy records are not available,
#        purely statistical extrapolations are used instead
#      - the extension varies by gas,
#        (we use optimisation, linear regression and spline interpolation),
#        aiming to make use of as much information as is possible
#        e.g. hemisphere specific ice core information
#        and the latitudinal gradient
#        over the period covered by ground-based observations
# 8. combine the extended global-mean, seasonality and latitudinal gradients
#    to create a dataset that extends over the period
#    year 1 to 2022 (the last year available for some observational networks
#    at the time the data was compiled)
#      - this dataset is on our binned grid,
#        which we choose to be a grid comprised of latitudinal bins 15-degrees in size
#      - it is not trivial to infer the global-means,
#        seasonality and latitudinal gradient used to construct the dataset
#        from the output dataset. For this reason,
#        we include these components separately
#        in the zenodo record[^zenodo-record]
#        that archives the output dataset,
#        all its inputs and intermediate data prdoucts
# 9. calculate annual-, hemispheric- and global-means
#    to produce our lower resolution data products
#      - we can also produce higher spatial resolution data products,
#        but have not done so at the moment to save processing and storage space
#        given that there has been no demand for these products from modelling teams
#
# The updated dataset (produced by adding satellite data information, hereafter referred to as EO-CMIP7) is produced with the same methodology, save from the fact that pre-processed satellite data is also used ats input in the first step. The pre-processing consists in scaling the satellite data using a function obtained by fitting the satellite data to the ground based data. The assumption is that this simple approach already provides us with a new dataset, therefore referred as 'scaled satellite data', compatible with ground based data, but with a wider coverage. Evalutation of results will focus on determining the fairness of this assumption.
#
# The input datasets, including the satellite datasets, and associated references
# are documented in the `references*` attributes of each netCDF file.
# This documentation is limited, so cannot document how each input dataset is used
# (that is the role of the manuscript),
# but does provide machine-readable provenance information
# (which is used to support links between all the input data
# e.g. linking of the Zenodo archive underpinning this dataset).
#
# [^zenodo-record]: https://doi.org/10.5281/zenodo.14892947
# [^gh-code]: https://github.com/climate-resource/CMIP-GHG-Concentration-Generation
#
#
# ## Satellide data preprocessing
#
# The satellite data used in this work is taken from OBS4MIPs L3 total column (XCH4/XCO2) product. This consists in a harmonized gridded multi-satellite merged product using EMMA, calibrated using the TCCON network, using an a-priori profile generated with SLIM.
#
# [TODO: update references]
#
# The satellite data is then pre-processed with the objective of producing a new dataset with the coverage of the satellite product but with values that can be treated in the same way as data from a ground based network. To achieve this, we scale the satellite data with a factor/a function derived by the relationship between satellite and ground based.
#
# We therefore perform a series of fits between matching (in space and time) datapoints for satellite and ground-based data. These fits include dependencies on concentrations, latitude, and seasonality, since those variables are identified as the key components for the PC decomposition performed in the GHG concentrations generation pipeline. Satellite data is then scaled using the function given by the best fit for each gas, which is respectively:
#
# - CO{raw-latex}`\textsubscript{2}`: linear dependency on concentration, latitude, seasonality, using inverse variance weighting
# - CH{raw-latex}`\textsubscript{4}`: non-linear (quadratic) dependency on concentration, linear dependency on latitude, using inverce variance weighting.
#
# An uncertainty estimation is then performed combining a Monte Carlo approach and fit uncertainty. We first run the fit 200 times perturbing the satellite and ground data within their error bars, taking then the variance of the distribution for each gridpoint and time. We then multiply all resulting uncertainty values by the squared root of the reduced {raw-latex}`\ensuremath{\chi^2}` (aka the standard error of the regression) wherever the reduced {raw-latex}`\ensuremath{\chi^2}` is > 1.
#
# The fitting and therefore pre-prossing of the data is kept simple for interpretability and to allow straightforward uncertainty estimation, but more sophisticated approaches including machine learning algorithms or complex generative approaches remain of interest for future works.
#
# ## Adding the satellite to the GHG generation pipeline
#
# To reduce the impact of satellite data on bins that already have robust information from the ground-based networks, we also weight the satellite data using the uncertainty we described above. We do not apply any weighting to the ground-based data, to remain as close as possible to the current pipeline. This results in the scaled satellite data being treated as an additional ground-based network where no other data is available, but being de facto ignored wherever we already have strong coverage.
#

# %% [markdown]
# # Finding and accessing the data

# %% [markdown] jp-MarkdownHeadingCollapsed=true
# ## ESGF
#
# [TODO: clarify if theese data will be stored on esgf, I assume not?]
#
# The **Earth System Grid Federation** {raw-latex}`\parencite{esgf_docs}`
# provides access to a range of climate data.
# The historical data of interest here,
# which is the data to be used
# for historical and piControl simulations within CMIP
# {raw-latex}`\parencite{dunne2025evolving}`,
# can be found under the "source ID", `CR-CMIP-1-0-0`.
# The concept of a "source ID" is a bit of a perculiar one
# to CMIP forcings data.
# In practice, it is simply a unique identifier for a collection of datasets
# (and it's best not to read more than that into it).
#
# It is possible to filter searches on ESGF
# via the user interface (see ESGF user guides[^esgf-user-guides-url]).
# Alternatively, searches can be encoded in URLs. However, a caveat with this
# approach is that URLs sometimes move, so we make no guarantee that this link
# will always be live. The following link provides an example of a search
# that is encoded in a URL:
#
# > [esgf-node.ornl.gov/search?project=input4MIPs&activeFacets=%7B%2ource_id%22%3A%22CR-CMIP-1-0-0%22%7D](https://esgf-node.ornl.gov/search?project=input4MIPs&activeFacets=%7B%2ource_id%22%3A%22CR-CMIP-1-0-0%22%7D)
#
# To download the data, we recommend accessing it directly via the ESGF user interfaces
# via links like the one above.
# Alternately, there are tools dedicated to accessing ESGF data,
# with two prominent examples being **esgpull**[^esgpull-url]
# and **intake-esgf**[^intake-esgf-url].
# Please refer to the tools' docs for usage instructions.
#
# [^esgf-user-guides-url]: https://esgf.github.io/esgf-user-support/user_guide.html#data-search-and-download
# [^esgpull-url]: https://esgf.github.io/esgf-download
# [^intake-esgf-url]: https://intake-esgf.readthedocs.io

# %% [markdown]
# ## Zenodo
#
# While it aims to be, the ESGF is technically not a permanent archive
# and does not issue DOIs.
# In order to provide more reliable, citable access to the data,
# we also provide it on **Zenodo** {raw-latex}`\parencite{zenodo}`.
# The data, as well as all the source code and input data used to process it,
# can be found at [TODO: add data on zenodo, put link here]

# %% [markdown] editable=true slideshow={"slide_type": ""}
# # Data description

# %% [markdown]
# ## Format
#
# The data is provided in **netCDF format** {raw-latex}`\parencite{zenodo}`.
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

# %% [markdown]
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
# The species provided are:
#
# <!-- Note: scripts/generate-ghg-listing.py cannot produce this subset (see below), so this list is maintained by hand --->
# - CH{raw-latex}`\textsubscript{4}`, in ppb
# - CO{raw-latex}`\textsubscript{2}`, in ppm
#
# This means that, compared with the original CMIP7 forcing dataset, we do not provide data for N{raw-latex}`\textsubscript{2}`O, ozone depleting substances, HCFCs, Halons, nor ozone flourinated compounds.
#
# Adding satellite data to this pipeline for other gases is non-trivial and heavily depends on the ground-based data quality and availability, physics and chemistry of the total column of gas and their impacts on the satellite retrieval algorithms, and on the available satellite data offering.

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
# ## Differences compared to CMIP6 and original
#
# [TODO: update with correct analysis]
# At present, the changes from CMIP6 are minor,
# with the maximum difference in effective radiative forcing terms
# being 0.05 W / m{raw-latex}`\textsuperscript{2}`
# (and generally much smaller than this, particularly after 1850).
# For more details, see the plots in the user guide below
# and the forthcoming manuscript.

# %% [markdown] editable=true slideshow={"slide_type": ""}
# # User guide

# %% [markdown] editable=true slideshow={"slide_type": ""}
# Having downloaded the data, using it is quite straightforward.

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
fetch_and_load = partial(
    fetch_and_load_ghg_dataset,
    local_data_root_dir=local_data_root_dir,
    # index_node=KnownIndexNode.DKRZ,
    # cmip_era="CMIP6",
    # source_id="UoM-CMIP-1-2-0",
    index_node=KnownIndexNode.ORNL,
)

# Get file paths for the EO-updated data
# (read directly from the local dev-run bundle, see EO_UPDATE_DATA_ROOT above;
# there is no ESGF fetch for this data yet)
co2_yearly_global_fps = get_eo_update_local_files(
    ghg="co2", time_sampling="yr", grid="gm"
)

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ## Annual-, global-mean data
#
# [TODO: format of the files needs to be finalised, so the metadata shown here will most likely change]
#
# We start with the annual-, global-mean data.
# Like all our datasets, this is composed of three files,
# each covering a different time period:
#
# 1. year 1 to year 999
# 2. year 1000 to year 1749
# 3. year 1750 to year 2022

# %% [markdown]
# For yearly data, the time labels in the filename are years
# (for months, the month is included e.g. you will see `000101-09912`
# rather than `0001-0999` in the filename,
# the files also have different values for the `frequency` attribute).
# Global-mean data is identified by the 'grid label' `gm`,
# which appears in the filename.
# Below we show the filenames for the CO{raw-latex}`\textsubscript{2}` output.

# %% editable=true slideshow={"slide_type": ""} tags=["remove_input"]
for fp in co2_yearly_global_fps:
    print(f"- {fp.name}")

# %% [markdown] editable=true slideshow={"slide_type": ""}
# For methane, similarly, the filenames are:

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
# Get file paths for the EO-updated data (see EO_UPDATE_DATA_ROOT above)
ch4_yearly_global_fps = get_eo_update_local_files(
    ghg="ch4", time_sampling="yr", grid="gm"
)

# %% editable=true slideshow={"slide_type": ""} tags=["remove_input"]
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
# Using a tool like [xarray](https://github.com/pydata/xarray)[^xarray-url],
# loading and working with the data is trivial.
# The resulting plot is shown in {numref}`Figure %s <ds-co2-yearly-global-fig>`.
#
# [^xarray-url]: https://github.com/pydata/xarray

# %% editable=true slideshow={"slide_type": ""}
import xarray as xr

time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)
ds_co2_yearly_global = xr.open_mfdataset(co2_yearly_global_fps, decode_times=time_coder)

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
# Force values to compute to avoid dask getting involved
ds_co2_yearly_global = ds_co2_yearly_global.compute()

# %% editable=true slideshow={"slide_type": ""}
ds_co2_yearly_global

# %% editable=true slideshow={"slide_type": ""} tags=["remove_output"]
fig, ax = plt.subplots(figsize=(6, 3))
ds_co2_yearly_global["co2"].plot.scatter(alpha=0.4, ax=ax)

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
glue("ds-co2-yearly-global-fig", fig, display=False)
plt.show()

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ```{glue:figure} ds-co2-yearly-global-fig
# ---
# width: 500px
# name: "ds-co2-yearly-global-fig"
# ---
#
# Atmospheric CO{raw-latex}`\textsubscript{2}` concentrations from year 1 to 2022
# in our dataset.
# ```

# %% [markdown] editable=true jp-MarkdownHeadingCollapsed=true slideshow={"slide_type": ""}
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
# as a scatter or a step plot
# ({numref}`Figure %s <ds_co2_yearly_step_fig>`).
# (The same logic applies to any spatial plots
# which could be created from our datasets
# that include spatial dimensions).

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
ds_plt = ds_co2_yearly_global.isel(time=slice(-5, None))

fig, ax = plt.subplots(figsize=(8, 4))
ds_plt["co2"].plot.scatter(ax=ax)

for bounds, val in zip(ds_plt["time_bnds"].values, ds_plt["co2"].values):
    ax.plot(bounds, [val, val], color="tab:blue", linewidth=1.0, alpha=0.7)

xticks = [cftime.DatetimeGregorian(y, 1, 1) for y in range(2018, 2024)]
ax.set_xticks(xticks)
ax.set_xlim(xticks[0], xticks[-1])
ax.grid()

glue("ds_co2_yearly_step_fig", fig, display=False)
plt.show()

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ```{glue:figure} ds_co2_yearly_step_fig
# ---
# width: 500px
# name: "ds_co2_yearly_step_fig"
# ---
#
# Illustration of the fact that each data point
# represents the average over its time bounds, not instantaneous values.
# As a result, it should be plotted with steps or scatters
# rather than an interpolated line.
# ```

# %% [markdown] editable=true jp-MarkdownHeadingCollapsed=true slideshow={"slide_type": ""}
# ## Monthly-, global-mean data
#
# If you want to have information at a finer level
# of temporal detail, we also provide monthly files.
# Like the global datasets, these come in three files.
#
# For monthly data, the time labels in the filename are months.
# Below we show the filenames for the CO{raw-latex}`\textsubscript{2}` output.

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
# Get file paths for the EO-updated data (see EO_UPDATE_DATA_ROOT above)
co2_monthly_global_fps = get_eo_update_local_files(
    ghg="co2", time_sampling="mon", grid="gm"
)

# %% editable=true slideshow={"slide_type": ""} tags=["remove_input"]
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

# %% [markdown] editable=true slideshow={"slide_type": ""}
# As above, as a result of the time average that the data represents,
# it is inappropriate to plot this data using a line plot.
# Scatter or step plots should be used instead
# ({numref}`Figure %s <ds-co2-monthly-global-fig>`).

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
ds_plt = ds_co2_monthly_global.isel(time=slice(-5 * 12, None))

fig, ax = plt.subplots(figsize=(8, 4))
ds_plt["co2"].plot.scatter(ax=ax)

for bounds, val in zip(ds_plt["time_bnds"].values, ds_plt["co2"].values):
    ax.plot(bounds, [val, val], color="tab:blue", linewidth=1.0, alpha=0.7)

xticks = [cftime.DatetimeGregorian(y, 1, 1) for y in range(2018, 2024)]
ax.set_xticks(xticks)
ax.set_xlim(xticks[0], xticks[-1])
ax.grid()

glue("ds-co2-monthly-global-fig", fig, display=False)
plt.show()

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ```{glue:figure} ds-co2-monthly-global-fig
# ---
# width: 500px
# name: "ds-co2-monthly-global-fig"
# ---
#
# Illustration of the monthly mean nature of the monthly datasets.
# Each value represents the average over its time bounds, not instantaneous values.
# As a result, they should be plotted with steps or scatters
# rather than an interpolated line.
# ```

# %% [markdown] editable=true slideshow={"slide_type": ""}
# The monthly data includes seasonality.
# Plotting the monthly and yearly data
# on the same axes makes particularly clear
# why a line plot is inappropriate
# ({numref}`Figure %s <ds-co2-monthly-global-vs-yearly-fig>`).

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
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

glue("ds-co2-monthly-global-vs-yearly-fig", fig, display=False)
plt.show()

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ```{glue:figure} ds-co2-monthly-global-vs-yearly-fig
# ---
# width: 500px
# name: "ds-co2-monthly-global-vs-yearly-fig"
# ---
#
# Monthly-mean compared to annual-mean datasets.
# Here illustrated with the CO{raw-latex}`\textsubscript{2}` dataset,
# but the same idea applies to all greenhouse gases.
# ```

# %% [markdown]
# At present, we do not provide data at a higher temporal resolution than monthly.
# In theory, such a dataset is possible to compile,
# however this requires careful consideration of daily
# and potentially sub-daily trends (e.g. the diurnal cycle).

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ## Monthly-, latitudinally-resolved data
#
# We also provide data with spatial,
# specifically latituindal, resolution.
# This data comes on a 15-degree latituindal grid
# (see below for details of the grid and latitudinal bounds).
# These files are identified by the grid label `gnz`.
# We only provide these files with monthly resolution.
#
# For completeness, we note that we also provide hemispheric means.
# These are not shown here,
# but are identified by the grid label `gr1z`.
#
# Below we show the filenames for the latitudinally-resolved data
# for CO{raw-latex}`\textsubscript{2}`

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
# Get file paths for the EO-updated data (see EO_UPDATE_DATA_ROOT above)
co2_monthly_lat_fps = get_eo_update_local_files(
    ghg="co2", time_sampling="mon", grid="gnz"
)

# %% editable=true slideshow={"slide_type": ""} tags=["remove_input"]
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

# %% [markdown] editable=true slideshow={"slide_type": ""}
# As above, but this time for the spatial axis,
# it is inappropriate to plot this data using a line plot.
# Scatter or step plots should be used instead
# ({numref}`Figure %s <ds-co2-monthly-lat-fig>`).

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
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

glue("ds-co2-monthly-lat-fig", fig, display=False)
plt.show()

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ```{glue:figure} ds-co2-monthly-lat-fig
# ---
# width: 600px
# name: "ds-co2-monthly-lat-fig"
# ---
#
# Illustration of the spatial mean nature of the latitudinally-resolved datasets
# (here shown for the year 2022 for CO{raw-latex}`\textsubscript{2}`
# but the same idea applies to all latitudinally-resolved datasets).
# Each value represents the average over its latitude bounds, not point values.
# As a result, they should be plotted with steps or scatters
# rather than an interpolated line.
# ```

# %% [markdown] editable=true slideshow={"slide_type": ""}
# We can compare the global-mean data
# to the data at each latitude
# ({numref}`Figure %s <ds-co2-monthly-lat-v-global-fig>`).
# The strength of the latitudinal gradient varies also by gas (not shown).

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
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

glue("ds-co2-monthly-lat-v-global-fig", fig, display=False)
plt.show()

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ```{glue:figure} ds-co2-monthly-lat-v-global-fig
# ---
# width: 550px
# name: "ds-co2-monthly-lat-v-global-fig"
# ---
#
# Comparison of global-, monthly-mean data with latitudinally-resolved, monthly-mean data
# for CO{raw-latex}`\textsubscript{2}`.
# The latitudinal variation, particularly the inverted seasonality in the two hemispheres,
# is a notable feature of the dataset.
# ```

# %% [markdown] editable=true slideshow={"slide_type": ""}
# The data can also be plotted in a so-called "magic carpet"
# to see the variation in space and time simultaneously
# ({numref}`Figure %s <ds-co2-magic-carpet-fig>`).


# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
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
glue("ds-co2-magic-carpet-fig", fig, display=False)
plt.show()

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ```{glue:figure} ds-co2-magic-carpet-fig
# ---
# width: 500px
# name: "ds-co2-magic-carpet-fig"
# ---
#
# So-called 'magic carpet' plot.
# This illustrates the variation in time and space simultaneously.
# Here this is illustrated with the CO{raw-latex}`\textsubscript{2}` dataset.
# ```

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
# Get file paths for the EO-updated data (see EO_UPDATE_DATA_ROOT above)
ch4_monthly_lat_fps = get_eo_update_local_files(
    ghg="ch4", time_sampling="mon", grid="gnz"
)
ds_ch4_monthly_lat = xr.open_mfdataset(
    ch4_monthly_lat_fps, decode_times=time_coder, data_vars=None, compat="no_conflicts"
)
ds_ch4_monthly_lat = ds_ch4_monthly_lat.compute()

ch4_monthly_global_fps = get_eo_update_local_files(
    ghg="ch4", time_sampling="mon", grid="gm"
)
ds_ch4_monthly_global = xr.open_mfdataset(
    ch4_monthly_global_fps, decode_times=time_coder
)
ds_ch4_monthly_global = ds_ch4_monthly_global.compute()

# %% [markdown] editable=true slideshow={"slide_type": ""}
# As for CO{raw-latex}`\textsubscript{2}`,
# the latitudinally-resolved data should be plotted
# with steps or scatters rather than an interpolated line
# ({numref}`Figure %s <ds-ch4-monthly-lat-fig>`).

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
ds_plt = ds_ch4_monthly_lat.isel(time=slice(-12, None))

mosaic_flat = [get_label_for_month(ds_plt.sel(time=time)) for time in ds_plt["time"]]

mosaic = [mosaic_flat[3 * i : 3 * (i + 1)] for i in range(len(mosaic_flat) // 3)]

fig, axes_d = plt.subplot_mosaic(mosaic, figsize=(8, 9), sharey=True, sharex=True)

for time in ds_plt["time"]:
    ds_plt_time = ds_plt.sel(time=time)
    label = get_label_for_month(ds_plt_time)

    axes_d[label].scatter(
        x=ds_plt_time["ch4"].values,
        y=ds_plt_time["lat"].values,
        s=10,
        label=label,
    )

    for bounds, val in zip(ds_plt_time["lat_bnds"].values, ds_plt_time["ch4"].values):
        axes_d[label].plot(
            [val, val], bounds, color="tab:blue", linewidth=1.0, alpha=0.7
        )

    yticks = np.arange(-90, 91, 15.0)
    axes_d[label].set_yticks(yticks)
    axes_d[label].set_ylim(yticks[0], yticks[-1])
    axes_d[label].grid()
    axes_d[label].set_title(label, fontsize="small")

for month in [1, 4, 7, 10]:
    axes_d[f"2022 - {calendar.month_name[month]}"].set_ylabel(
        "Latitude (degrees north)"
    )

for month in range(10, 13):
    axes_d[f"2022 - {calendar.month_name[month]}"].set_xlabel("ch4 [ppb]")

plt.tight_layout()

glue("ds-ch4-monthly-lat-fig", fig, display=False)
plt.show()

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ```{glue:figure} ds-ch4-monthly-lat-fig
# ---
# width: 600px
# name: "ds-ch4-monthly-lat-fig"
# ---
#
# Illustration of the spatial mean nature of the latitudinally-resolved datasets
# (here shown for the year 2022 for CH{raw-latex}`\textsubscript{4}`
# but the same idea applies to all latitudinally-resolved datasets).
# Each value represents the average over its latitude bounds, not point values.
# As a result, they should be plotted with steps or scatters
# rather than an interpolated line.
# ```

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
fig, ax = plt.subplots(figsize=(8, 4))

time_slice = slice(-5 * 12, None)

ds_plt = ds_ch4_monthly_global.isel(time=time_slice)
ds_plt["ch4"].plot.scatter(
    ax=ax, label="global-mean", color="tab:blue", s=30, zorder=10.0
)

for bounds, val in zip(ds_plt["time_bnds"].values, ds_plt["ch4"].values):
    ax.plot(bounds, [val, val], color="tab:blue", linewidth=1.0, alpha=0.7)

ds_all_lats = ds_ch4_monthly_lat.isel(time=time_slice)

for i, lat in enumerate(sorted(ds_all_lats["lat"])[::-1]):
    ds_plt = ds_all_lats.sel(lat=lat)
    colour = matplotlib.colormaps["magma"](i / len(ds_all_lats["lat"]))

    ds_plt["ch4"].plot.scatter(
        ax=ax, label=f"{float(lat)}", color=colour, marker="x", s=10
    )

    for bounds, val in zip(ds_plt["time_bnds"].values, ds_plt["ch4"].values):
        ax.plot(bounds, [val, val], color=colour, linewidth=1.0, alpha=0.7)

ax.legend(loc="center left", bbox_to_anchor=(1.05, 0.5))

xticks = [cftime.DatetimeGregorian(y, 1, 1) for y in range(2018, 2024)]
ax.set_xticks(xticks)
ax.set_xlim(xticks[0], xticks[-1])
ax.grid()

glue("ds-ch4-monthly-lat-v-global-fig", fig, display=False)
plt.show()

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ```{glue:figure} ds-ch4-monthly-lat-v-global-fig
# ---
# width: 550px
# name: "ds-ch4-monthly-lat-v-global-fig"
# ---
#
# Comparison of global-, monthly-mean data with latitudinally-resolved, monthly-mean data
# for CH{raw-latex}`\textsubscript{4}`.
# The latitudinal variation, particularly the inverted seasonality in the two hemispheres,
# is a notable feature of the dataset here too.
# ```

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(projection="3d")

tmp = ds_ch4_monthly_lat["ch4"].isel(time=range(-10 * 12, 0)).copy()
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
)

ax.view_init(15, -135, 0)

plt.tight_layout()
glue("ds-ch4-magic-carpet-fig", fig, display=False)
plt.show()

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ```{glue:figure} ds-ch4-magic-carpet-fig
# ---
# width: 500px
# name: "ds-ch4-magic-carpet-fig"
# ---
#
# So-called 'magic carpet' plot.
# This illustrates the variation in time and space simultaneously.
# Here this is illustrated with the CH{raw-latex}`\textsubscript{4}` dataset.
# ```

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ## Differences from CMIP6 and CMIP7
#
# ### File formats and naming
#
# The file formats are the same as those of the officially provided CMIP7 concentration files, which are generally close to CMIP6. We stress again the important difference that only data for CH{raw-latex}`\textsubscript{4}` and CO{raw-latex}`\textsubscript{2}` are provided.
# As in CMIP6, we do not provide any vertical profiles.
# For users who require such profiles,
# we refer to the 'The vertical dimension' sub-header
# in Section 4 of {raw-latex}`\textcite{meinshausen_historical_2017}`.
# There are three key changes:
#
# 1. we have split the global-mean and hemispheric-mean data into separate files.
#    In CMIP6, this data was in the same file (with a grid label of `GMNHSH`).
#    We have split this for two reasons:
#    a) `GMNHSH` is not a grid label recognised in the CMIP CVs
#       {raw-latex}`\parencite{wcrp_cmip_cvs_mip}` and
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
# }
# ```
#
# {raw-latex}`\newpage`

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ### Data comparisons with CMIP6
#
# Comparing the data from CMIP6 and from this data (CMIP7-like, with adedd EO information, referred to from now on as EO-CMIP7) shows minor changes
# (although doing this comparison requires a bit of care
# because of the changes in file formats).

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
gases_to_show = ["co2", "ch4"]


def load_eo_update_dataset(ghg, time_sampling, grid):
    """Load the EO-updated (satellite-extended) dataset for a single gas"""
    fps = get_eo_update_local_files(ghg=ghg, time_sampling=time_sampling, grid=grid)
    ds = xr.open_mfdataset(fps, decode_times=time_coder)

    # Unify time axis days to simplify
    ds["time"] = [
        cftime.DatetimeProlepticGregorian(v.year, v.month, 15)
        for v in ds["time"].values
    ]

    return ds.compute()


ds_gases_full_d = {}
for gas in gases_to_show:
    # "EO-CMIP7" is our local, EO-updated data (see EO_UPDATE_DATA_ROOT above)
    ds_gases_full_d[gas] = {"EO-CMIP7": load_eo_update_dataset(gas, "yr", "gm")}

    query_kwargs = {
        "ghg": gas,
        "time_sampling": "yr",
        "grid": "gm",
        "source_id": "UoM-CMIP-1-2-0",
        "cmip_era": "CMIP6",
        "engine": engine,
    }
    ds = fetch_and_load(**query_kwargs)

    # Unify time axis days to simplify
    ds["time"] = [
        cftime.DatetimeProlepticGregorian(v.year, v.month, 15)
        for v in ds["time"].values
    ]

    # compute to avoid dask weirdness
    ds_gases_full_d[gas]["CMIP6"] = ds.compute()

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
    era_a: str,
    era_b: str,
):
    """
    Plot overviews of timeseries and deltas between era_b and era_a
    """
    for ax_name, ax in axes_d.items():
        if ax_name.endswith("_delta"):
            continue

        gas = ax_name

        for era, ds in ds_d[gas].items():
            label = f"{era} ({ds.attrs['source_id']})"
            ds[gas].plot.scatter(
                ax=axes_d[gas], label=label, alpha=0.7, edgecolors="none"
            )

        ax.legend()
        ax.set_title(gas)
        ax.xaxis.set_tick_params(labelbottom=True)

        ax_delta = axes_d[f"{gas}_delta"]

        da_a = ds_d[gas][era_a][gas]
        da_b = ds_d[gas][era_b][gas]
        overlapping_times = np.intersect1d(da_a["time"], da_b["time"])
        delta = da_b.sel(time=overlapping_times) - da_a.sel(time=overlapping_times)
        # xarray drops attrs (including units) on arithmetic by default,
        # so the delta plot would otherwise show an unlabelled axis
        delta.attrs["units"] = da_a.attrs["units"]
        delta.attrs["long_name"] = f"{gas} difference"
        ax_delta.set_title(f"{era_b} - {era_a}", fontsize="small")
        delta.plot.scatter(
            ax=ax_delta,
            color="tab:grey",
            edgecolors="none",
            s=10,
        )
        ax_delta.axhline(0.0, color="k", linestyle="--")

        ax_delta.xaxis.set_tick_params(labelbottom=True)


plt_mosaic = [
    ["co2", "ch4"],
    ["co2", "ch4"],
    ["co2_delta", "ch4_delta"],
]
get_default_delta_mosaic = partial(
    plt.subplot_mosaic,
    mosaic=plt_mosaic,
    figsize=(8, 5),
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

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
fig, axes_d = get_default_delta_mosaic()
axes_d = remove_empty_axes(axes_d)

plot_overview_and_deltas(
    ds_gases_full_d,
    axes_d,
    era_a="CMIP6",
    era_b="EO-CMIP7",
)

plt.tight_layout()
glue("cmip6-v-cmip7-year-1-2022-fig", fig, display=False)
plt.show()

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ```{glue:figure} cmip6-v-cmip7-year-1-2022-fig
# ---
# width: 600px
# name: "cmip6-v-cmip7-year-1-2022-fig"
# ---
#
# CMIP6 vs. EO-CMIP7 for the time period from year 1 to 2022.
# For each pair of plots, the top panel shows the datasets' absolute values,
# the bottom panel shows the difference (EO-CMIP7 minus CMIP6).
# We show CO{raw-latex}`\textsubscript{2}` and
# CH{raw-latex}`\textsubscript{4}`,
# the two gases covered by this EO-updated dataset.
# ```
#
# {raw-latex}`\newpage`

# %% [markdown] editable=true slideshow={"slide_type": ""}
# #### Atmospheric concentrations: Year 1750 - 2022

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
fig, axes_d = get_default_delta_mosaic()
axes_d = remove_empty_axes(axes_d)

min_year = 1750
plot_overview_and_deltas(
    sel_times(ds_gases_full_d, lambda x: x.dt.year >= min_year),
    axes_d,
    era_a="CMIP6",
    era_b="EO-CMIP7",
)

plt.tight_layout()
glue("cmip6-v-cmip7-year-1750-2022-fig", fig, display=False)
plt.show()

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ```{glue:figure} cmip6-v-cmip7-year-1750-2022-fig
# ---
# width: 600px
# name: "cmip6-v-cmip7-year-1750-2022-fig"
# ---
#
# CMIP6 vs. EO-CMIP7 for the time period from year 1750 to 2022.
# For each pair of plots, the top panel shows the datasets' absolute values,
# the bottom panel shows the difference (EO-CMIP7 minus CMIP6).
# We show CO{raw-latex}`\textsubscript{2}` and
# CH{raw-latex}`\textsubscript{4}`,
# the two gases covered by this EO-updated dataset.
# ```
#
# {raw-latex}`\newpage`

# %% [markdown] editable=true slideshow={"slide_type": ""}
# #### Atmospheric concentrations: Year 1957 - 2022
#
# 1957 is the start of the Scripps ground-based record.
# Before this, data is based on ice cores alone.

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
fig, axes_d = get_default_delta_mosaic()
axes_d = remove_empty_axes(axes_d)

min_year = 1957
plot_overview_and_deltas(
    sel_times(ds_gases_full_d, lambda x: x.dt.year >= min_year),
    axes_d,
    era_a="CMIP6",
    era_b="EO-CMIP7",
)

plt.tight_layout()
glue("cmip6-v-cmip7-year-1957-2022-fig", fig, display=False)
plt.show()

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ```{glue:figure} cmip6-v-cmip7-year-1957-2022-fig
# ---
# width: 600px
# name: "cmip6-v-cmip7-year-1957-2022-fig"
# ---
#
# CMIP6 vs. EO-CMIP7 for the time period from year 1957 to 2022.
# For each pair of plots, the top panel shows the datasets' absolute values,
# the bottom panel shows the difference (EO-CMIP7 minus CMIP6).
# We show CO{raw-latex}`\textsubscript{2}` and
# CH{raw-latex}`\textsubscript{4}`,
# the two gases covered by this EO-updated dataset.
# ```
#
# {raw-latex}`\newpage`

# %% [markdown] editable=true slideshow={"slide_type": ""}
# #### Approximate radiative effect: Year 1 - 2022
#
# As seen above, in atmospheric concentration terms
# the differences are small.
# However, this can be put on a common scale
# by comparing the differences in radiative effect terms
# ({numref}`Figure %s <cmip6-v-cmip7-year-1-2022-re-fig>`
# and {numref}`Figure %s <cmip6-v-cmip7-year-1750-2022-re-fig>`).
# This gives an approximation of the size of the difference
# that would be seen by an Earth System Model's (ESM's) radiation code.
# This uses basic linear approximations,
# assuming that the radiative effect of each gas
# is simply its atmospheric concentration multiplied by a constant.
# This isn't the same as effective radiative forcing (ERF).
# For that comparison, see the later sections focussed on ERF.

# %% [markdown]
# Values below come from Table 7.SM.7 of
# IPCC AR6 WG1 Ch. 7 Supplementary Material
# {raw-latex}`\parencite{IPCC_2021_WGI_Ch_7_SM}`.

# %% editable=true slideshow={"slide_type": ""}
from openscm_units import unit_registry

Q = unit_registry.Quantity

RADIATIVE_EFFICIENCIES = {
    "co2": Q(1.33e-5, "W / m^2 / ppb"),
    "ch4": Q(3.88e-4, "W / m^2 / ppb"),
}

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
ds_gases_full_radiative_effect_d = {}
target_units = "W / m^2"
for gas, gas_ds in ds_gases_full_d.items():
    ds_gases_full_radiative_effect_d[gas] = {}
    for mip_era, ds in gas_ds.items():
        tmp = ds.copy(deep=True)

        tmp[gas][:] = (
            (Q(tmp[gas].values, tmp[gas].attrs["units"]) * RADIATIVE_EFFICIENCIES[gas])
            .to(target_units)
            .m
        )
        tmp[gas].attrs["units"] = target_units
        tmp[gas].attrs["long_name"] = "approx. radiative effect"

        ds_gases_full_radiative_effect_d[gas][mip_era] = tmp

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
fig, axes_d = plt.subplot_mosaic(
    mosaic=plt_mosaic,
    figsize=(8, 6),
    sharex=True,
)
axes_d = remove_empty_axes(axes_d)

plot_overview_and_deltas(
    ds_gases_full_radiative_effect_d,
    axes_d,
    era_a="CMIP6",
    era_b="EO-CMIP7",
)

for name, ax in axes_d.items():
    if name.endswith("_delta"):
        continue

    ax.set_ylim([0, 6.0])

plt.tight_layout()
glue("cmip6-v-cmip7-year-1-2022-re-fig", fig, display=False)
plt.show()

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ```{glue:figure} cmip6-v-cmip7-year-1-2022-re-fig
# ---
# width: 600px
# name: "cmip6-v-cmip7-year-1-2022-re-fig"
# ---
#
# CMIP6 vs. EO-CMIP7 for the time period from year 1 to 2022
# in radiative effect terms.
# For each pair of plots, the top panel shows the datasets
# in radiative effect terms i.e. the product of the concentration
# and its radiative efficiency.
# The bottom panel shows the difference (EO-CMIP7 minus CMIP6).
# We show CO{raw-latex}`\textsubscript{2}` and
# CH{raw-latex}`\textsubscript{4}`,
# the two gases covered by this EO-updated dataset.
# ```
#
# {raw-latex}`\newpage`

# %% [markdown]
# #### Approximate radiative effect: Year 1750 - 2022
#
# This is the period relevant for historical simulations in CMIP.

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
fig, axes_d = plt.subplot_mosaic(
    mosaic=plt_mosaic,
    figsize=(8, 6),
    sharex=True,
)
axes_d = remove_empty_axes(axes_d)

min_year = 1750
plot_overview_and_deltas(
    sel_times(ds_gases_full_radiative_effect_d, lambda x: x.dt.year >= min_year),
    axes_d,
    era_a="CMIP6",
    era_b="EO-CMIP7",
)

for name, ax in axes_d.items():
    if name.endswith("_delta"):
        continue

    ax.set_ylim([0, 6.0])

plt.tight_layout()
glue("cmip6-v-cmip7-year-1750-2022-re-fig", fig, display=False)
plt.show()

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ```{glue:figure} cmip6-v-cmip7-year-1750-2022-re-fig
# ---
# width: 600px
# name: "cmip6-v-cmip7-year-1750-2022-re-fig"
# ---
#
# CMIP6 vs. EO-CMIP7 for the time period from year 1750 to 2022
# in radiative effect terms
# (see caption of
# {numref}`Figure %s <cmip6-v-cmip7-year-1-2022-re-fig>`
# for details).
# ```
#
# {raw-latex}`\newpage`

# %% [markdown] editable=true slideshow={"slide_type": ""}
# #### Approximate effective radiative forcing: Year 1750 - 2022
#
# The above isn't effective radiative forcing.
# For that, you have to normalise the data to some reference year.
# There are a few different choices for this reference year.
# In IPCC reports, it is 1750 so that is what we show here.
# It should be noted that some ESMs may make other choices,
# but these would not have a great effect on the interpretation
# of the difference between the CMIP6 and EO-CMIP7 datasets.
#
# Note that this approximation is linear,
# which is a particularly strong approximation for CO{raw-latex}`\textsubscript{2}`
# because of its logarithmic forcing nature.
# We show this approximation here
# ({numref}`Figure %s <cmip6-v-cmip7-year-1750-2022-erf-fig>`)
# nonetheless because it provides an order of magnitude estimate
# for the change from CMIP6 in ERF terms.

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
ds_gases_full_erf_d = {}
reference_year = 1750
for gas, gas_ds in ds_gases_full_radiative_effect_d.items():
    ds_gases_full_erf_d[gas] = {}
    for mip_era, ds in gas_ds.items():
        tmp = ds.copy(deep=True)

        tmp[gas][:] = (
            tmp[gas][:]
            - tmp.sel(time=ds["time"].dt.year == reference_year)[gas][:].values
        )
        tmp[gas].attrs["long_name"] = "approx. ERF"
        # tmp.attrs = ds.attrs

        ds_gases_full_erf_d[gas][mip_era] = tmp

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
fig, axes_d = get_default_delta_mosaic()
axes_d = remove_empty_axes(axes_d)

min_year = 1750
plot_overview_and_deltas(
    sel_times(ds_gases_full_erf_d, lambda x: x.dt.year >= min_year),
    axes_d,
    era_a="CMIP6",
    era_b="EO-CMIP7",
)

for name, ax in axes_d.items():
    if name.endswith("_delta"):
        continue

    ax.set_ylim([0, 2.0])

plt.tight_layout()
glue("cmip6-v-cmip7-year-1750-2022-erf-fig", fig, display=False)
plt.show()

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ```{glue:figure} cmip6-v-cmip7-year-1750-2022-erf-fig
# ---
# width: 600px
# name: "cmip6-v-cmip7-year-1750-2022-erf-fig"
# ---
#
# CMIP6 vs. EO-CMIP7 for the time period from year 1759 to 2022
# in approximate effective radiative forcing terms.
# For each pair of plots, the top panel shows the datasets
# in approximate effective radiative forcing terms.
# The bottom panel shows the difference (EO-CMIP7 minus CMIP6).
# We show CO{raw-latex}`\textsubscript{2}` and
# CH{raw-latex}`\textsubscript{4}`,
# the two gases covered by this EO-updated dataset.
# ```
#
# {raw-latex}`\newpage`

# %% [markdown]
# In summary, in ERF terms, the differences from CMIP6 are very small.
# For all gases, they are less than around 0.03 W / m{raw-latex}`\textsuperscript{2}`.
# Compared to the estimated total greenhouse gas forcing and uncertainty in IPCC AR6
# {raw-latex}`\parencite[see Section 7.3.5.2 of AR6 WG1 Chapter 7,][]{IPCC_2021_WGI_Ch_7}`
# estimated to be 3.84 W / m{raw-latex}`\textsuperscript{2}`
# (very likely range of 3.46 to 4.22 W / m{raw-latex}`\textsuperscript{2}`),
# such differences are particularly small.

# %% [markdown]
# #### Atmospheric concentrations including seasonality: Year 2000 - 2022
#
# The final comparisons we show are atmospheric concentrations including seasonality
# ({numref}`Figure %s <cmip6-v-cmip7-year-2000-2022-seasonality-fig>`).
# Given that most greenhouse gases
# are well-mixed with lifetimes much greater than a year,
# these differences are unlikely to be of huge interest to ESMs.
# However, for other applications, such seasonality differences may matter more.

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
ds_gases_full_monthly_d = {}
for gas in gases_to_show:
    # "EO-CMIP7" is our local, EO-updated data (see EO_UPDATE_DATA_ROOT above)
    ds_gases_full_monthly_d[gas] = {
        "EO-CMIP7": load_eo_update_dataset(gas, "mon", "gm")
    }

    query_kwargs = {
        "ghg": gas,
        "time_sampling": "mon",
        "grid": "gm",
        "source_id": "UoM-CMIP-1-2-0",
        "cmip_era": "CMIP6",
        "engine": engine,
    }
    ds = fetch_and_load(**query_kwargs)

    # Unify time axis days to simplify
    ds["time"] = [
        cftime.DatetimeProlepticGregorian(v.year, v.month, 15)
        for v in ds["time"].values
    ]

    # compute to avoid dask weirdness
    ds_gases_full_monthly_d[gas]["CMIP6"] = ds.compute()

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
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
    era_a="CMIP6",
    era_b="EO-CMIP7",
)

plt.tight_layout()
glue("cmip6-v-cmip7-year-2000-2022-seasonality-fig", fig, display=False)
plt.show()

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ```{glue:figure} cmip6-v-cmip7-year-2000-2022-seasonality-fig
# ---
# width: 600px
# name: "cmip6-v-cmip7-year-2000-2022-seasonality-fig"
# ---
#
# CMIP6 vs. EO-CMIP7 for the time period from year 2000 to 2022.
# The shown dataset is the global-, monthly-mean dataset
# i.e. includes seasonality.
# For each pair of plots, the top panel shows the datasets' absolute values,
# the bottom panel shows the difference (EO-CMIP7 minus CMIP6).
# We show CO{raw-latex}`\textsubscript{2}` and
# CH{raw-latex}`\textsubscript{4}`,
# the two gases covered by this EO-updated dataset.
# ```

# %% [markdown] editable=true slideshow={"slide_type": ""}
# Like the annual-means,
# the atmospheric concentrations including seasonality
# are reasonably consistent between CMIP6 and EO-CMIP7.
# There are some areas of change.
# Full details of these changes will be provided
# in future work .

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ### Data comparison with original CMIP7 concentrations
#
# Comparing the original CMIP7 data (without satellite information)
# and EO-CMIP7 (this dataset, with added EO information) shows minor changes in the case of CO{raw-latex}`\textsubscript{2}`. In case of CH{raw-latex}`\textsubscript{4}`, even though the magnitude of the changes is not high, some unphysical artefacts are introduced, specifically in the year 2003 (when satellite measurements start) and in 1948 (when the firn dataset stops).
# For an in-detail exploration and discussion of these differences, feel free to contact the authors of this user guide, but for now we do not recommend using the EO-CMIP7 concentrations for CH{raw-latex}`\textsubscript{4}` in your work.

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
ds_gases_full_cmip7_d = {}
for gas in gases_to_show:
    # Reuse the EO-updated data already loaded above
    ds_gases_full_cmip7_d[gas] = {
        "EO-CMIP7": ds_gases_full_d[gas]["EO-CMIP7"].copy(deep=True)
    }

    query_kwargs = {
        "ghg": gas,
        "time_sampling": "yr",
        "grid": "gm",
        "source_id": "CR-CMIP-1-0-0",
        "cmip_era": "CMIP7",
        "engine": engine,
    }
    ds = fetch_and_load(**query_kwargs)

    # Unify time axis days to simplify
    ds["time"] = [
        cftime.DatetimeProlepticGregorian(v.year, v.month, 15)
        for v in ds["time"].values
    ]

    # compute to avoid dask weirdness
    ds_gases_full_cmip7_d[gas]["CMIP7"] = ds.compute()

# %% [markdown] editable=true slideshow={"slide_type": ""}
# #### Atmospheric concentrations: Year 1 - 2022

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
fig, axes_d = get_default_delta_mosaic()
axes_d = remove_empty_axes(axes_d)

plot_overview_and_deltas(
    ds_gases_full_cmip7_d,
    axes_d,
    era_a="CMIP7",
    era_b="EO-CMIP7",
)

plt.tight_layout()
glue("cmip7-v-eo-cmip7-year-1-2022-fig", fig, display=False)
plt.show()

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ```{glue:figure} cmip7-v-eo-cmip7-year-1-2022-fig
# ---
# width: 600px
# name: "cmip7-v-eo-cmip7-year-1-2022-fig"
# ---
#
# Original CMIP7 vs. EO-CMIP7 for the time period from year 1 to 2022.
# For each pair of plots, the top panel shows the datasets' absolute values,
# the bottom panel shows the difference (EO-CMIP7 minus original CMIP7).
# We show CO{raw-latex}`\textsubscript{2}` and
# CH{raw-latex}`\textsubscript{4}`,
# the two gases covered by this EO-updated dataset.
# ```
#
# {raw-latex}`\newpage`

# %% [markdown] editable=true slideshow={"slide_type": ""}
# #### Atmospheric concentrations: Year 1750 - 2022

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
fig, axes_d = get_default_delta_mosaic()
axes_d = remove_empty_axes(axes_d)

min_year = 1750
plot_overview_and_deltas(
    sel_times(ds_gases_full_cmip7_d, lambda x: x.dt.year >= min_year),
    axes_d,
    era_a="CMIP7",
    era_b="EO-CMIP7",
)

plt.tight_layout()
glue("cmip7-v-eo-cmip7-year-1750-2022-fig", fig, display=False)
plt.show()

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ```{glue:figure} cmip7-v-eo-cmip7-year-1750-2022-fig
# ---
# width: 600px
# name: "cmip7-v-eo-cmip7-year-1750-2022-fig"
# ---
#
# Original CMIP7 vs. EO-CMIP7 for the time period from year 1750 to 2022.
# For each pair of plots, the top panel shows the datasets' absolute values,
# the bottom panel shows the difference (EO-CMIP7 minus original CMIP7).
# We show CO{raw-latex}`\textsubscript{2}` and
# CH{raw-latex}`\textsubscript{4}`,
# the two gases covered by this EO-updated dataset.
# ```
#
# {raw-latex}`\newpage`

# %% [markdown] editable=true slideshow={"slide_type": ""}
# #### Atmospheric concentrations: Year 1957 - 2022
#
# 1957 is the start of the Scripps ground-based record.
# Before this, data is based on ice cores alone.

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
fig, axes_d = get_default_delta_mosaic()
axes_d = remove_empty_axes(axes_d)

min_year = 1957
plot_overview_and_deltas(
    sel_times(ds_gases_full_cmip7_d, lambda x: x.dt.year >= min_year),
    axes_d,
    era_a="CMIP7",
    era_b="EO-CMIP7",
)

plt.tight_layout()
glue("cmip7-v-eo-cmip7-year-1957-2022-fig", fig, display=False)
plt.show()

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ```{glue:figure} cmip7-v-eo-cmip7-year-1957-2022-fig
# ---
# width: 600px
# name: "cmip7-v-eo-cmip7-year-1957-2022-fig"
# ---
#
# Original CMIP7 vs. EO-CMIP7 for the time period from year 1957 to 2022.
# For each pair of plots, the top panel shows the datasets' absolute values,
# the bottom panel shows the difference (EO-CMIP7 minus original CMIP7).
# We show CO{raw-latex}`\textsubscript{2}` and
# CH{raw-latex}`\textsubscript{4}`,
# the two gases covered by this EO-updated dataset.
# ```
#
# {raw-latex}`\newpage`

# %% [markdown] editable=true slideshow={"slide_type": ""}
# #### Approximate radiative effect: Year 1 - 2022
#
# As seen above, in atmospheric concentration terms
# the differences are small.
# However, this can be put on a common scale
# by comparing the differences in radiative effect terms
# ({numref}`Figure %s <cmip7-v-eo-cmip7-year-1-2022-re-fig>`
# and {numref}`Figure %s <cmip7-v-eo-cmip7-year-1750-2022-re-fig>`).
# This gives an approximation of the size of the difference
# that would be seen by an Earth System Model's (ESM's) radiation code.
# This uses basic linear approximations,
# assuming that the radiative effect of each gas
# is simply its atmospheric concentration multiplied by a constant.
# This isn't the same as effective radiative forcing (ERF).
# For that comparison, see the later sections focussed on ERF.

# %% [markdown]
# Values below come from Table 7.SM.7 of
# IPCC AR6 WG1 Ch. 7 Supplementary Material
# {raw-latex}`\parencite{IPCC_2021_WGI_Ch_7_SM}`.

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
ds_gases_full_cmip7_radiative_effect_d = {}
target_units = "W / m^2"
for gas, gas_ds in ds_gases_full_cmip7_d.items():
    ds_gases_full_cmip7_radiative_effect_d[gas] = {}
    for mip_era, ds in gas_ds.items():
        tmp = ds.copy(deep=True)

        tmp[gas][:] = (
            (Q(tmp[gas].values, tmp[gas].attrs["units"]) * RADIATIVE_EFFICIENCIES[gas])
            .to(target_units)
            .m
        )
        tmp[gas].attrs["units"] = target_units
        tmp[gas].attrs["long_name"] = "approx. radiative effect"

        ds_gases_full_cmip7_radiative_effect_d[gas][mip_era] = tmp

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
fig, axes_d = plt.subplot_mosaic(
    mosaic=plt_mosaic,
    figsize=(8, 6),
    sharex=True,
)
axes_d = remove_empty_axes(axes_d)

plot_overview_and_deltas(
    ds_gases_full_cmip7_radiative_effect_d,
    axes_d,
    era_a="CMIP7",
    era_b="EO-CMIP7",
)

for name, ax in axes_d.items():
    if name.endswith("_delta"):
        continue

    ax.set_ylim([0, 6.0])

plt.tight_layout()
glue("cmip7-v-eo-cmip7-year-1-2022-re-fig", fig, display=False)
plt.show()

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ```{glue:figure} cmip7-v-eo-cmip7-year-1-2022-re-fig
# ---
# width: 600px
# name: "cmip7-v-eo-cmip7-year-1-2022-re-fig"
# ---
#
# Original CMIP7 vs. EO-CMIP7 for the time period from year 1 to 2022
# in radiative effect terms.
# For each pair of plots, the top panel shows the datasets
# in radiative effect terms i.e. the product of the concentration
# and its radiative efficiency.
# The bottom panel shows the difference (EO-CMIP7 minus original CMIP7).
# We show CO{raw-latex}`\textsubscript{2}` and
# CH{raw-latex}`\textsubscript{4}`,
# the two gases covered by this EO-updated dataset.
# ```
#
# {raw-latex}`\newpage`

# %% [markdown]
# #### Approximate radiative effect: Year 1750 - 2022
#
# This is the period relevant for historical simulations in CMIP.

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
fig, axes_d = plt.subplot_mosaic(
    mosaic=plt_mosaic,
    figsize=(8, 6),
    sharex=True,
)
axes_d = remove_empty_axes(axes_d)

min_year = 1750
plot_overview_and_deltas(
    sel_times(ds_gases_full_cmip7_radiative_effect_d, lambda x: x.dt.year >= min_year),
    axes_d,
    era_a="CMIP7",
    era_b="EO-CMIP7",
)

for name, ax in axes_d.items():
    if name.endswith("_delta"):
        continue

    ax.set_ylim([0, 6.0])

plt.tight_layout()
glue("cmip7-v-eo-cmip7-year-1750-2022-re-fig", fig, display=False)
plt.show()

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ```{glue:figure} cmip7-v-eo-cmip7-year-1750-2022-re-fig
# ---
# width: 600px
# name: "cmip7-v-eo-cmip7-year-1750-2022-re-fig"
# ---
#
# Original CMIP7 vs. EO-CMIP7 for the time period from year 1750 to 2022
# in radiative effect terms
# (see caption of
# {numref}`Figure %s <cmip7-v-eo-cmip7-year-1-2022-re-fig>`
# for details).
# ```
#
# {raw-latex}`\newpage`

# %% [markdown] editable=true slideshow={"slide_type": ""}
# #### Approximate effective radiative forcing: Year 1750 - 2022
#
# The above isn't effective radiative forcing.
# For that, you have to normalise the data to some reference year.
# There are a few different choices for this reference year.
# In IPCC reports, it is 1750 so that is what we show here.
# It should be noted that some ESMs may make other choices,
# but these would not have a great effect on the interpretation
# of the difference between the original CMIP7 and EO-CMIP7 datasets.
#
# Note that this approximation is linear,
# which is a particularly strong approximation for CO{raw-latex}`\textsubscript{2}`
# because of its logarithmic forcing nature.
# We show this approximation here
# ({numref}`Figure %s <cmip7-v-eo-cmip7-year-1750-2022-erf-fig>`)
# nonetheless because it provides an order of magnitude estimate
# for the change from original CMIP7 in ERF terms.
# The forthcoming manuscripts will explore the subtleties
# of this quantification in more detail.

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
ds_gases_full_cmip7_erf_d = {}
reference_year = 1750
for gas, gas_ds in ds_gases_full_cmip7_radiative_effect_d.items():
    ds_gases_full_cmip7_erf_d[gas] = {}
    for mip_era, ds in gas_ds.items():
        tmp = ds.copy(deep=True)

        tmp[gas][:] = (
            tmp[gas][:]
            - tmp.sel(time=ds["time"].dt.year == reference_year)[gas][:].values
        )
        tmp[gas].attrs["long_name"] = "approx. ERF"

        ds_gases_full_cmip7_erf_d[gas][mip_era] = tmp

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
fig, axes_d = get_default_delta_mosaic()
axes_d = remove_empty_axes(axes_d)

min_year = 1750
plot_overview_and_deltas(
    sel_times(ds_gases_full_cmip7_erf_d, lambda x: x.dt.year >= min_year),
    axes_d,
    era_a="CMIP7",
    era_b="EO-CMIP7",
)

for name, ax in axes_d.items():
    if name.endswith("_delta"):
        continue

    ax.set_ylim([0, 2.0])

plt.tight_layout()
glue("cmip7-v-eo-cmip7-year-1750-2022-erf-fig", fig, display=False)
plt.show()

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ```{glue:figure} cmip7-v-eo-cmip7-year-1750-2022-erf-fig
# ---
# width: 600px
# name: "cmip7-v-eo-cmip7-year-1750-2022-erf-fig"
# ---
#
# Original CMIP7 vs. EO-CMIP7 for the time period from year 1750 to 2022
# in approximate effective radiative forcing terms.
# For each pair of plots, the top panel shows the datasets
# in approximate effective radiative forcing terms.
# The bottom panel shows the difference (EO-CMIP7 minus original CMIP7).
# We show CO{raw-latex}`\textsubscript{2}` and
# CH{raw-latex}`\textsubscript{4}`,
# the two gases covered by this EO-updated dataset.
# ```
#
# {raw-latex}`\newpage`

# %% [markdown]
# In summary, in ERF terms, the differences from original CMIP7 are very small.
# For all gases, they are less than around 0.04 W / m{raw-latex}`\textsuperscript{2}`.
# Compared to the estimated total greenhouse gas forcing and uncertainty in IPCC AR6
# {raw-latex}`\parencite[see Section 7.3.5.2 of AR6 WG1 Chapter 7,][]{IPCC_2021_WGI_Ch_7}`
# estimated to be 3.84 W / m{raw-latex}`\textsuperscript{2}`
# (very likely range of 3.46 to 4.22 W / m{raw-latex}`\textsuperscript{2}`),
# such differences are particularly small.

# %% [markdown]
# #### Atmospheric concentrations including seasonality: Year 2000 - 2022
#
# The final comparisons we show are atmospheric concentrations including seasonality
# ({numref}`Figure %s <cmip7-v-eo-cmip7-year-2000-2022-seasonality-fig>`).
# Given that most greenhouse gases
# are well-mixed with lifetimes much greater than a year,
# these differences are unlikely to be of huge interest to ESMs.
# However, for other applications, such seasonality differences may matter more.

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
ds_gases_full_cmip7_monthly_d = {}
for gas in gases_to_show:
    # Reuse the EO-updated data already loaded above
    ds_gases_full_cmip7_monthly_d[gas] = {
        "EO-CMIP7": ds_gases_full_monthly_d[gas]["EO-CMIP7"].copy(deep=True)
    }

    query_kwargs = {
        "ghg": gas,
        "time_sampling": "mon",
        "grid": "gm",
        "source_id": "CR-CMIP-1-0-0",
        "cmip_era": "CMIP7",
        "engine": engine,
    }
    ds = fetch_and_load(**query_kwargs)

    # Unify time axis days to simplify
    ds["time"] = [
        cftime.DatetimeProlepticGregorian(v.year, v.month, 15)
        for v in ds["time"].values
    ]

    # compute to avoid dask weirdness
    ds_gases_full_cmip7_monthly_d[gas]["CMIP7"] = ds.compute()

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
fig, axes_d = get_default_delta_mosaic()
axes_d = remove_empty_axes(axes_d)

min_year = 2000
max_year = 2022
plot_overview_and_deltas(
    sel_times(
        ds_gases_full_cmip7_monthly_d,
        lambda x: (x.dt.year >= min_year) & (x.dt.year <= max_year),
    ),
    axes_d,
    era_a="CMIP7",
    era_b="EO-CMIP7",
)

plt.tight_layout()
glue("cmip7-v-eo-cmip7-year-2000-2022-seasonality-fig", fig, display=False)
plt.show()

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ```{glue:figure} cmip7-v-eo-cmip7-year-2000-2022-seasonality-fig
# ---
# width: 600px
# name: "cmip7-v-eo-cmip7-year-2000-2022-seasonality-fig"
# ---
#
# Original CMIP7 vs. EO-CMIP7 for the time period from year 2000 to 2022.
# The shown dataset is the global-, monthly-mean dataset
# i.e. includes seasonality.
# For each pair of plots, the top panel shows the datasets' absolute values,
# the bottom panel shows the difference (EO-CMIP7 minus original CMIP7).
# We show CO{raw-latex}`\textsubscript{2}` and
# CH{raw-latex}`\textsubscript{4}`,
# the two gases covered by this EO-updated dataset.
# ```

# %% [markdown] editable=true slideshow={"slide_type": ""}
# Like the annual-means,
# the atmospheric concentrations including seasonality
# are reasonably consistent for CO{raw-latex}`\textsubscript{2}` between the original CMIP7 and EO-CMIP7 datasets.
# There are some areas of change. For CH{raw-latex}`\textsubscript{4}`, the magnitute of change is also limited but non-physical artifacts are introduced by the addition of satellite data, so we do not recommend using CH{raw-latex}`\textsubscript{4}` concentrations from the EO-CMIP7 dataset at its current state.
#
