# ---
# jupyter:
#   authors:
#   - name: Anna Lanteri
#   - name: Zebedee Nicholls
#   - name: Florence Bockting
#   - name: "Mika Pfl\xFCger"
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
#
#     Scenario Dataset: Data Description and User Guide'
# ---

# %% [markdown]
# ```{role} raw-latex(raw)
# :format: latex
# ```

# %% [markdown] editable=true slideshow={"slide_type": ""} jp-MarkdownHeadingCollapsed=true
# # Overview
#
# Here we provide a short description of the scenario dataset of concentrations
# for CO{raw-latex}`\textsubscript{2}` and CH{raw-latex}`\textsubscript{4}`
# updated with satellite data (EO-Scenario), and a guide for users.
# This is intended to provide a short introduction for users of the data,
# and a comparison with the baseline CMIP7 Scenarios forcing datasets.
#
# The dataset is an extension on the forcings for CMIP7 scenarios datasets for
# CO{raw-latex}`\textsubscript{2}` and CH{raw-latex}`\textsubscript{4}`,
# achieved by pre-processing and including satellite data to the historical GHG concentration
# generation pipeline, and then using those concentrations in the CMIP7 scenario concentration
# generation pipeline.
#
# In this document, we will describe both the
# original dataset and the updated version, to have a complete and self-contained description. This means this document contains
# duplicated information from the user guide: ‘CMIP7 Greenhouse Gas (GHG) Concentration Forcing Scenarios Dataset’

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
from myst_nb import glue

from local.data_loading import (
    fetch_and_load_ghg_dataset_scenarios,
    get_scenario,
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

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
# TODO: update this to point at the final, published EO-update scenario dataset
# once it exists (e.g. once it's on ESGF and can be fetched the same way
# as the original dataset). For now, this points directly at a local
# output bundle from cmip7-scenariomip-ghg-concentrations, so there is
# no ESGF fetch step for the EO-updated data (CO2 and CH4 only).
import xarray as xr

EO_UPDATE_DATA_ROOT = (
    REPO_ROOT.parent
    / "cmip7-scenariomip-ghg-concentrations"
    / "output-bundles"
    / "local-fit-satweighted"
    / "data"
    / "processed"
    / "esgf-ready"
    / "input4MIPs"
    / "CMIP7"
    / "ScenarioMIP"
    / "CR"
)
EO_UPDATE_DATA_VERSION = "v20261001"
EO_UPDATE_GASES = ("co2", "ch4")

# Historical data: the EO-updated historical dataset ("EO-CMIP7"),
# i.e. the one used in the historical EO-update user guide.
# TODO: also update this once the dataset is published.
EO_UPDATE_HISTORY_ROOT = (
    REPO_ROOT.parent
    / "CMIP-GHG-Concentration-Generation"
    / "output-bundles"
    / "eo-update"
    / "data"
    / "processed"
    / "esgf-ready"
    / "EO-update-input4MIPs"
    / "CMIP7Plus"
    / "CMIP"
    / "CR"
    / "CR-CMIP7Plus"
    / "atmos"
)
EO_UPDATE_HISTORY_VERSION = "v20261008"


def load_eo_update_history(ghg, time_sampling, grid, with_scenario_dim=False):
    """
    Load the EO-updated historical dataset

    If `with_scenario_dim` is `True`,
    the result has a `scenario` dimension (of length one, labelled "CMIP")
    to match what `fetch_and_load_ghg_dataset_scenarios` returns
    for the original history.
    """
    d = EO_UPDATE_HISTORY_ROOT / time_sampling / ghg / grid / EO_UPDATE_HISTORY_VERSION
    ds = xr.open_mfdataset(
        sorted(d.glob("*.nc")),
        decode_times=xr.coders.CFDatetimeCoder(use_cftime=True),
    )
    if with_scenario_dim:
        ds = xr.concat([ds.assign_coords(scenario="CMIP")], "scenario")

    return ds


def get_eo_update_local_files(ghg, time_sampling, grid):
    """Get the local files for the EO-updated (satellite-extended) scenario dataset"""
    return tuple(
        sorted(
            EO_UPDATE_DATA_ROOT.glob(
                f"CR-*-1-1-0/atmos/{time_sampling}/{ghg}/{grid}/"
                f"{EO_UPDATE_DATA_VERSION}/*.nc"
            )
        )
    )


def load_eo_update_scenarios(ghg, time_sampling, grid, include_ext=True):
    """
    Load all scenarios of the EO-updated dataset

    Mirrors the output of `fetch_and_load_ghg_dataset_scenarios`
    i.e. one dataset with a `scenario` dimension.
    If `include_ext` is `False`, the extension scenarios are skipped
    (saves a lot of time if they are going to be dropped anyway).
    """
    time_coder_eo = xr.coders.CFDatetimeCoder(use_cftime=True)

    fps_by_source_id = {}
    for fp in get_eo_update_local_files(ghg, time_sampling, grid):
        source_id = fp.relative_to(EO_UPDATE_DATA_ROOT).parts[0]
        if not include_ext and "-ext-" in source_id:
            continue

        fps_by_source_id.setdefault(source_id, []).append(fp)

    res_l = []
    for fps in fps_by_source_id.values():
        tmp = xr.open_mfdataset(fps, decode_times=time_coder_eo)
        tmp = tmp.assign_coords(scenario=get_scenario(tmp))
        res_l.append(tmp)

    return xr.concat(res_l, "scenario", join="outer")


# %% [markdown] editable=true slideshow={"slide_type": ""} jp-MarkdownHeadingCollapsed=true
# # Dataset construction
#
# The original CMIP7 scenarios GHG concentration dataset (produced without using satellite data)
# is constructed following a similar methodology
# to {raw-latex}`\textcite{meinshausen_shared_2020}`.
# The full method will be described in a forthcoming paper.
# In brief, the method is:
#
# 1. **retrieve concentrations** of CO{raw-latex}`\textsubscript{2}` and CH{raw-latex}`\textsubscript{4}`
#
#     1. start with harmonised, complete emissions sets for each scenario
#        of interest.
#
#        - These are not publicly available yet, but will be released soon
#          as part of the ongoing ScenarioMIP process
#
#     1. run MAGICC
#        {raw-latex}`\parencite{meinshausen2011emulating1,meinshausen2011emulating2}`
#        to translate these emissions into global-mean concentrations
#
#        1. MAGICC is run in the same configuration which was used in AR6
#          {raw-latex}`\parencite[as described/evaluated in Cross-Chapter Box 7.1][]{IPCC_2021_WGI_Ch_7}`.
#          This represents our best estimate, in line with the last IPCC report,
#          of the concentrations that result from the emissions.
#          Note, given that CMIP7 models have not yet been run, this will,
#          almost by definition, produce different concentrations
#          than those outputed by the ESMs run in emissions-driven mode.
#     1. harmonise the global-mean concentrations to the historical
#        concentrations using gradient-aware harmonisation[^gah] to ensure a smooth
#        transition in both the absolute values and the gradient, improving on the
#        CMIP6 data which had an abrupt jump in the gradient (most notable for
#        methane, further details below)
#
# 3. **seasonality and latitudinal gradients**: We use the same statistical
#    models for seasonality and latitudinal gradients as were used over the
#    historical period
#    (these models are based on orthogonal function methods, linear regression
#    and basic scaling arguments).
#    These are then applied to the future, based on the same input drivers.
#     - one exception is CO{raw-latex}`\textsubscript{2}`, which uses GPP as the input driver
#       rather than a regression on a combined global-mean surface temperature and
#       CO{raw-latex}`\textsubscript{2}` concentration metric (as is used in the product that covers
#       history)
#     - arguably we should harmonise the inputs to these models too
#       (where here, 'harmonise' means
#       ensure a smooth transition from history to the future).
#       We don't do this step as such a harmonisation would only provide
#       second-order corrections and harmonisation choices
#       for quantities such as gross primary productivity are highly subjective.
# 5. **combine the global-means, seasonality and latitudinal gradients** to
#    produce our complete, gridded projections
# 6. **calculate lower resolution products** from the gridded product
# 7. **release** all the data
#
#
# The updated dataset (produced by adding satellite data information, hereafter referred to as EO-Scenarios) is produced with the same methodology, save from the fact that the historical global means, seasonality, and latitudinal gradients of the historical period are taken from the EO-CMIP7 dataset, described in the user guide: 'Adding Earth Observation data to the generation of CMIP7 CO₂ and CH₄ Concentration Historical Dataset'. In short, it is a concentration dataset produced with the same methods as the GHG CMIP7 concentrations, but with an added observation dataset coprised of satellite observations pre-processed to be consistent with ground based observations.
#

# %% [markdown] editable=true slideshow={"slide_type": ""}
# # Finding and accessing the data

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ## Zenodo
#
# In order to provide more
# reliable, citable access to the data, we also provide it on Zenodo {raw-latex}`\parencite{zenodo}`. The data, as well as all the source code
# and input data used to process it, can be found at
#
# [TODO: add data on zenodo, put link here]
#
# There is no current plan to store data on ESGF, but it can be done in case this dataset is of high importance for enough users

# %% [markdown] editable=true slideshow={"slide_type": ""} jp-MarkdownHeadingCollapsed=true
# # Data description

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ## Format
#
# The data is provided in netCDF format {raw-latex}`\parencite{unidata_netcdf}`.
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
# In EO-Scenarios data, the scenario information follows the same logic as with the baseline CMIP7 Scenarios dataset.
#
# Determining the scenario to which each dataset applies
# is not trivial
# (there was [discussion](https://github.com/PCMDI/input4MIPs_CVs/discussions/64)[^gh-discussion-input4mips-cvs-64]
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
# For example, for the source ID `CR-ml-1-1-0`,
# the scenario identifier is `ml`.
# A Python function for doing this extraction is below.
#
# [^gh-discussion-input4mips-cvs-64]: https://github.com/PCMDI/input4MIPs_CVs/discussions/64


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


print(f"{extract_scenario_id('CR-ml-1-1-0')=}")
print(f"{extract_scenario_id('CR-l-1-1-0')=}")

# %% [markdown] editable=true slideshow={"slide_type": ""}
# These scenario IDs can then be used to find details of the complete scenario.
# These details will be provided both via the CMIP CVs
# (see https://github.com/WCRP-CMIP/CMIP7-CVs)
# and the final ScenarioMIP paper
# (revisions of {raw-latex}`\textcite{van_vuuren_scenariomip_2025}`
#  are expected soon).
# Scenario IDs of the final datasets can be confirmed
# [here](https://github.com/WCRP-CMIP/CMIP7-CVs/discussions/1#discussioncomment-14585785)[^scenario-id-final-comment].
#
# [^scenario-id-final-comment]: https://github.com/WCRP-CMIP/CMIP7-CVs/discussions/1#discussioncomment-14585785

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
# The species provided are:
#
# - CH{raw-latex}`\textsubscript{4}`, in ppb
# - CO{raw-latex}`\textsubscript{2}`, in ppm
#
# This means that, compared with the original CMIP7 forcing dataset, we do not provide data for N{raw-latex}`\textsubscript{2}`O, ozone depleting substances, HCFCs, Halons, nor ozone flourinated compounds.
#
# Adding satellite data to this pipeline for other gases is non-trivial and heavily depends on the ground-based data quality and availability, physics and chemistry of the total column of gas and their impacts on the satellite retrieval algorithms, and on the available satellite data offering.

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
# (using the [gradient-aware-harmonisation](https://github.com/climate-resource/gradient-aware-harmonisation)[^gah]
# package developed as part of the ESA project).
# This leads to smoother, more realistic transitions
# between the history, observation-based period
# and the scenario, model-based projections.
#
# [^gah]: https://github.com/climate-resource/gradient-aware-harmonisation
#
# ### File formats and naming
#
# The file formats are generally close to CMIP6.
# As in CMIP6, we do not provide any vertical profiles.
# For users who require such profiles,
# we refer to the 'The vertical dimension' sub-header
# in Section 4 of {raw-latex}`\cite{meinshausen_historical_2017}`.
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
# }
# ```

# %% [markdown] editable=true slideshow={"slide_type": ""}
# # User guide

# %% [markdown] editable=true slideshow={"slide_type": ""}
# Having downloaded the data, using it is relatively straightforward
# (scenario identification issue discussed above notwithstanding).

# %% [markdown] editable=true slideshow={"slide_type": ""} jp-MarkdownHeadingCollapsed=true
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
# Below we show the filenames for the CO{raw-latex}`\textsubscript{2}` output.

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
# Ensure data is downloaded
query_kwargs_co2_yearly_global = dict(
    ghg="co2",
    time_sampling="yr",
    grid="gm",
    cmip_era="CMIP7",
    source_version="1.1.0",
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


# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
def fetch_and_load_eo(**query_kwargs):
    """
    Load data, using the local EO-updated data for CO2 and CH4

    Both the scenarios and the history are read from the local bundles.
    Only the CMIP6 data comes from ESGF.
    The extension scenarios are not loaded (every caller drops them anyway).
    """
    if (
        query_kwargs.get("target_mip") == "CMIP"
        and query_kwargs.get("cmip_era") == "CMIP7"
        and query_kwargs.get("ghg") in EO_UPDATE_GASES
    ):
        return load_eo_update_history(
            ghg=query_kwargs["ghg"],
            time_sampling=query_kwargs["time_sampling"],
            grid=query_kwargs["grid"],
            with_scenario_dim=True,
        )

    if (
        query_kwargs.get("target_mip") == "ScenarioMIP"
        and query_kwargs.get("cmip_era") == "CMIP7"
        and query_kwargs.get("source_version") == "1.1.0"
        and query_kwargs.get("ghg") in EO_UPDATE_GASES
    ):
        return load_eo_update_scenarios(
            ghg=query_kwargs["ghg"],
            time_sampling=query_kwargs["time_sampling"],
            grid=query_kwargs["grid"],
            include_ext=False,
        )

    return fetch_and_load(**query_kwargs)


# (EO-updated data is read from the local bundle, see EO_UPDATE_DATA_ROOT)

# Get file paths
co2_yearly_global_fps = get_eo_update_local_files(
    ghg="co2", time_sampling="yr", grid="gm"
)

# %% editable=true slideshow={"slide_type": ""} tags=["remove_input"]
for fp in sorted(co2_yearly_global_fps):
    if "ext" in fp.name:
        continue
    print(f"- {fp.name}")

# %% [markdown] editable=true slideshow={"slide_type": ""}
# Output for methane is named identically,
# with `co2` being replaced by `ch4`.
# For methane the filenames are:

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
# Ensure data is downloaded
query_kwargs_ch4_yearly_global = {
    **query_kwargs_co2_yearly_global,
    "ghg": "ch4",
}
# (EO-updated data is read from the local bundle, see EO_UPDATE_DATA_ROOT)

# %% editable=true slideshow={"slide_type": ""} tags=["remove_input"]
# Get file paths
ch4_yearly_global_fps = get_eo_update_local_files(
    ghg="ch4", time_sampling="yr", grid="gm"
)
for fp in sorted(ch4_yearly_global_fps):
    if "ext" in fp.name:
        continue
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
# !ncdump -h {next(v for v in co2_yearly_global_fps if "ext" not in v.name)} | fold -w 80 -s

# %% [markdown] editable=true slideshow={"slide_type": ""}
# Using a tool like [xarray](https://github.com/pydata/xarray)[^xarray-url],
# loading and working with the data is trivial.
# The resulting plot is shown in {numref}`Figure %s <ds-co2-yearly-global-fig>`.
#
# [^xarray-url]: https://github.com/pydata/xarray

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
file_to_open = next(v for v in co2_yearly_global_fps if "ext" not in str(v))

# %% editable=true slideshow={"slide_type": ""}
file_to_open.name

# %% editable=true slideshow={"slide_type": ""}

time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)
ds_example_co2_yearly_global = xr.open_dataset(
    file_to_open,
    decode_times=time_coder,
)

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
# Force values to compute to avoid dask getting involved
ds_example_co2_yearly_global = ds_example_co2_yearly_global.compute()

# %% editable=true slideshow={"slide_type": ""}
ds_example_co2_yearly_global

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
fig, ax = plt.subplots(figsize=(6, 3))
ds_example_co2_yearly_global["co2"].plot.scatter(alpha=0.4)
glue("ds-co2-yearly-global-fig", fig, display=False)
plt.show()

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ```{glue:figure} ds-co2-yearly-global-fig
# ---
# width: 450px
# name: "ds-co2-yearly-global-fig"
# ---
#
# Atmospheric CO{raw-latex}`\textsubscript{2}` concentrations from year 2022 to 2100
# in one of our future forcing datasets.
# ```

# %% [markdown] editable=true slideshow={"slide_type": ""} jp-MarkdownHeadingCollapsed=true
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
# as a scatter or a step plot
# ({numref}`Figure %s <ds_co2_yearly_step_fig>`).
# (The same logic applies to any spatial plots
# which could be created from our datasets
# that include spatial dimensions).

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
ds_plt = ds_example_co2_yearly_global.isel(time=slice(None, 5))

fig, ax = plt.subplots(figsize=(8, 4))
ds_plt["co2"].plot.scatter(ax=ax)

for bounds, val in zip(ds_plt["time_bnds"].values, ds_plt["co2"].values):
    ax.plot(bounds, [val, val], color="tab:blue", linewidth=1.0, alpha=0.7)

xticks = [cftime.DatetimeGregorian(y, 1, 1) for y in range(2022, 2028)]
ax.set_xticks(xticks)
ax.set_xlim(xticks[0], xticks[-1])
ax.grid()

glue("ds_co2_yearly_step_fig", fig, display=False)
plt.show()

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ```{glue:figure} ds_co2_yearly_step_fig
# ---
# width: 450px
# name: "ds_co2_yearly_step_fig"
# ---
#
# Illustration of the fact that each data point
# represents the average over its time bounds, not instantaneous values.
# As a result, it should be plotted with steps or scatters
# rather than an interpolated line.
# ```

# %% [markdown] editable=true slideshow={"slide_type": ""} jp-MarkdownHeadingCollapsed=true
# ## Monthly-, global-mean data
#
# If you want to have information at a finer level
# of temporal detail, we also provide monthly files.
#
# For monthly data, the time labels in the filename are months.
# Below we show the filenames for the CO{raw-latex}`\textsubscript{2}` output.

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
# Ensure data is downloaded
query_kwargs_co2_monthly_global = {
    **query_kwargs_co2_yearly_global,
    "time_sampling": "mon",
}
# (EO-updated data is read from the local bundle, see EO_UPDATE_DATA_ROOT)

# %% editable=true slideshow={"slide_type": ""} tags=["remove_input"]
# Get file paths
co2_monthly_global_fps = get_eo_update_local_files(
    ghg="co2", time_sampling="mon", grid="gm"
)
for fp in sorted(co2_monthly_global_fps):
    if "ext" in fp.name:
        continue

    print(f"- {fp.name}")

# %% [markdown] editable=true slideshow={"slide_type": ""}
# Again, the data can be trivially loaded with [xarray](https://github.com/pydata/xarray).

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
file_to_open = next(v for v in co2_monthly_global_fps if "ext" not in str(v))

# %% editable=true slideshow={"slide_type": ""}
file_to_open.name

# %% editable=true slideshow={"slide_type": ""}
ds_example_co2_monthly_global = xr.open_mfdataset(
    file_to_open,
    decode_times=time_coder,
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
# Scatter or step plots should be used instead
# ({numref}`Figure %s <ds-co2-monthly-global-fig>`).

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
ds_plt = ds_example_co2_monthly_global.isel(time=slice(None, 5 * 12))

fig, ax = plt.subplots(figsize=(8, 4))
ds_plt["co2"].plot.scatter(ax=ax)

for bounds, val in zip(ds_plt["time_bnds"].values, ds_plt["co2"].values):
    ax.plot(bounds, [val, val], color="tab:blue", linewidth=1.0, alpha=0.7)

xticks = [cftime.DatetimeGregorian(y, 1, 1) for y in range(2022, 2028)]
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

xticks = [cftime.DatetimeGregorian(y, 1, 1) for y in range(2022, 2028)]
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
# but the same idea applies to methane.
# ```

# %% [markdown] editable=true slideshow={"slide_type": ""} jp-MarkdownHeadingCollapsed=true
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
# for CO{raw-latex}`\textsubscript{2}`

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
# Ensure data is downloaded
query_kwargs_co2_monthly_lat = {
    **query_kwargs_co2_yearly_global,
    "time_sampling": "mon",
    "grid": "gnz",
}
# (EO-updated data is read from the local bundle, see EO_UPDATE_DATA_ROOT)

# %% editable=true slideshow={"slide_type": ""} tags=["remove_input"]
# Get file paths
co2_monthly_lat_fps = get_eo_update_local_files(
    ghg="co2", time_sampling="mon", grid="gnz"
)
for fp in sorted(co2_monthly_lat_fps):
    if "ext" in fp.name:
        continue

    print(f"- {fp.name}")

# %% [markdown] editable=true slideshow={"slide_type": ""}
# Again, the data can be trivially loaded with [xarray](https://github.com/pydata/xarray).

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
file_to_open = next(v for v in co2_monthly_lat_fps if "ext" not in str(v))

# %% editable=true slideshow={"slide_type": ""}
file_to_open.name

# %% editable=true slideshow={"slide_type": ""}
ds_example_co2_monthly_lat = xr.open_mfdataset(
    file_to_open,
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
# Scatter or step plots should be used instead
# ({numref}`Figure %s <ds-co2-monthly-lat-fig>`).

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
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

xticks = [cftime.DatetimeGregorian(y, 1, 1) for y in range(2022, 2028)]
ax.set_xticks(xticks)
ax.set_xlim(xticks[0], xticks[-1])
ax.grid()
ax.set_title("")

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

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ## Transition from history
#
# Just like for the baseline CMIP7 Scenarios forcings datasets and their respective CMIP7 historical forcings, each dataset of EO-Scenarios transitions smoothly from the historical data updated with EO observations, EO-CMIP7.

# %% [markdown] editable=true slideshow={"slide_type": ""} jp-MarkdownHeadingCollapsed=true
# ### Annual-, global-mean data
#
# We start with the annual- , global-mean data
# ({numref}`Figure %s <ds-co2-transition-from-history-yearly-fig>`).

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
query_kwargs_co2_yearly_global_history = {
    **query_kwargs_co2_yearly_global,
    "cmip_era": "CMIP7",
    "target_mip": "CMIP",
}


def fetch_and_load_history(**query_kwargs):
    """Load the EO-updated historical data (no ESGF fetch needed)"""
    return load_eo_update_history(
        ghg=query_kwargs["ghg"],
        time_sampling=query_kwargs["time_sampling"],
        grid=query_kwargs["grid"],
    )


ds_history_co2_yearly_global = fetch_and_load_history(
    **query_kwargs_co2_yearly_global_history
).compute()

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
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
glue("ds-co2-transition-from-history-yearly-fig", fig, display=False)
plt.show()

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ```{glue:figure} ds-co2-transition-from-history-yearly-fig
# ---
# width: 500px
# name: "ds-co2-transition-from-history-yearly-fig"
# ---
#
# Transition from the historical dataset to the future
# in the global-, annual-mean dataset.
# The transition is smooth, preserving both the absolute value and gradient
# at the transition point (also known as the 'harmonisation point').
# Here this is illustrated with the CO{raw-latex}`\textsubscript{2}` dataset
# for a single scenario, but the same applies to methane and all scenarios.
# ```
#
# {raw-latex}`\newpage`

# %% [markdown] editable=true slideshow={"slide_type": ""} jp-MarkdownHeadingCollapsed=true
# ### Monthly-, global-mean data
#
# Next we consider the monthly, global-mean data
# ({numref}`Figure %s <ds-co2-transition-from-history-monthly-fig>`).

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
query_kwargs_co2_monthly_global_history = {
    **query_kwargs_co2_yearly_global_history,
    "time_sampling": "mon",
}

ds_history_co2_monthly_global = fetch_and_load_history(
    **query_kwargs_co2_monthly_global_history
).compute()

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
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
glue("ds-co2-transition-from-history-monthly-fig", fig, display=False)
plt.show()

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ```{glue:figure} ds-co2-transition-from-history-monthly-fig
# ---
# width: 500px
# name: "ds-co2-transition-from-history-monthly-fig"
# ---
#
# Transition from the historical dataset to the future
# in the global-, monthly-mean dataset.
# The transition is smooth,
# preserving both the absolute value and gradient
# at the transition point (also known as the 'harmonisation point').
# Here this is illustrated with the CO{raw-latex}`\textsubscript{2}` dataset
# for a single scenario.
# The same idea applies to methane and all scenarios.
# ```
#
# {raw-latex}`\newpage`

# %% [markdown] editable=true slideshow={"slide_type": ""} jp-MarkdownHeadingCollapsed=true
# ### Monthly-, 15-degree latitudinally binned data
#
# Next we consider the monthly, latitudinally-resolved data
# ({numref}`Figure %s <ds-co2-transition-from-history-monthly-lat-fig>`).

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
# pdf

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
start_year = 2015
end_year = 2030

fg = sns.relplot(
    data=pdf.loc[:, start_year:end_year].openscm.to_long_data(),
    x="time",
    y="value",
    col="lat",
    col_wrap=3,
    col_order=sorted(pdf.index.get_level_values("lat").unique())[::-1],
    kind="scatter",
    hue="experiment",
    height=2.5,
    aspect=1.0,
)

unit_l = pdf.index.get_level_values("unit").unique()
if len(unit_l) > 1:
    raise AssertionError(unit_l)
unit = unit_l[0]
for ax in fg.axes.flatten():
    if ax.get_ylabel():
        ax.set_ylabel(unit)


glue("ds-co2-transition-from-history-monthly-lat-fig", fg.fig, display=False)
plt.show()

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ```{glue:figure} ds-co2-transition-from-history-monthly-lat-fig
# ---
# width: 600px
# name: "ds-co2-transition-from-history-monthly-lat-fig"
# ---
#
# Transition from the historical dataset to the future
# in the global-, latitudinally-resolved dataset.
# The transition is smooth,
# preserving both the absolute value and gradient
# at the transition point (also known as the 'harmonisation point').
# Here this is illustrated with the CO{raw-latex}`\textsubscript{2}` dataset
# for a single scenario.
# The same idea applies to methane and all scenarios.
# ```
#
# {raw-latex}`\newpage`

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ## Full scenario set
#
# Having seen the transition for a single scenario,
# we now show the full scenario set
# including the transition from history.

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ### Annual-, global-mean data
#
# We start with the annual-, global-mean data
# ({numref}`Figure %s <scenarios-yearly-fig>`).

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
gases_to_show = ["co2", "ch4"]
ds_gases_full_yearly_d = {}
for gas in gases_to_show:
    ds_gases_full_yearly_d[gas] = {}
    for key, target_mip, source_version, institution_id, cmip_era in (
        ("history", "CMIP", "1.0.0", "CR", "CMIP7"),
        ("scenarios", "ScenarioMIP", "1.1.0", "CR", "CMIP7"),
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
        ds = fetch_and_load_eo(**query_kwargs)

        # Unify time axis days to simplify
        ds["time"] = [
            cftime.DatetimeProlepticGregorian(v.year, v.month, 15)
            for v in ds["time"].values
        ]

        # compute to avoid dask weirdness
        ds = ds.compute()

        # avoid extensions
        ds = ds.sel(
            scenario=[v for v in ds["scenario"].values if "ext" not in v]
        ).dropna("time")

        ds_gases_full_yearly_d[gas][key] = ds

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
palette = {
    "history": "k",
    "history-cmip6": "tab:grey",
    "historical": "k",
    "historical-cmip6": "tab:grey",
    "vl": "#24a4ff",
    "ln": "#4a0daf",
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
    "ln",
    "vl",
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

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
col = "ghg"

fg = sns.relplot(
    data=pdf.openscm.to_long_data(),
    x="time",
    y="value",
    hue="experiment",
    palette={k: v for k, v in palette.items() if k in pdf.pix.unique("experiment")},
    hue_order=[v for v in hue_order if v in pdf.pix.unique("experiment")],
    kind="scatter",
    col=col,
    col_wrap=2,
    facet_kws=dict(sharey=False),
    height=3.5,
    aspect=1.1,
)

for ax in fg.axes.flatten():
    col_val = ax.get_title().split(f"{col} = ")[1]
    unit_l = (
        pdf.loc[pix.ismatch(**{col: col_val})].index.get_level_values("unit").unique()
    )
    if len(unit_l) > 1:
        raise AssertionError(unit_l)
    unit = unit_l[0]
    ax.set_ylabel(unit)

glue("scenarios-yearly-fig", fg.fig, display=False)
plt.show()

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ```{glue:figure} scenarios-yearly-fig
# ---
# width: 600px
# name: "scenarios-yearly-fig"
# ---
#
# Scenarios including the
# transition from the historical period to the future/projection period
# in the global-, annual-mean datasets.
# ```
#
# {raw-latex}`\newpage`

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ### Monthly-, global-mean data
#
# Then we consider the data including seasonality
# i.e. the monthly-, global-mean data
# ({numref}`Figure %s <scenarios-monthly-fig>`).

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
ds_gases_full_monthly_d = {}
for gas in gases_to_show:
    ds_gases_full_monthly_d[gas] = {}
    for key, target_mip, source_version, institution_id, cmip_era in (
        ("history", "CMIP", "1.0.0", "CR", "CMIP7"),
        ("scenarios", "ScenarioMIP", "1.1.0", "CR", "CMIP7"),
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
        ds = fetch_and_load_eo(**query_kwargs)

        # Unify time axis days to simplify
        ds["time"] = [
            cftime.DatetimeProlepticGregorian(v.year, v.month, 15)
            for v in ds["time"].values
        ]

        # compute to avoid dask weirdness
        ds = ds.compute()

        # avoid extensions
        ds = ds.sel(
            scenario=[v for v in ds["scenario"].values if "ext" not in v]
        ).dropna("time")

        ds_gases_full_monthly_d[gas][key] = ds

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
start_year = 2020
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

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
col = "ghg"

fg = sns.relplot(
    data=pdf.openscm.to_long_data(),
    x="time",
    y="value",
    hue="experiment",
    palette={k: v for k, v in palette.items() if k in pdf.pix.unique("experiment")},
    hue_order=[v for v in hue_order if v in pdf.pix.unique("experiment")],
    kind="scatter",
    col=col,
    col_wrap=2,
    facet_kws=dict(sharey=False),
    height=3.5,
    aspect=1.1,
)

for ax in fg.axes.flatten():
    col_val = ax.get_title().split(f"{col} = ")[1]
    unit_l = (
        pdf.loc[pix.ismatch(**{col: col_val})].index.get_level_values("unit").unique()
    )
    if len(unit_l) > 1:
        raise AssertionError(unit_l)
    unit = unit_l[0]
    ax.set_ylabel(unit)

# plt.tight_layout()
glue("scenarios-monthly-fig", fg.fig, display=False)
plt.show()

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ```{glue:figure} scenarios-monthly-fig
# ---
# width: 600px
# name: "scenarios-monthly-fig"
# ---
#
# Scenarios including the
# transition from the historical period to the future/projection period
# in the monthly-, annual-mean datasets.
# ```
#
# {raw-latex}`\newpage`

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ### Monthly-, 15-degree latitudinal data
#
# Then we consider the data including seasonality
# and latitudinal-resolution
# i.e. the monthly-, latitudinally-resolved data
# ({numref}`Figure %s <scenarios-monthly-lat-fig>`).

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
ds_gases_full_monthly_lat_d = {}
for gas in gases_to_show:
    ds_gases_full_monthly_lat_d[gas] = {}
    for key, target_mip, source_version, institution_id, cmip_era in (
        ("history", "CMIP", "1.0.0", "CR", "CMIP7"),
        ("scenarios", "ScenarioMIP", "1.1.0", "CR", "CMIP7"),
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
        ds = fetch_and_load_eo(**query_kwargs)

        # Unify time axis days to simplify
        ds["time"] = [
            cftime.DatetimeProlepticGregorian(v.year, v.month, 15)
            for v in ds["time"].values
        ]

        # compute to avoid dask weirdness
        ds = ds.compute()

        # avoid extensions
        ds = ds.sel(
            scenario=[v for v in ds["scenario"].values if "ext" not in v]
        ).dropna("time")

        ds_gases_full_monthly_lat_d[gas][key] = ds

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
start_year = 2020
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

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
col = "lat"
row = "ghg"

fg = sns.relplot(
    data=pdf.loc[pix.isin(lat=[-82.5, 7.5, 82.5])].openscm.to_long_data(),
    x="time",
    y="value",
    hue="experiment",
    palette={k: v for k, v in palette.items() if k in pdf.pix.unique("experiment")},
    hue_order=[v for v in hue_order if v in pdf.pix.unique("experiment")],
    kind="scatter",
    row=row,
    col=col,
    facet_kws=dict(sharey=False),
    height=2.5,
    aspect=1.0,
)

for ax in fg.axes.flatten():
    if ax.get_ylabel():
        t_split = ax.get_title().split("|")
        row_val = t_split[0].split(f"{row} = ")[1].strip()

        unit_l = (
            pdf.loc[(pdf.index.get_level_values(row) == row_val)]
            .index.get_level_values("unit")
            .unique()
        )
        if len(unit_l) > 1:
            raise AssertionError(unit_l)
        unit = unit_l[0]
        ax.set_ylabel(unit)

# plt.tight_layout()
glue("scenarios-monthly-lat-fig", fg.fig, display=False)
plt.show()

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ```{glue:figure} scenarios-monthly-lat-fig
# ---
# width: 500px
# name: "scenarios-monthly-lat-fig"
# ---
#
# Scenarios including the
# transition from the historical period to the future/projection period
# in the monthly-, latitudinally-resolved datasets.
# ```
#
# {raw-latex}`\newpage`

# %% [markdown]
# One thing to note is that for methane, when looking at the high latitude plot, we see a mismatch of seasonality oscillation between the historical data with the EO-update and the scenario projections. This is weird and we should look into it [TODO: look into it]
#

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ## Differences from the baseline CMIP7 Scenario datasets

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ### Data comparisons
#
# Comparing the data from EO-Scenarios and the baseline CMIP7 Scenarios shows the following changes:
#
# 1. Most obviously, the onlu gases provided here are CO{raw-latex}`\textsubscript{2}` and CH{raw-latex}`\textsubscript{4}`
# [TODO: add analysis]
#

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ### Annual-, global-mean data
#
# We start with the annual-, global-mean data comparison of CO{raw-latex}`\textsubscript{2}` and CH{raw-latex}`\textsubscript{4}` for all seven scenarios.
#
# ({numref}`Figure %s <scenarios-cmip-comparison-yearly-co2-fig>` for CO{raw-latex}`\textsubscript{2}`
# and {numref}`Figure %s <scenarios-cmip-comparison-yearly-ch4-fig>` for CH{raw-latex}`\textsubscript{4}`).

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
    "vl": "low",
    "ln": "low",
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
    "ln",
    "vl",
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
        ("scenarios-cmip7", "ScenarioMIP", "1.1.0", "CR", "CMIP7"),
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
        ds = fetch_and_load_eo(**query_kwargs)

        # Unify time axis days to simplify
        ds["time"] = [
            cftime.DatetimeProlepticGregorian(v.year, v.month, 15)
            for v in ds["time"].values
        ]

        # compute to avoid dask weirdness
        ds = ds.compute()

        # avoid extensions
        ds = ds.sel(
            scenario=[v for v in ds["scenario"].values if "ext" not in v]
        ).dropna("time")

        ds_gases_full_yearly_multi_phase_d[gas][key] = ds

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

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
# Original (baseline) CMIP7 data, to compare against the EO-updated data.
# `fetch_and_load` is the plain ESGF loader (no local EO-update routing).
ds_baseline_cmip7_yearly_d = {}
for gas in gases_to_show:
    ds_baseline_cmip7_yearly_d[gas] = {}
    for key, target_mip, source_version in (
        ("history", "CMIP", "1.0.0"),
        ("scenarios", "ScenarioMIP", "1.1.0"),
    ):
        ds = fetch_and_load(
            ghg=gas,
            time_sampling="yr",
            grid="gm",
            target_mip=target_mip,
            source_version=source_version,
            institution_id="CR",
            cmip_era="CMIP7",
            engine=engine,
        )

        # Unify time axis days to simplify
        ds["time"] = [
            cftime.DatetimeProlepticGregorian(v.year, v.month, 15)
            for v in ds["time"].values
        ]

        # compute to avoid dask weirdness
        ds = ds.compute()

        # avoid extensions
        ds = ds.sel(
            scenario=[v for v in ds["scenario"].values if "ext" not in v]
        ).dropna("time")

        ds_baseline_cmip7_yearly_d[gas][key] = ds


# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
def get_history_plus_scenario_series(  # noqa: PLR0913
    ds_history: xr.Dataset,
    ds_scenarios: xr.Dataset,
    gas: str,
    scenario: str,
    start_year: int,
    end_year: int,
    monthly: bool,
) -> pd.Series:
    """
    Get the history followed by the scenario as one series

    The index is the year (annual data) or the fractional year (monthly data).
    The scenario takes over from the history in the first time step of the scenario
    (the two overlap in 2022).
    """

    def get_index(da: xr.DataArray) -> np.ndarray:
        if monthly:
            return (da["time"].dt.year + (da["time"].dt.month - 0.5) / 12).values

        return da["time"].dt.year.values

    history = ds_history[gas].isel(scenario=0).dropna("time")
    scen = ds_scenarios[gas].sel(scenario=scenario).dropna("time")

    history_s = pd.Series(history.values, index=get_index(history))
    scen_s = pd.Series(scen.values, index=get_index(scen))

    history_s = history_s[history_s.index < scen_s.index.min()]
    res = pd.concat([history_s, scen_s]).sort_index()

    return res[(res.index >= start_year) & (res.index <= end_year + 1)]


def plot_eo_vs_baseline(  # noqa: PLR0913
    gas: str,
    ds_eo_d: dict[str, dict[str, xr.Dataset]],
    ds_baseline_d: dict[str, dict[str, xr.Dataset]],
    start_year: int,
    end_year: int,
    monthly: bool,
) -> matplotlib.figure.Figure:
    """
    Plot EO-updated vs baseline CMIP7 data (left) and their difference (right)

    One row per scenario.
    The left-hand y-axes are free, the right-hand (difference) y-axes are shared.
    """
    scenarios = [
        v for v in hue_order if v in ds_baseline_d[gas]["scenarios"]["scenario"].values
    ]
    unit = ds_baseline_d[gas]["scenarios"][gas].attrs["units"]

    fig, axes = plt.subplots(
        len(scenarios),
        2,
        figsize=(8.0, 1.35 * len(scenarios) + 0.8),
        sharex=True,
        squeeze=False,
        gridspec_kw=dict(width_ratios=[1.4, 1.0]),
    )
    for ax in axes[1:, 1]:
        ax.sharey(axes[0, 1])

    for (ax_levels, ax_diff), scenario in zip(axes, scenarios):
        series_kwargs = dict(
            gas=gas,
            scenario=scenario,
            start_year=start_year,
            end_year=end_year,
            monthly=monthly,
        )
        eo = get_history_plus_scenario_series(
            ds_eo_d[gas]["history-cmip7"],
            ds_eo_d[gas]["scenarios-cmip7"],
            **series_kwargs,
        )
        baseline = get_history_plus_scenario_series(
            ds_baseline_d[gas]["history"],
            ds_baseline_d[gas]["scenarios"],
            **series_kwargs,
        )
        colour = palette[scenario]

        ax_levels.plot(
            baseline.index, baseline.values, color=colour, lw=3.0, alpha=0.35
        )
        ax_levels.plot(eo.index, eo.values, color=colour, lw=0, marker="o", ms=2.5)
        ax_levels.set_ylabel(unit)
        ax_levels.set_title(f"scenario = {scenario}", loc="left", fontsize=9)

        diff = (eo - baseline).dropna()
        ax_diff.axhline(0.0, color="0.5", lw=0.8)
        ax_diff.plot(diff.index, diff.values, color=colour, lw=1.5)
        ax_diff.set_title("EO-update minus baseline", loc="left", fontsize=9)
        ax_diff.set_ylabel(unit)

    axes[-1, 0].set_xlabel("time")
    axes[-1, 1].set_xlabel("time")

    legend_handles = [
        matplotlib.lines.Line2D(
            [], [], color="0.3", lw=3.0, alpha=0.35, label="baseline CMIP7"
        ),
        matplotlib.lines.Line2D(
            [], [], color="0.3", lw=0, marker="o", ms=3.0, label="EO-update"
        ),
    ]
    fig.legend(handles=legend_handles, loc="upper center", ncol=2, frameon=False)
    fig.tight_layout(rect=(0, 0, 1, 0.97))

    return fig


# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
fig = plot_eo_vs_baseline(
    "co2",
    ds_eo_d=ds_gases_full_yearly_multi_phase_d,
    ds_baseline_d=ds_baseline_cmip7_yearly_d,
    start_year=2000,
    end_year=2100,
    monthly=False,
)
glue("scenarios-cmip-comparison-yearly-co2-fig", fig, display=False)
plt.show()

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
fig = plot_eo_vs_baseline(
    "ch4",
    ds_eo_d=ds_gases_full_yearly_multi_phase_d,
    ds_baseline_d=ds_baseline_cmip7_yearly_d,
    start_year=2000,
    end_year=2100,
    monthly=False,
)
glue("scenarios-cmip-comparison-yearly-ch4-fig", fig, display=False)
plt.show()

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ```{glue:figure} scenarios-cmip-comparison-yearly-co2-fig
# ---
# width: 560px
# name: "scenarios-cmip-comparison-yearly-co2-fig"
# ---
#
# Comparison of the EO-updated and the baseline CMIP7 history and scenarios
# in the global-, annual-mean CO{raw-latex}`\textsubscript{2}` datasets.
# Left: the EO-updated data (markers) and the baseline CMIP7 data (line),
# one row per scenario (free y-axes).
# Right: the difference between the two (EO-updated minus baseline),
# on a common y-axis across scenarios.
# ```
#
# ```{glue:figure} scenarios-cmip-comparison-yearly-ch4-fig
# ---
# width: 560px
# name: "scenarios-cmip-comparison-yearly-ch4-fig"
# ---
#
# As {numref}`Figure %s <scenarios-cmip-comparison-yearly-co2-fig>`,
# but for CH{raw-latex}`\textsubscript{4}`.
# ```
#
# {raw-latex}`\newpage`

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ### Monthly-, global-mean data
#
# Finally, we consider the monthly-, global-mean data
# ({numref}`Figure %s <scenarios-cmip-comparison-monthly-co2-fig>` for CO{raw-latex}`\textsubscript{2}`
# and {numref}`Figure %s <scenarios-cmip-comparison-monthly-ch4-fig>` for CH{raw-latex}`\textsubscript{4}`).

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
ds_gases_full_monthly_multi_phase_d = {}
for gas in gases_to_show:
    ds_gases_full_monthly_multi_phase_d[gas] = {}
    for key, target_mip, source_version, institution_id, cmip_era in (
        ("history-cmip7", "CMIP", "1.0.0", "CR", "CMIP7"),
        ("scenarios-cmip7", "ScenarioMIP", "1.1.0", "CR", "CMIP7"),
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
        ds = fetch_and_load_eo(**query_kwargs)

        # Unify time axis days to simplify
        ds["time"] = [
            cftime.DatetimeProlepticGregorian(v.year, v.month, 15)
            for v in ds["time"].values
        ]

        # compute to avoid dask weirdness
        ds = ds.compute()

        # avoid extensions
        ds = ds.sel(
            scenario=[v for v in ds["scenario"].values if "ext" not in v]
        ).dropna("time")

        ds_gases_full_monthly_multi_phase_d[gas][key] = ds

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

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
# Original (baseline) CMIP7 monthly data, to compare against the EO-updated data.
ds_baseline_cmip7_monthly_d = {}
for gas in gases_to_show:
    ds_baseline_cmip7_monthly_d[gas] = {}
    for key, target_mip, source_version in (
        ("history", "CMIP", "1.0.0"),
        ("scenarios", "ScenarioMIP", "1.1.0"),
    ):
        ds = fetch_and_load(
            ghg=gas,
            time_sampling="mon",
            grid="gm",
            target_mip=target_mip,
            source_version=source_version,
            institution_id="CR",
            cmip_era="CMIP7",
            engine=engine,
        )

        # Unify time axis days to simplify
        ds["time"] = [
            cftime.DatetimeProlepticGregorian(v.year, v.month, 15)
            for v in ds["time"].values
        ]

        # compute to avoid dask weirdness
        ds = ds.compute()

        # avoid extensions
        ds = ds.sel(
            scenario=[v for v in ds["scenario"].values if "ext" not in v]
        ).dropna("time")

        ds_baseline_cmip7_monthly_d[gas][key] = ds

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
fig = plot_eo_vs_baseline(
    "co2",
    ds_eo_d=ds_gases_full_monthly_multi_phase_d,
    ds_baseline_d=ds_baseline_cmip7_monthly_d,
    start_year=2012,
    end_year=2025,
    monthly=True,
)
glue("scenarios-cmip-comparison-monthly-co2-fig", fig, display=False)
plt.show()

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
fig = plot_eo_vs_baseline(
    "ch4",
    ds_eo_d=ds_gases_full_monthly_multi_phase_d,
    ds_baseline_d=ds_baseline_cmip7_monthly_d,
    start_year=2012,
    end_year=2025,
    monthly=True,
)
glue("scenarios-cmip-comparison-monthly-ch4-fig", fig, display=False)
plt.show()

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ```{glue:figure} scenarios-cmip-comparison-monthly-co2-fig
# ---
# width: 560px
# name: "scenarios-cmip-comparison-monthly-co2-fig"
# ---
#
# Comparison of the EO-updated and the baseline CMIP7 history and scenarios
# in the global-, monthly-mean CO{raw-latex}`\textsubscript{2}` datasets
# over the time period from 2012 to 2025.
# Left: the EO-updated data (markers) and the baseline CMIP7 data (line),
# one row per scenario (free y-axes).
# Right: the difference between the two (EO-updated minus baseline),
# on a common y-axis across scenarios.
# ```
#
# ```{glue:figure} scenarios-cmip-comparison-monthly-ch4-fig
# ---
# width: 560px
# name: "scenarios-cmip-comparison-monthly-ch4-fig"
# ---
#
# As {numref}`Figure %s <scenarios-cmip-comparison-monthly-co2-fig>`,
# but for CH{raw-latex}`\textsubscript{4}`.
# ```

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ### Monthly-, 15-degree latitudinally binned data
#
# Finally, we consider the latitudinal structure of the differences,
# using the monthly-, 15-degree latitudinally binned data
# for the `m` scenario
# ({numref}`Figure %s <scenarios-cmip-comparison-monthly-lat-co2-fig>` for CO{raw-latex}`\textsubscript{2}`
# and {numref}`Figure %s <scenarios-cmip-comparison-monthly-lat-ch4-fig>` for CH{raw-latex}`\textsubscript{4}`),
# followed by the same data averaged to annual means
# ({numref}`Figure %s <scenarios-cmip-comparison-yearly-lat-co2-fig>` for CO{raw-latex}`\textsubscript{2}`
# and {numref}`Figure %s <scenarios-cmip-comparison-yearly-lat-ch4-fig>` for CH{raw-latex}`\textsubscript{4}`).

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
# Monthly, 15-degree latitudinally binned data for the `m` scenario:
# EO-updated (local bundle) and baseline (ESGF)
ds_m_monthly_lat_d = {}
for gas in gases_to_show:
    ds_m_monthly_lat_d[gas] = {}
    for key, ds in (
        ("eo", load_eo_update_scenarios(gas, "mon", "gnz", include_ext=False)),
        (
            "baseline",
            fetch_and_load(
                ghg=gas,
                time_sampling="mon",
                grid="gnz",
                target_mip="ScenarioMIP",
                source_version="1.1.0",
                institution_id="CR",
                cmip_era="CMIP7",
                engine=engine,
            ),
        ),
    ):
        # Unify time axis days to simplify
        ds["time"] = [
            cftime.DatetimeProlepticGregorian(v.year, v.month, 15)
            for v in ds["time"].values
        ]

        ds_m_monthly_lat_d[gas][key] = (
            ds.sel(scenario="m")[gas].dropna("time", how="all").compute()
        )


def plot_lat_diff(gas: str, yearly: bool) -> matplotlib.figure.Figure:
    """
    Plot the EO-updated minus baseline difference as a function of time and latitude

    If `yearly`, the monthly data is first averaged to annual means.
    """
    eo_da = ds_m_monthly_lat_d[gas]["eo"]
    baseline_da = ds_m_monthly_lat_d[gas]["baseline"]
    if yearly:
        eo_da = eo_da.resample(time="YS").mean()
        baseline_da = baseline_da.resample(time="YS").mean()

    diff_da = eo_da - baseline_da
    unit = ds_m_monthly_lat_d[gas]["baseline"].attrs.get("units", gas)

    fig, ax = plt.subplots(figsize=(8, 3.5))
    diff_da.plot(
        ax=ax,
        x="time",
        y="lat",
        cmap="RdBu_r",
        cbar_kwargs={"label": f"EO-update minus baseline ({unit})"},
    )
    freq_label = "yearly" if yearly else "monthly"
    ax.set_title(
        f"{gas.upper()} - scenario m - {freq_label} - EO-update minus baseline"
    )
    ax.set_xlabel("Time")
    plt.tight_layout()

    return fig


# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
fig = plot_lat_diff("co2", yearly=False)
glue("scenarios-cmip-comparison-monthly-lat-co2-fig", fig, display=False)
plt.show()

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
fig = plot_lat_diff("ch4", yearly=False)
glue("scenarios-cmip-comparison-monthly-lat-ch4-fig", fig, display=False)
plt.show()

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
fig = plot_lat_diff("co2", yearly=True)
glue("scenarios-cmip-comparison-yearly-lat-co2-fig", fig, display=False)
plt.show()

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
fig = plot_lat_diff("ch4", yearly=True)
glue("scenarios-cmip-comparison-yearly-lat-ch4-fig", fig, display=False)
plt.show()

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ```{glue:figure} scenarios-cmip-comparison-monthly-lat-co2-fig
# ---
# width: 560px
# name: "scenarios-cmip-comparison-monthly-lat-co2-fig"
# ---
#
# Difference between the EO-updated and the baseline CMIP7 data
# (EO-updated minus baseline)
# in the monthly-, 15-degree latitudinally binned CO{raw-latex}`\textsubscript{2}` data
# for the `m` scenario.
# ```
#
# ```{glue:figure} scenarios-cmip-comparison-monthly-lat-ch4-fig
# ---
# width: 560px
# name: "scenarios-cmip-comparison-monthly-lat-ch4-fig"
# ---
#
# As {numref}`Figure %s <scenarios-cmip-comparison-monthly-lat-co2-fig>`,
# but for CH{raw-latex}`\textsubscript{4}`.
# ```
#
# ```{glue:figure} scenarios-cmip-comparison-yearly-lat-co2-fig
# ---
# width: 560px
# name: "scenarios-cmip-comparison-yearly-lat-co2-fig"
# ---
#
# As {numref}`Figure %s <scenarios-cmip-comparison-monthly-lat-co2-fig>`,
# but averaged to annual means.
# ```
#
# ```{glue:figure} scenarios-cmip-comparison-yearly-lat-ch4-fig
# ---
# width: 560px
# name: "scenarios-cmip-comparison-yearly-lat-ch4-fig"
# ---
#
# As {numref}`Figure %s <scenarios-cmip-comparison-monthly-lat-ch4-fig>`,
# but averaged to annual means.
# ```
