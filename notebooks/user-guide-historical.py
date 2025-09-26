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
from local.data_loading import fetch_and_load_ghg_dataset

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
# Please see these tools docs for usage instructions.

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

# %% [markdown]
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
#
# Having downloaded the data, using it is quite straightforward.

# %% [markdown] editable=true slideshow={"slide_type": ""}
# 1. ncdump a file for each grid-frequency combo to show dimensions and axes
# 2. load a file using xarray
# 3. load a dataset using xarray
# 4. plot
# 5. comparisons with CMIP6 (just do global- and hemispheric-mean)

# %%
fetch_and_load_ghg_dataset(
    ghg="co2", grid="gm", frequency="yr", cmip_era="CMIP7", source_id="CR-CMIP-1-0-0"
)

# %% [markdown]
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

# %%
