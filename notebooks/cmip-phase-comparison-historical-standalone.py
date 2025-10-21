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
# Here we make some standalone plots
# that compare the historical concentrations over CMIP phases.

# %% [markdown] editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
# ## Imports

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
from functools import partial

import cftime
import matplotlib
import matplotlib.pyplot as plt
import nc_time_axis  # noqa: F401
import numpy as np
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

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ### Data comparisons
#
# Comparing the data from CMIP6 and CMIP7 shows minor changes
# (although doing this comparison requires a bit of care
# because of the changes in file formats).

# %%
fetch_and_load = partial(
    fetch_and_load_ghg_dataset,
    local_data_root_dir=local_data_root_dir,
    # index_node=KnownIndexNode.DKRZ,
    # cmip_era="CMIP6",
    # source_id="UoM-CMIP-1-2-0",
    index_node=KnownIndexNode.ORNL,
)

# %% editable=true slideshow={"slide_type": ""} tags=["remove_cell"]
gases_to_show = ["co2", "ch4", "n2o", "cfc12eq", "hfc134aeq"]
ds_gases_full_d = {}
for gas in gases_to_show:
    ds_gases_full_d[gas] = {}
    for source_id, cmip_era in (
        ("CR-CMIP-1-0-0", "CMIP7"),
        ("UoM-CMIP-1-2-0", "CMIP6"),
        (None, "CMIP5"),
    ):
        query_kwargs = {
            "ghg": gas,
            "time_sampling": "yr",
            "grid": "gm",
            "target_mip": "CMIP",
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

# %%
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
from typing import Callable

import numpy.typing as npt
import xarray as xr


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
        if "_delta" in ax_name:
            continue

        gas = ax_name

        target_unit_conc = ds_d[gas]["CMIP7"][gas].attrs["units"]
        target_unit_re = "W / m^2"

        for cmip_era, ds in ds_d[gas].items():
            label = f"{cmip_era} ({ds.attrs['source_id']})"
            tmp = ds[gas].copy()
            tmp.values = Q(tmp.values, tmp.attrs["units"]).to(target_unit_conc).m
            ds[gas].plot.scatter(
                ax=axes_d[gas], label=label, alpha=0.7, edgecolors="none"
            )

        ax.legend()
        ax.set_title(gas)
        ax.xaxis.set_tick_params(labelbottom=False)
        ax.set_ylabel(target_unit_conc)
        ax.set_xlabel(None)

        ax_delta = axes_d[f"{gas}_delta"]

        da_cmip7 = ds_d[gas]["CMIP7"][gas]

        for cmip_era, ds in ds_d[gas].items():
            if cmip_era == "CMIP7":
                continue

            da_other = ds_d[gas][cmip_era][gas]
            overlapping_times = np.intersect1d(da_other["time"], da_cmip7["time"])

            da_cmip7_st = da_cmip7.sel(time=overlapping_times)
            da_other_st = da_other.sel(time=overlapping_times)

            delta = da_cmip7_st.copy()
            tmp = Q(da_cmip7_st.values, da_cmip7_st.attrs["units"]) - Q(
                da_other_st.values, da_other_st.attrs["units"]
            )
            delta.values = tmp.to(target_unit_conc).m

            delta.plot.scatter(
                ax=ax_delta,
                label=f"CMIP7 - {cmip_era}",
                edgecolors="none",
                s=10,
            )
            ax_delta.axhline(0.0, color="k", linestyle="--")
            ax_delta.legend()

            ax_delta.xaxis.set_tick_params(labelbottom=False)
            ax_delta.set_ylabel(target_unit_conc)
            ax_delta.set_xlabel(None)

            ax_delta_re = axes_d[f"{gas}_delta_re"]

            tmp = RADIATIVE_EFFICIENCIES[gas] * Q(delta.values, delta.attrs["units"])
            delta_re = delta.copy()
            target_unit = "W / m^2"
            delta_re.values = tmp.to(target_unit_re).m
            delta_re.attrs["units"] = target_unit_re

            delta_re.plot.scatter(
                ax=ax_delta_re,
                label=f"CMIP7 - {cmip_era}",
                edgecolors="none",
                s=10,
            )
            ax_delta_re.axhline(0.0, color="k", linestyle="--")

            ax_delta_re.xaxis.set_tick_params(labelbottom=True)
            ax_delta_re.set_ylabel(target_unit_re)
            ax_delta_re.legend()


plt_mosaic = [
    ["co2", "ch4", "n2o"],
    ["co2", "ch4", "n2o"],
    ["co2_delta", "ch4_delta", "n2o_delta"],
    ["co2_delta_re", "ch4_delta_re", "n2o_delta_re"],
    ["cfc12eq", "hfc134aeq", ""],
    ["cfc12eq", "hfc134aeq", ""],
    ["cfc12eq_delta", "hfc134aeq_delta", ""],
    ["cfc12eq_delta_re", "hfc134aeq_delta_re", ""],
]
get_default_delta_mosaic = partial(
    plt.subplot_mosaic,
    mosaic=plt_mosaic,
    figsize=(16, 14),
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

# %%

# %%
# 1. add CMIP5 ingestion and loading
# 2. load yearly, global-mean, make comparisons across CMIP5, CMIP6, CMIP7
#    - plots on same axes
#    - difference plots absolute
#    - difference plots in ERF terms
# 3. load monthly, global-mean, make comparisons across CMIP6 and CMIP7
#    - plots on same axes
#    - difference plots absolute
#    - difference plots in ERF terms
# 4. load monthly, latitudinal, make comparisons across CMIP6 and CMIP7
#    - plots on same axes - faceted by latitude
#    - difference plots absolute
#    - difference plots in ERF terms
assert False, "up to here"

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
# min_year = 1990
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

# %% [markdown]
# Values below come from Table 7.SM.7 of
# IPCC AR7 WG1 Ch. 7 Supplementary Material[^4].
#
# [^4]: https://www.ipcc.ch/report/ar6/wg1/downloads/report/IPCC_AR6_WGI_Chapter07_SM.pdf

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
# (see Section 7.3.5.2 of AR6 WG1 Chapter 7[^5]),
# estimated to be 3.84 W / m<sup>2</sup>
# (very likely range of 3.46 to 4.22 W / m<sup>2</sup>),
# such differences are particularly small.
#
# [^5]: https://www.ipcc.ch/report/ar6/wg1/chapter/chapter-7/

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

# %% [markdown]
# Like the annual-means,
# the atmospheric concentrations including seasonality
# are reasonably consistent between CMIP6 and CMIP7.
# There are some areas of change.
# Full details of these changes will be provided
# in the forthcoming manuscripts.
