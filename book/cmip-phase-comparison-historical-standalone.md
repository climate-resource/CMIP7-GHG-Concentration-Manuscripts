---
authors:
- name: Zebedee Nicholls
- name: Florence Bockting
- name: Mika Plf{\"u}ger
jupytext:
  notebook_metadata_filter: title,authors
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.17.3
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
title: 'CMIP Greenhouse Gas (GHG) Concentration Historical Dataset:

  Data Description and User Guide'
---

+++ {"editable": true, "slideshow": {"slide_type": ""}}

# Comparison of CMIP Phases
## Overview

Here we make some standalone plots
that compare the historical concentrations over CMIP phases.

+++ {"editable": true, "slideshow": {"slide_type": ""}, "tags": ["remove_cell"]}

## Imports

```{code-cell}
---
editable: true
slideshow:
  slide_type: ''
tags: [remove_cell]
---
from functools import partial

import cftime
import matplotlib
import matplotlib.pyplot as plt
import nc_time_axis  # noqa: F401
import numpy as np
import tqdm.auto

from local.data_loading import fetch_and_load_ghg_dataset
from local.esgf.db_helpers import create_all_tables, get_sqlite_engine
from local.esgf.search.search_query import KnownIndexNode
from local.paths import REPO_ROOT
```

```{code-cell}
---
editable: true
slideshow:
  slide_type: ''
tags: [remove_cell]
---
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
```

+++ {"editable": true, "slideshow": {"slide_type": ""}}

## Data comparisons

Comparing the data from CMIP6 and CMIP7 shows minor changes
(although doing this comparison requires a bit of care
because of the changes in file formats).

```{code-cell}
fetch_and_load = partial(
    fetch_and_load_ghg_dataset,
    local_data_root_dir=local_data_root_dir,
    # index_node=KnownIndexNode.DKRZ,
    # cmip_era="CMIP6",
    # source_id="UoM-CMIP-1-2-0",
    index_node=KnownIndexNode.ORNL,
)
```

```{code-cell}
---
editable: true
slideshow:
  slide_type: ''
tags: [remove_cell]
---
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
```

Values below come from Table 7.SM.7 of
IPCC AR7 WG1 Ch. 7 Supplementary Material[^4].

[^4]: https://www.ipcc.ch/report/ar6/wg1/downloads/report/IPCC_AR6_WGI_Chapter07_SM.pdf

```{code-cell}
from openscm_units import unit_registry

Q = unit_registry.Quantity

RADIATIVE_EFFICIENCIES = {
    "co2": Q(1.33e-5, "W / m^2 / ppb"),
    "ch4": Q(3.88e-4, "W / m^2 / ppb"),
    "n2o": Q(3.2e-3, "W / m^2 / ppb"),
    "cfc12eq": Q(0.358, "W / m^2 / ppb"),
    "hfc134aeq": Q(0.167, "W / m^2 / ppb"),
}
```

```{code-cell}
---
editable: true
slideshow:
  slide_type: ''
tags: [remove_cell]
---
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


def sel_lat(
    ds_d: dict[str, dict[str, xr.Dataset]],
    sel_func: Callable[[xr.DataArray], npt.NDArray[bool]],
) -> dict[str, dict[str, xr.Dataset]]:
    """
    Select lat from our dictionary of [xr.Dataset][]'s
    """
    res = {
        gas: {key: value.sel(lat=sel_func(value["lat"])) for key, value in tmp.items()}
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
    # ["co2", "ch4", "n2o"],
    ["co2_delta", "ch4_delta", "n2o_delta"],
    ["co2_delta_re", "ch4_delta_re", "n2o_delta_re"],
    ["cfc12eq", "hfc134aeq", ""],
    # ["cfc12eq", "hfc134aeq", ""],
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
```

+++ {"editable": true, "slideshow": {"slide_type": ""}}

### Global, annual-mean concentrations: Year 1 - 2022

```{code-cell}
plt.rcParams["axes.xmargin"] = 0
```

```{code-cell}
---
editable: true
slideshow:
  slide_type: ''
tags: [remove_input]
---
fig, axes_d = get_default_delta_mosaic()
axes_d = remove_empty_axes(axes_d)

plot_overview_and_deltas(
    ds_gases_full_d,
    axes_d,
)

for ax in axes_d.values():
    xticks = [
        cftime.DatetimeProlepticGregorian(y, 1, 1) for y in np.arange(0, 2000, 500)
    ]
    ax.set_xticks(xticks)
    ax.set_xticklabels([v.year for v in xticks])

plt.tight_layout()
# plt.savefig("key-species-global-annual-changes-across-cmip-phases.png")
plt.show()
```

+++ {"editable": true, "slideshow": {"slide_type": ""}}

### Global, annual-mean concentrations: Year 1750 - 2022

```{code-cell}
# TODO: copy https://github.com/climate-resource/CMIP6-vs-CMIP7-GHG-Concentrations/blob/clean-up/notebooks/0101_demonstrate-cmip6-eq-issue.py
# into this repo to demonstrate the issue with the equivalent species
```

```{code-cell}
---
editable: true
slideshow:
  slide_type: ''
tags: [remove_input]
---
fig, axes_d = get_default_delta_mosaic()
axes_d = remove_empty_axes(axes_d)

min_year = 1750
plot_overview_and_deltas(
    sel_times(ds_gases_full_d, lambda x: x.dt.year >= min_year),
    axes_d,
)
for ax in axes_d.values():
    xticks = [
        cftime.DatetimeProlepticGregorian(y, 1, 1) for y in np.arange(1750, 2050, 50)
    ]
    ax.set_xticks(xticks)
    ax.set_xticklabels([v.year for v in xticks])

plt.tight_layout()
plt.show()
```

+++ {"editable": true, "slideshow": {"slide_type": ""}}

### Global, annual-mean concentrations: Year 1957 - 2022

1957 is the start of the Scripps ground-based record.
Before this, data is based on ice cores alone.

```{code-cell}
---
editable: true
slideshow:
  slide_type: ''
tags: [remove_input]
---
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
```

### Global, monthly-mean concentrations: Year 1 - 2022

```{code-cell}
ds_gases_full_monthly_d = {}
for gas in gases_to_show:
    ds_gases_full_monthly_d[gas] = {}
    for source_id, cmip_era in (
        ("CR-CMIP-1-0-0", "CMIP7"),
        ("UoM-CMIP-1-2-0", "CMIP6"),
        # (None, "CMIP5"),
    ):
        query_kwargs = {
            "ghg": gas,
            "time_sampling": "mon",
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
        ds_gases_full_monthly_d[gas][cmip_era] = ds.compute()
```

```{code-cell}
fig, axes_d = get_default_delta_mosaic()
axes_d = remove_empty_axes(axes_d)

# min_year = 1957
min_year = 1990
min_year = 1
# min_year = 1750
plot_overview_and_deltas(
    sel_times(ds_gases_full_monthly_d, lambda x: x.dt.year >= min_year),
    axes_d,
)
for ax in axes_d.values():
    xticks = [
        cftime.DatetimeProlepticGregorian(y, 1, 1)
        # for y in np.arange(1750, 2050, 50)
        # for y in np.arange(1750, 1760, 1)
        for y in np.arange(1, 2050, 500)
        # for y in np.arange(1, 20, 1)
    ]
    ax.set_xticks(xticks)
    # ax.set_xlim(xticks[0], xticks[-1])
    ax.set_xticklabels([v.year for v in xticks])

plt.tight_layout()
plt.show()
```

### Latitudinally-resolved, monthly-mean concentrations: Year 1 - 2022

```{code-cell}
gases_to_show = ["co2", "ch4"]
ds_gases_full_monthly_lat_d = {}
for gas in gases_to_show:
    ds_gases_full_monthly_lat_d[gas] = {}
    for source_id, cmip_era in (
        ("CR-CMIP-1-0-0", "CMIP7"),
        ("UoM-CMIP-1-2-0", "CMIP6"),
        # (None, "CMIP5"),
    ):
        query_kwargs = {
            "ghg": gas,
            "time_sampling": "mon",
            "grid": "gnz",
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
        ds_gases_full_monthly_lat_d[gas][cmip_era] = ds.compute()
```

```{code-cell}
def plot_lat_selection(
    gas: str,
    ds_d: dict[str, dict[str, xr.Dataset]],
    ax: matplotlib.axes.Axes,
    ax_delta: matplotlib.axes.Axes,
    ax_delta_re: matplotlib.axes.Axes,
) -> None:
    """
    Plot selection for a latitude-specific dataset
    """
    target_unit_conc = ds_d[gas]["CMIP7"][gas].attrs["units"]
    target_unit_re = "W / m^2"

    for cmip_era, ds in ds_d[gas].items():
        label = f"{cmip_era} ({ds.attrs['source_id']})"
        tmp = ds[gas].copy()
        tmp.values = Q(tmp.values, tmp.attrs["units"]).to(target_unit_conc).m
        ds[gas].plot.scatter(ax=ax, label=label, alpha=0.7, edgecolors="none")

    ax.legend()
    ax.set_title(
        f"lat: {float(lat)}",
        # fontsize="small",
    )
    ax.xaxis.set_tick_params(labelbottom=False)
    ax.set_ylabel(target_unit_conc)
    ax.set_xlabel(None)

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
        ax_delta.set_title(None)

        tmp = RADIATIVE_EFFICIENCIES[gas] * Q(delta.values, delta.attrs["units"])
        delta_re = delta.copy()
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
        ax_delta_re.set_title(None)
```

```{code-cell}
gas = "co2"
# gas = "ch4"
min_year = 1
# min_year = 1750
# min_year = 1850
# min_year = 2000
sel_times_func = lambda x: (x.dt.year >= min_year)  # noqa: E731
# sel_times_func = lambda x: (x.dt.year >= min_year) & (x.dt.year <= min_year + 2)

ncols = 4
fig, axes = plt.subplots(ncols=ncols, nrows=9, figsize=(14, 16), sharex=True)
ax_flat = axes.flatten()

for i, lat in tqdm.auto.tqdm(
    enumerate(ds_gases_full_monthly_lat_d[gas]["CMIP7"]["lat"][::-1]), leave=False
):
    ax_idx = i % ncols + 3 * ncols * (i // ncols)
    # print(ax_idx)
    ax = ax_flat[ax_idx]
    ax_delta = ax_flat[ax_idx + ncols]
    ax_delta_re = ax_flat[ax_idx + 2 * ncols]

    plot_lat_selection(
        gas=gas,
        ds_d=sel_times(
            sel_lat(ds_gases_full_monthly_lat_d, lambda x: x == lat),
            sel_times_func,
        ),
        ax=ax,
        ax_delta=ax_delta,
        ax_delta_re=ax_delta_re,
    )
    # ax_flat[ax_idx].legend().remove()
    if gas == "co2":
        ax.set_ylim([250, 420])
        ax_delta.set_ylim([-2.5, 4.5])
        ax_delta_re.set_ylim([-0.03, 0.071])
    elif gas == "ch4":
        ax.set_ylim([600, 1900])
        ax_delta.set_ylim([-70, 35])
        ax_delta_re.set_ylim([-0.028, 0.02])
    # # break

plt.tight_layout()
# plt.savefig(f"{gas}_lat-monthly.png")
plt.suptitle(gas, y=1.0)
plt.show()
```
