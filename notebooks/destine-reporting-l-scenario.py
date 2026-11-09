# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.17.3
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %%
from functools import partial
from typing import Any

import cftime
import matplotlib.pyplot as plt
import nc_time_axis  # noqa: F401
import pandas as pd
import pandas_indexing as pix
import pandas_openscm
import seaborn as sns
import tqdm.auto
import xarray as xr

from local.data_loading import (
    fetch_and_load_ghg_dataset_scenarios,
)
from local.esgf.db_helpers import create_all_tables, get_sqlite_engine
from local.esgf.search.search_query import KnownIndexNode
from local.paths import REPO_ROOT

# %%
pandas_openscm.register_pandas_accessors()

# %%
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

# %%
fetch_and_load = partial(
    fetch_and_load_ghg_dataset_scenarios,
    local_data_root_dir=local_data_root_dir,
    # index_node=KnownIndexNode.DKRZ,
    # cmip_era="CMIP6",
    # source_id="UoM-CMIP-1-2-0",
    index_node=KnownIndexNode.ORNL,
)

# %%
gases_to_show = ["co2", "ch4", "n2o", "cfc12eq", "hfc134aeq"]
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
        ds = fetch_and_load(**query_kwargs)

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


# %%
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


# %%
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

# %%
start_year = 2000

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

# %%
from matplotlib.lines import Line2D

# %%
col = "ghg"

data = pdf.openscm.to_long_data()
data = data[data["experiment"].isin(["historical", "l"])]

fg = sns.relplot(
    data=data,
    x="time",
    y="value",
    hue="experiment",
    palette={k: v for k, v in palette.items() if k in data["experiment"].unique()},
    hue_order=[v for v in hue_order if v in data["experiment"].unique()],
    kind="scatter",
    col=col,
    col_wrap=3,
    facet_kws=dict(sharey=False),
    height=2.5,
    aspect=1.0,
    legend=False,
)

plt.tight_layout()
for i, ax in enumerate(fg.axes.flatten()):
    col_val = ax.get_title().split(f"{col} = ")[1]
    unit_l = (
        pdf.loc[pix.ismatch(**{col: col_val})].index.get_level_values("unit").unique()
    )
    if len(unit_l) > 1:
        raise AssertionError(unit_l)
    unit = unit_l[0]
    ax.set_ylabel(unit)
    if i == len(fg.axes.flatten()) - 1:
        custom_handles = [
            Line2D(
                [0],
                [0],
                color=palette["historical"],
                marker="o",
                linestyle="None",
                label="historical",
            ),
            Line2D(
                [0], [0], color=palette["l"], marker="o", linestyle="None", label="l"
            ),
            # Line2D([0], [0], color='red', linestyle='--', label='Dashed Line'),
            # Line2D([0], [0], color='green', marker='o', markersize=10, linestyle='None', label='Marker Only')
        ]
        ax.legend(
            handles=custom_handles,
            bbox_to_anchor=(1.05, 0.5),
            loc="center left",
            title="Experiment",
        )

plt.show()

# %%
gases_to_show = ["co2", "ch4", "n2o", "cfc12eq", "hfc134aeq"]
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
        ds = ds.compute()

        # avoid extensions
        ds = ds.sel(
            scenario=[v for v in ds["scenario"].values if "ext" not in v]
        ).dropna("time")

        ds_gases_full_monthly_d[gas][key] = ds

# %%
for gas, ds_d in ds_gases_full_monthly_d.items():
    pdf_l = [
        to_data_frame(
            ds_d["history"].sel(scenario="CMIP"),
            unstack_col="time",
            assign_metadata={"experiment": "history"},
        ),
        to_data_frame(
            ds_d["scenarios"].sel(scenario="l"),
            unstack_col="time",
            assign_metadata={"experiment": "l"},
        ),
    ]

    pdf = pd.concat([v.reorder_levels(pdf_l[0].index.names) for v in pdf_l]).sort_index(
        axis="columns"
    )
    pdf.columns = [v.year + (v.month * 2 - 1) / 24 for v in pdf.columns]
    # pdf

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
        palette=palette,
        height=2.0,
        aspect=1.0,
    )

    unit_l = pdf.index.get_level_values("unit").unique()
    if len(unit_l) > 1:
        raise AssertionError(unit_l)
    unit = unit_l[0]
    for ax in fg.axes.flatten():
        if ax.get_ylabel():
            ax.set_ylabel(unit)

    # glue("ds-co2-transition-from-history-monthly-lat-fig", fg.fig, display=False)
    plt.suptitle(gas, y=1.03)
    plt.show()

# %%
