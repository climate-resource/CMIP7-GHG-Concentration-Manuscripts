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
import numpy as np
import pandas as pd
import pandas_indexing as pix
import pandas_openscm
import seaborn as sns
import tqdm.auto
import xarray as xr

from local.data_loading import (
    fetch_and_load_ghg_dataset_scenarios,
    get_ghg_dataset_local_files,
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
query_kwargs_co2_yearly_global = dict(
    ghg="co2",
    time_sampling="yr",
    grid="gm",
    cmip_era="CMIP7",
    source_id="CR-CMIP-1-0-0",
    engine=engine,
)

time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)

fetch_and_load = partial(
    fetch_and_load_ghg_dataset_scenarios,
    local_data_root_dir=local_data_root_dir,
    # index_node=KnownIndexNode.DKRZ,
    # cmip_era="CMIP6",
    # source_id="UoM-CMIP-1-2-0",
    index_node=KnownIndexNode.ORNL,
)

# %%
gases_to_show = ["co2"]

ds_gases_full_yearly_multi_phase_d = {}
for gas in gases_to_show:
    ds_gases_full_yearly_multi_phase_d[gas] = {}
    for key, target_mip, source_version, institution_id, cmip_era in (
        ("history-cmip7", "CMIP", "1.0.0", "CR", "CMIP7"),
        ("scenarios-cmip7", "ScenarioMIP", "1.0.0", "CR", "CMIP7"),
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
    # "vllo": "low",
    # "vlho": "low",
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
    # "vlho",
    "vl",
    # "vllo",
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

# %%
start_year = 1850
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

# %%
palette = {
    "history": "k",
    "history-cmip6": "tab:grey",
    "historical": "k",
    "historical-cmip6": "tab:grey",
    "vl": "#24a4ff",
    # "vllo": "#24a4ff",
    "ln": "#4a0daf",
    # "vlho": "#4a0daf",
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
    # "l",
    "ln",
    # "vlho",
    "vl",
    "vllo",
    "historical",
]

# %%
start_year = 1850
end_year = 2100


pdfh = pdf.loc[
    pix.ismatch(cmip_era="CMIP7") | pix.ismatch(experiment="historical"),
    (pdf.columns >= start_year) & (pdf.columns <= end_year + 1),
]
fg = sns.relplot(
    data=pdfh.openscm.to_long_data(),
    x="time",
    y="value",
    hue="experiment",
    palette={k: v for k, v in palette.items() if k in pdfh.pix.unique("experiment")},
    hue_order=[v for v in hue_order_incl_cmip6 if v in pdfh.pix.unique("experiment")],
    style="cmip_era",
    markers={"CMIP6": "+", "CMIP7": "x"},
    # edgecolor="none",
    alpha=0.6,
    s=40,
    kind="scatter",
    # row=row,
    # row="scenario_group",
    row_order=["low", "continuing-trends", "high"],
    facet_kws=dict(sharey=False),
    height=3.0,
    aspect=2.0,
)

for h in fg._legend.legend_handles:
    if h.get_label() == "CMIP6":
        h.set_marker("+")
        h.set_alpha(1.0)
        h.set_markeredgewidth(1.0)
    if h.get_label() == "CMIP7":
        h.set_marker("x")
        h.set_alpha(1.0)
        h.set_markeredgewidth(1.0)

for ax in fg.axes.flatten():
    ax.set_xlabel("")
    ax.set_ylabel("co2 [ppm]")
    fg.figure.suptitle("CMIP7 GHG forcing (available scenarios)")
    ax.set_title("CMIP6 historical is shown too", size="small", y=0.85)

# plt.tight_layout()

# %%
query_kwargs_co2_monthly_lat = {
    **query_kwargs_co2_yearly_global,
    "time_sampling": "mon",
    "grid": "gnz",
}

# Ensure data is downloaded
_ = fetch_and_load(**query_kwargs_co2_monthly_lat)

# Get file paths
co2_monthly_lat_fps = get_ghg_dataset_local_files(**query_kwargs_co2_monthly_lat)

ds_co2_monthly_lat = xr.open_mfdataset(
    co2_monthly_lat_fps, decode_times=time_coder, data_vars=None, compat="no_conflicts"
)
# Force values to compute to avoid dask getting involved
ds_co2_monthly_lat = ds_co2_monthly_lat.compute()

# %%
fig = plt.figure(figsize=(12, 6))
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

plt.title("Illustration of seasonality and latitudinal gradient")
plt.tight_layout()
plt.show()

# %% [markdown]
# Unused plots

# %%
start_year = 2000
end_year = 2100

fg = sns.relplot(
    data=pdf.loc[
        :, (pdf.columns >= start_year) & (pdf.columns <= end_year + 1)
    ].openscm.to_long_data(),
    x="time",
    y="value",
    hue="experiment",
    palette={k: v for k, v in palette.items() if k in pdf.pix.unique("experiment")},
    hue_order=[v for v in hue_order_incl_cmip6 if v in pdf.pix.unique("experiment")],
    style="cmip_era",
    markers={"CMIP6": "+", "CMIP7": 5},
    # edgecolor="none",
    alpha=0.6,
    s=50,
    kind="scatter",
    # row=row,
    row="scenario_group",
    row_order=["low", "continuing-trends", "high"],
    facet_kws=dict(sharey=False),
    height=2.0,
    aspect=2.0,
)

for h in fg._legend.legend_handles:
    if h.get_label() == "CMIP6":
        h.set_marker("+")
        h.set_alpha(1.0)
        h.set_markeredgewidth(1.0)

for ax in fg.axes.flatten():
    if ax.get_ylabel():
        gas = "co2"
        unit_l = (
            pdf.loc[(pdf.index.get_level_values("ghg") == gas)]
            .index.get_level_values("unit")
            .unique()
        )
        if len(unit_l) > 1:
            raise AssertionError(unit_l)
        unit = unit_l[0]
        ax.set_ylabel(unit)

# plt.tight_layout()
plt.show()
