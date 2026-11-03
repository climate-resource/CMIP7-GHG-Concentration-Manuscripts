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

# %% [markdown] editable=true slideshow={"slide_type": ""}
# # Scale seasonality based on MAGICC NPP

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ## Imports

# %% editable=true slideshow={"slide_type": ""}
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import openscm_units
import pandas_indexing as pix
import pandas_openscm
import pandas_openscm.db
import pint_xarray
import seaborn as sns
import tqdm.auto
import xarray as xr

from cmip7_scenariomip_ghg_generation.scenario_info import ScenarioInfo

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ## Parameters

# %% editable=true slideshow={"slide_type": ""} tags=["parameters"]
ghg: str = "co2"
scenario_info_markers: str = (
    "WITCH 6.0;SSP5 - Medium-Low Emissions_a;hl;;"
    "REMIND-MAgPIE 3.5-4.10;SSP1 - Very Low Emissions;vl;;"
    "MESSAGEix-GLOBIOM-GAINS 2.1-M-R12;SSP2 - Low Emissions;l;;"
    "IMAGE 3.4;SSP2 - Medium Emissions;m;;"
    "GCAM 7.1 scenarioMIP;SSP3 - High Emissions;h;;"
    "AIM 3.0;SSP2 - Low Overshoot;ln;;"
    "COFFEE 1.6;SSP2 - Medium-Low Emissions;ml"
)
harmonisation_year: int = 2023
magicc_output_db_dir: str = "../output-bundles/dev-test/data/interim/magicc-output/db"
magicc_db_backend_str: str = "feather"
historical_data_seasonality_lat_gradient_info_root: str = (
    "../output-bundles/dev-test/data/raw/historical-ghg-data-interim"
)
out_file: str = (
    "../output-bundles/dev-test/data/interim/seasonality/modelling-based-projection_co2_seasonality-all-time.nc"
)


# %% tags=["injected-parameters"]
# Parameters
ghg = "co2"
scenario_info_markers = "REMIND-MAgPIE 3.5-4.11;VL-cf;vl-cf"
harmonisation_year = 2015
magicc_output_db_dir = "data/interim/magicc-output/db"
magicc_db_backend_str = "feather"
historical_data_seasonality_lat_gradient_info_root = "../../cmip-ghg-concentration-generation/v1.0.0/data/interim"
out_file = "manuscript-outputs/rerun-markers/seasonality/modelling-based-projection_co2_seasonality-all-time_vl-cf.nc"

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ## Parse parameters

# %%
siml = []
for v in scenario_info_markers.split(";;"):
    model, scenario, cmip_scenario_name = v.split(";")
    siml.append(
        ScenarioInfo(
            cmip_scenario_name=cmip_scenario_name,
            model=model,
            scenario=scenario,
        )
    )

scenario_info_markers_p = tuple(siml)
historical_data_seasonality_lat_gradient_info_root_p = Path(historical_data_seasonality_lat_gradient_info_root)
magicc_output_db = pandas_openscm.db.OpenSCMDB(
    backend_data=pandas_openscm.db.DATA_BACKENDS.get_instance(magicc_db_backend_str),
    backend_index=pandas_openscm.db.INDEX_BACKENDS.get_instance(magicc_db_backend_str),
    db_dir=Path(magicc_output_db_dir),
)
out_file_p = Path(out_file)

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ## Set up

# %% editable=true slideshow={"slide_type": ""}
ur = pint_xarray.setup_registry(openscm_units.unit_registry)
Q = ur.Quantity
pandas_openscm.register_pandas_accessor()

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ## Load data

# %% [markdown]
# ### MAGICC NPP

# %%
magicc_output_db_reader = magicc_output_db.create_reader()
# magicc_output_db_reader

# %%
magiccc_output_l = []
for si in tqdm.auto.tqdm(scenario_info_markers_p):
    magiccc_output_l.append(
        magicc_output_db_reader.load(
            pix.isin(model=si.model, scenario=si.scenario, climate_model="MAGICCv7.6.0a3", variable="CO2_CURRENT_NPP")
            # progress=True,
        ).pix.assign(scenario=si.cmip_scenario_name)
    )

magiccc_output = pix.concat(magiccc_output_l)
magiccc_output_median = magiccc_output.openscm.groupby_except("run_id").median()
# magiccc_output_median.pix.project(["model", "scenario"]).T.plot()
# magiccc_output_median

# %%
# magiccc_output_median.pix.project(["model", "scenario"]).T.plot()

# Take any scenario, doesn't matter as all the same pre-harmonisation
magiccc_output_median_historical = magiccc_output_median.loc[:, :harmonisation_year].iloc[:1, :]
magiccc_output_median_historical


# %% [markdown]
# ### CMIP7 historical seasonality


# %%
def load_file_from_glob(glob: str, base_dir: Path) -> xr.Dataset:
    """
    Load a single file based on a glob pattern
    """
    file_l = list(base_dir.rglob(glob))
    if len(file_l) != 1:
        raise AssertionError(file_l)

    ds = xr.load_dataset(file_l[0])

    return ds


# %%
cmip7_seasonality_base = (
    load_file_from_glob(
        f"{ghg}_observational-network_seasonality.nc", historical_data_seasonality_lat_gradient_info_root_p
    )
    .pint.quantify(unit_registry=ur)
    .rename({"__xarray_dataarray_variable__": "seasonality"})
)
cmip7_seasonality_base

# %%
cmip7_seasonality_pieces_ds = load_file_from_glob(
    f"{ghg}_observational-network_seasonality-change-eofs.nc", historical_data_seasonality_lat_gradient_info_root_p
).pint.quantify(unit_registry=ur)
cmip7_seasonality_pieces_ds["principal-components"].plot(hue="eof")
# cmip7_seasonality_pieces_ds

# %% [markdown]
# ## Scale seasonality pc

# %%
if cmip7_seasonality_pieces_ds["eof"].size != 1:
    raise AssertionError

pc0 = cmip7_seasonality_pieces_ds["principal-components"].sel(eof=0)
# pc0

# %%
regression_years = np.intersect1d(
    pc0["year"],
    magiccc_output_median_historical.columns,
)
regression_years

# %%
pc0_to_regress = pc0.sel(year=regression_years)
# pc0_to_regress

# %%
magiccc_output_median_historical_u_l = magiccc_output_median_historical.pix.unique("unit")
if len(magiccc_output_median_historical_u_l) != 1:
    raise AssertionError(magiccc_output_median_historical_u_l)

magiccc_output_median_historical_u = magiccc_output_median_historical_u_l[0]
# magiccc_output_median_historical_u

# %%
magiccc_output_median_historical_to_regress_q = Q(
    magiccc_output_median_historical.loc[:, regression_years].values.squeeze(),
    magiccc_output_median_historical_u,
)
# magiccc_output_median_historical_to_regress_q

# %%
x = magiccc_output_median_historical_to_regress_q
A = np.vstack([x.m, np.ones(x.size)]).T
y = pc0_to_regress.data

res = np.linalg.lstsq(A, y.m, rcond=None)
m, c = res[0]
m = Q(m, (y / x).units)
# c = Q(c, y.units)

fig, ax = plt.subplots()
ax.scatter(x.m, y.m, label="raw data")
ax.plot(x.m, (m * x + c).m, color="tab:orange", label="regression")
ax.set_ylabel("PC0")
ax.set_xlabel("NPP")
ax.legend()

# %% [markdown]
# ## Scale PC

# %%
last_hist_year = int(cmip7_seasonality_pieces_ds["year"].max().values)
# last_hist_year

# %%
delta_npp = (
    magiccc_output_median.loc[:, last_hist_year:].subtract(magiccc_output_median[last_hist_year], axis="rows")
).pix.assign(variable="change_in_npp")
# delta_npp

# %%
delta_npp_unit_l = delta_npp.pix.unique("unit")
if len(delta_npp_unit_l) != 1:
    raise AssertionError(delta_npp_unit_l)

delta_npp_unit = delta_npp_unit_l[0]

delta_npp_xr = xr.DataArray(
    delta_npp.values,
    dims=["scenario", "year"],
    coords=dict(scenario=delta_npp.index.get_level_values("scenario"), year=delta_npp.columns.values),
).pint.quantify(delta_npp_unit, unit_registry=ur)
# delta_npp_xr

# %%
delta_pc = delta_npp_xr * m
# delta_pc

# %%
pc0_extended = delta_pc + pc0.sel(year=last_hist_year).data
pc0_extended = xr.concat(
    [
        pc0.sel(year=pc0["year"] < pc0_extended["year"].min()),
        pc0_extended,
    ],
    dim="year",
)

fig, axes = plt.subplots(nrows=2, figsize=(8, 8))
magiccc_output_median.pix.project("scenario").T.plot(ax=axes[0])
pc0_extended.pint.to("dimensionless").plot(ax=axes[1], hue="scenario")
pc0.pint.to("dimensionless").plot(ax=axes[1], linestyle="--", alpha=0.5, marker="x")
for ax in axes:
    sns.move_legend(ax, loc="center left", bbox_to_anchor=(1.05, 0.5))

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ## Prepare output

# %%
out_seasonality_delta = (cmip7_seasonality_pieces_ds["eofs"] * pc0_extended).sum("eof")
out = out_seasonality_delta + cmip7_seasonality_base["seasonality"]
out = out.pint.dequantify()
# out

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ## Save

# %% editable=true slideshow={"slide_type": ""}
out_file_p.parent.mkdir(exist_ok=True, parents=True)
out.to_netcdf(out_file_p)
out_file_p

# %%
# Added for the CMIP7 ScenarioMIP GHG manuscript.
# The original run only plots the regression, we need its data.
import json

import pandas as pd

manuscript_regression_file = Path("manuscript-outputs/rerun-markers/regressions/co2-seasonality-change-vl-cf_regression.csv")
manuscript_regression_file.parent.mkdir(exist_ok=True, parents=True)
pd.DataFrame({"year": regression_years, "x": x.m, "y": y.m}).to_csv(
    manuscript_regression_file, index=False
)
Path("manuscript-outputs/rerun-markers/regressions/co2-seasonality-change-vl-cf_regression.json").write_text(
    json.dumps(
        {
            "x_unit": str(x.units),
            "y_unit": str(y.units),
            "slope": float(m.m),
            "slope_unit": str(m.units),
            "intercept": float(c),
        },
        indent=2,
    )
)
manuscript_regression_file

# %%
# Added for the CMIP7 ScenarioMIP GHG manuscript.
# The original run only saves the monthly PCs, we want the annual ones too.
import pandas as pd

manuscript_pcs_l = []
for manuscript_eof, manuscript_pc in (
    (0, pc0_extended),
):
    manuscript_pc = manuscript_pc.pint.to("dimensionless").pint.dequantify()
    manuscript_pc = manuscript_pc.sel(year=manuscript_pc["year"] >= last_hist_year)
    manuscript_pc_df = (
        manuscript_pc.rename("value").to_dataframe()[["value"]].reset_index()
    )
    manuscript_pc_df["eof"] = manuscript_eof
    manuscript_pcs_l.append(manuscript_pc_df[["eof", "scenario", "year", "value"]])

manuscript_pcs_file = Path("manuscript-outputs/rerun-markers/seasonality/co2_seasonality-change-pc_vl-cf.csv")
manuscript_pcs_file.parent.mkdir(exist_ok=True, parents=True)
pd.concat(manuscript_pcs_l).to_csv(manuscript_pcs_file, index=False)
manuscript_pcs_file
