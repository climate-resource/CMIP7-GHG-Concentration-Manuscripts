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
# # Scale latitudinal gradient EOFs

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ## Imports

# %% editable=true slideshow={"slide_type": ""}
from functools import partial
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import openscm_units
import pandas as pd
import pandas_indexing as pix  # noqa: F401
import pandas_openscm
import pint
import pint_xarray
import seaborn as sns
import xarray as xr

from cmip7_scenariomip_ghg_generation.mean_preserving_interpolation import (
    LaiKaplanInterpolator,
    interpolate_annual_mean_to_monthly,
)
from cmip7_scenariomip_ghg_generation.mean_preserving_interpolation.lai_kaplan import (
    get_wall_control_points_y_cubic,
)

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ## Parameters

# %% editable=true slideshow={"slide_type": ""} tags=["parameters"]
ghg: str = "ch4"
annual_mean_emissions_file: str = (
    "../output-bundles/dev-test/data/interim/single-variable-files/ch4_eof-one-scaling.feather"
)
harmonisation_year: int = 2015
historical_data_root_dir: str = "../output-bundles/dev-test/data/raw/historical-ghg-concs"
historical_data_seasonality_lat_gradient_info_root: str = (
    "../output-bundles/dev-test/data/raw/historical-ghg-data-interim"
)
out_file: str = "../output-bundles/dev-test/data/interim/latitudinal-gradient/ch4_latitudinal-gradient-info.nc"


# %% tags=["injected-parameters"]
# Parameters
ghg = "co2"
annual_mean_emissions_file = "manuscript-outputs/rerun-markers/single-variable-files/co2_eof-one-scaling.feather"
harmonisation_year = 2023
historical_data_root_dir = "../../cmip-ghg-concentration-generation/v1.0.0/data/processed/esgf-ready/input4MIPs/CMIP7/CMIP/CR/CR-CMIP-1-0-0"
historical_data_seasonality_lat_gradient_info_root = "../../cmip-ghg-concentration-generation/v1.0.0/data/interim"
out_file = "manuscript-outputs/rerun-markers/latitudinal-gradient/co2_latitudinal-gradient-info.nc"

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ## Parse parameters

# %% editable=true slideshow={"slide_type": ""}
annual_mean_emissions_file_p = Path(annual_mean_emissions_file)
historical_data_root_dir_p = Path(historical_data_root_dir)
historical_data_seasonality_lat_gradient_info_root_p = Path(historical_data_seasonality_lat_gradient_info_root)
out_file_p = Path(out_file)

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ## Set up

# %% editable=true slideshow={"slide_type": ""}
ur = pint_xarray.setup_registry(openscm_units.unit_registry)
Q = ur.Quantity
pandas_openscm.register_pandas_accessor()

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ## Load data

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ### Annual-mean emissions

# %% editable=true slideshow={"slide_type": ""}
annual_mean_emissions = pd.read_feather(annual_mean_emissions_file_p)
ghg_search_for = ghg

annual_mean_emissions_emms_units = annual_mean_emissions.loc[
    annual_mean_emissions.index.get_level_values("unit").str.lower().str.contains(ghg_search_for)
]
if len(annual_mean_emissions_emms_units.pix.unique("unit")) != 1:
    raise AssertionError

annual_mean_emissions_emms_u = annual_mean_emissions_emms_units.pix.unique("unit")[0]

annual_mean_emissions_emms_units


# %%
# Take any scenario, doesn't matter as all the same pre-harmonisation
annual_mean_emissions_emms_units_historical = annual_mean_emissions_emms_units.loc[:, :harmonisation_year].iloc[:1, :]
annual_mean_emissions_emms_units_historical

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ### CMIP7 historical GHG concentrations


# %% editable=true slideshow={"slide_type": ""}
def load_file_from_glob(glob: str, base_dir: Path) -> xr.Dataset:
    """
    Load a single file based on a glob pattern
    """
    file_l = list(base_dir.rglob(glob))
    if len(file_l) != 1:
        raise AssertionError(file_l)

    ds = xr.load_dataset(file_l[0])

    return ds


# %% editable=true slideshow={"slide_type": ""}
cmip7_historical_gm_monthly_ds = load_file_from_glob(f"*{ghg}_*gm_175001-*.nc", historical_data_root_dir_p)
cmip7_historical_gm_monthly = cmip7_historical_gm_monthly_ds[ghg]
# cmip7_historical_gm_monthly

# %% [markdown]
# ### CMIP7 historical latitudinal gradient

# %%
cmip7_lat_gradient_pieces_ds = load_file_from_glob(
    f"{ghg}_allyears-lat-gradient-eofs-pcs.nc", historical_data_seasonality_lat_gradient_info_root_p
).pint.quantify(unit_registry=ur)
# cmip7_lat_gradient_pieces_ds

# %% [markdown]
# ## Regress emissions against PC0

# %%
pc0 = cmip7_lat_gradient_pieces_ds["principal-components"].sel(eof=0)
start_non_constant_pc0 = np.where(pc0.diff("year") != 0.0)[0][0]
start_non_constant_pc0

# %%
annual_mean_emissions_emms_non_zero_years = (
    annual_mean_emissions_emms_units_historical[annual_mean_emissions_emms_units_historical > 0.0]
    .dropna(axis="columns")
    .columns
)
annual_mean_emissions_emms_non_zero_years

# %%
regression_years = np.intersect1d(
    pc0.sel(year=pc0["year"] >= start_non_constant_pc0)["year"].values,
    annual_mean_emissions_emms_non_zero_years,
)
regression_years

# %%
pc0_to_regress = pc0.sel(year=regression_years)
# pc0_to_regress

# %%
annual_mean_emissions_emms_units_historical_to_regress_q = Q(
    annual_mean_emissions_emms_units_historical.loc[:, regression_years].values.squeeze(), annual_mean_emissions_emms_u
)
# annual_mean_emissions_emms_units_historical_to_regress_q

# %%
x = annual_mean_emissions_emms_units_historical_to_regress_q
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
ax.set_xlabel("emissions")
ax.legend()

# %% [markdown]
# ## Scale latitudinal gradient PCs

# %% [markdown]
# ### PC0

# %%
# So ugly.
# Basically, we need to keep things consistent
# at the join point.
# For ScenarioMIP, that's 2022: we want everything the same in 2022.
# For PolMIP, that's 2015: we want everything the same in 2015,
# divergence thereafter.
last_hist_year = min(
    cmip7_historical_gm_monthly["time"].dt.year.values[-1],
    harmonisation_year,
)
last_hist_year

# %%
delta_E = (
    annual_mean_emissions_emms_units.loc[:, last_hist_year:].subtract(
        annual_mean_emissions_emms_units[last_hist_year], axis="rows"
    )
).pix.assign(variable="change_in_emissions")
delta_E

# %%
delta_E_unit_l = delta_E.pix.unique("unit")
if len(delta_E_unit_l) != 1:
    raise AssertionError(delta_E_unit_l)

delta_E_unit = delta_E_unit_l[0]

delta_E_xr = xr.DataArray(
    delta_E.values,
    dims=["scenario", "year"],
    coords=dict(scenario=delta_E.index.get_level_values("scenario"), year=delta_E.columns),
).pint.quantify(delta_E_unit, unit_registry=ur)
# delta_E_xr

# %%
delta_pc = delta_E_xr * m
# delta_pc

# %%
pc0_extended = delta_pc + pc0.sel(year=last_hist_year)

fig, axes = plt.subplots(nrows=2, figsize=(8, 8))
annual_mean_emissions_emms_units.T.plot(ax=axes[0])
pc0_extended.pint.to("dimensionless").plot(ax=axes[1], hue="scenario")
pc0.pint.to("dimensionless").plot(ax=axes[1])
for ax in axes:
    sns.move_legend(ax, loc="center left", bbox_to_anchor=(1.05, 0.5))

# pc0_extended

# %% [markdown]
# ### PC1

# %%
pc1 = cmip7_lat_gradient_pieces_ds["principal-components"].sel(eof=1)
pc1_extended = pc1.pint.dequantify().interp(
    year=annual_mean_emissions.columns, kwargs={"fill_value": pc1.sel(year=last_hist_year).data.m}
)

pc1_extended_all_scenarios_l = []
for scenario in pc0_extended["scenario"]:
    pc1_extended_all_scenarios_l.append(pc1_extended.assign_coords(scenario=scenario.data))

pc1_extended_all_scenarios = xr.concat(pc1_extended_all_scenarios_l, dim="scenario")
pc1_extended_all_scenarios

pc1_extended_all_scenarios.plot(hue="scenario")

pc1_extended_all_scenarios = pc1_extended_all_scenarios.sel(
    year=annual_mean_emissions.loc[:, last_hist_year:].columns
).pint.quantify(unit_registry=ur)

# pc1_extended_all_scenarios

# %%
exp_n_eofs = 2
if cmip7_lat_gradient_pieces_ds["eof"].size != exp_n_eofs:
    raise AssertionError


# %% [markdown]
# ## Interpolate to monthly

# %% [markdown]
# We have to ensure that the scenarios stay harmonised until the last month of the historical data.
# Therefore, we use the same wall control point value for all scenarios.
# As a simple choice, we use the average wall control point across all scenarios.


# %%
def get_wall_control_points(
    intervals_x: pint.UnitRegistry.Quantity,
    intervals_y: pint.UnitRegistry.Quantity,
    control_points_wall_x: pint.UnitRegistry.Quantity,
    fixed_control_point: pint.UnitRegistry.Quantity,
) -> pint.UnitRegistry.Quantity:
    """
    Get wall control points including setting our fixed control point
    """
    # Start off with standard implementation
    control_points_wall_y = get_wall_control_points_y_cubic(
        intervals_x=intervals_x,
        intervals_y=intervals_y,
        control_points_wall_x=control_points_wall_x,
    )

    # Also fix the control point between overlap year month 12 and overlap year + 1 month 1
    fixed_control_point_idxr = np.where(
        control_points_wall_x == openscm_units.unit_registry.Quantity(last_hist_year + 1, "yr")
    )
    control_points_wall_y[fixed_control_point_idxr] = fixed_control_point

    return control_points_wall_y


# %%
pcs_extended_monthly_l = []
for pc_extended in [
    pc0_extended,
    pc1_extended_all_scenarios,
]:
    tmp_l = []
    fixed_control_point = (
        pc_extended.sel(year=[last_hist_year, last_hist_year + 1]).mean("year").mean("scenario").data.squeeze()
    )

    for scenario, sda in pc_extended.groupby("scenario"):
        scenario_monthly = interpolate_annual_mean_to_monthly(
            values=sda.data.m.squeeze(),
            values_units=pc_extended.data.u,
            years=pc_extended["year"].values,
            algorithm=LaiKaplanInterpolator(
                progress_bar=True,
                get_wall_control_points_y_from_interval_ys=partial(
                    get_wall_control_points, fixed_control_point=fixed_control_point
                ),
            ),
            unit_registry=openscm_units.unit_registry,
        ).assign_coords(scenario=scenario, eof=pc_extended["eof"].values)
        tmp_l.append(scenario_monthly)

    pcs_extended_monthly_l.append(xr.concat(tmp_l, dim="scenario"))

pcs_extended_monthly = xr.concat(pcs_extended_monthly_l, dim="eof")
pcs_extended_monthly.name = "principal-components-monthly"
pcs_extended_monthly.attrs["description"] = "principal component values on a monthly timestep"
pcs_extended_monthly

# %%
months_per_year = 12

for _, pc_extended_monthly in pcs_extended_monthly.groupby("eof"):
    last_overlap_month = pc_extended_monthly.sel(
        time=(pc_extended_monthly["time"].dt.year == last_hist_year)
        & (pc_extended_monthly["time"].dt.month == months_per_year)
    )

    np.testing.assert_allclose(
        last_overlap_month.data.m.squeeze(),
        last_overlap_month.data.m[0].squeeze(),
        rtol=1e-3,
    )

# %%
for years in (
    np.arange(last_hist_year + 1, 2500 + 1),
    np.arange(last_hist_year + 1, last_hist_year + 15),
):
    pcs_extended_monthly.sel(time=pcs_extended_monthly["time"].dt.year.isin(years)).pint.to("dimensionless").plot(
        alpha=0.6,
        hue="scenario",
        col="eof",
    )

    plt.show()

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ## Prepare output

# %%
eofs = cmip7_lat_gradient_pieces_ds["eofs"]

out = xr.merge(
    [
        eofs,
        pcs_extended_monthly.pint.to("dimensionless"),
    ],
    combine_attrs="drop_conflicts",
).pint.dequantify()

out

# %%
# Double check this can work
out_tmp = out.pint.quantify(unit_registry=ur)
tmp = (out_tmp["eofs"] * out_tmp["principal-components-monthly"]).sum("eof")
# tmp

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

manuscript_regression_file = Path("manuscript-outputs/rerun-markers/regressions/co2_regression.csv")
manuscript_regression_file.parent.mkdir(exist_ok=True, parents=True)
pd.DataFrame({"year": regression_years, "x": x.m, "y": y.m}).to_csv(
    manuscript_regression_file, index=False
)
Path("manuscript-outputs/rerun-markers/regressions/co2_regression.json").write_text(
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
    (1, pc1_extended_all_scenarios),
):
    manuscript_pc = manuscript_pc.pint.to("dimensionless").pint.dequantify()
    manuscript_pc = manuscript_pc.sel(year=manuscript_pc["year"] >= last_hist_year)
    manuscript_pc_df = (
        manuscript_pc.rename("value").to_dataframe()[["value"]].reset_index()
    )
    manuscript_pc_df["eof"] = manuscript_eof
    manuscript_pcs_l.append(manuscript_pc_df[["eof", "scenario", "year", "value"]])

manuscript_pcs_file = Path("manuscript-outputs/rerun-markers/latitudinal-gradient/co2_pcs-annual.csv")
manuscript_pcs_file.parent.mkdir(exist_ok=True, parents=True)
pd.concat(manuscript_pcs_l).to_csv(manuscript_pcs_file, index=False)
manuscript_pcs_file
