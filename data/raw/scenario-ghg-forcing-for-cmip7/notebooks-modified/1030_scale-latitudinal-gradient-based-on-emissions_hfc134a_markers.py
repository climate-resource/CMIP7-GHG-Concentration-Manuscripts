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
# # Scale latitudinal gradient based on emissions

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ## Imports

# %% editable=true slideshow={"slide_type": ""}
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
import yaml

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
ghg: str = "c8f18"
annual_mean_emissions_file: str = "../output-bundles/dev-test/data/interim/single-variable-files/c8f18_total.feather"
historical_data_root_dir: str = "../output-bundles/dev-test/data/raw/historical-ghg-concs"
historical_data_seasonality_lat_gradient_info_root: str = (
    "../output-bundles/dev-test/data/raw/historical-ghg-data-interim"
)
out_file: str = "../output-bundles/dev-test/data/interim/latitudinal-gradient/c8f18_latitudinal-gradient-info.nc"


# %% tags=["injected-parameters"]
# Parameters
ghg = "hfc134a"
annual_mean_emissions_file = "manuscript-outputs/rerun-markers/single-variable-files/hfc134a_total.feather"
historical_data_root_dir = "../../cmip-ghg-concentration-generation/v1.0.0/data/processed/esgf-ready/input4MIPs/CMIP7/CMIP/CR/CR-CMIP-1-0-0"
historical_data_seasonality_lat_gradient_info_root = "../../cmip-ghg-concentration-generation/v1.0.0/data/interim"
out_file = "manuscript-outputs/rerun-markers/latitudinal-gradient/hfc134a_latitudinal-gradient-info.nc"

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
if "hfc4310" in ghg:
    ghg_search_for = "hfc4310"
else:
    ghg_search_for = ghg

annual_mean_emissions_emms_units = annual_mean_emissions.loc[
    annual_mean_emissions.index.get_level_values("unit").str.lower().str.contains(ghg_search_for)
]
if len(annual_mean_emissions_emms_units.pix.unique("unit")) != 1:
    raise AssertionError

annual_mean_emissions_emms_units


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
if ghg != "c8f18":
    cmip7_lat_gradient_pieces_ds = load_file_from_glob(
        f"{ghg}_allyears-lat-gradient-eofs-pcs.nc", historical_data_seasonality_lat_gradient_info_root_p
    ).pint.quantify(unit_registry=ur)

else:
    # Assume constant latitudinal gradient
    # because the historical regression doesn't make sense
    # and we didn't save the data out.
    pass

# %% [markdown]
# ### CMIP7 historical latitudinal gradient regression info

# %%
if ghg != "c8f18":
    file_l = list(
        historical_data_seasonality_lat_gradient_info_root_p.rglob(f"{ghg}_pc0-total-emissions-regression.yaml")
    )
    if len(file_l) != 1:
        raise AssertionError(file_l)
    with open(file_l[0]) as fh:
        reg_info = yaml.safe_load(fh)
else:
    reg_info = {"m": (0.0, "yr / (kt C8F18)")}
# reg_info

# %% [markdown]
# ## Scale latitudinal gradient pc

# %%
last_hist_year = cmip7_historical_gm_monthly["time"].dt.year.values[-1]
# last_hist_year

# %%
if ghg != "c8f18":
    pc = cmip7_lat_gradient_pieces_ds["principal-components"]
else:
    years = cmip7_historical_gm_monthly.groupby("time.year").mean()["year"]
    pc = xr.DataArray(
        Q(np.ones((years.size, 1)), "dimensionless"),
        dims=("year", "eof"),
        coords=dict(year=years, eof=[0]),
        name="principal-components",
    ).transpose("eof", "year")

if pc["eof"].shape != (1,):
    msg = "Expect only one EOF"
    raise AssertionError(msg)

# pc

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
delta_pc = delta_E_xr * Q(reg_info["m"][0], reg_info["m"][1])
# delta_pc

# %%
pc_extended = delta_pc + pc.sel(year=last_hist_year).data

fig, axes = plt.subplots(nrows=2, figsize=(8, 8))
annual_mean_emissions_emms_units.T.plot(ax=axes[0])
pc_extended.pint.to("dimensionless").plot(ax=axes[1], hue="scenario")
pc.pint.to("dimensionless").plot(ax=axes[1])
for ax in axes:
    sns.move_legend(ax, loc="center left", bbox_to_anchor=(1.05, 0.5))

# %% [markdown]
# ## Interpolate to monthly

# %% [markdown]
# We have to ensure that the scenarios stay harmonised until the last month of the historical data.
# Therefore, we use the same wall control point value for all scenarios.
# As a simple choice, we use the average wall control point across all scenarios.

# %%
fixed_control_point = (
    pc_extended.sel(year=[last_hist_year, last_hist_year + 1]).mean("year").mean("scenario").data.squeeze()
)
# fixed_control_point


# %%
def get_wall_control_points(
    intervals_x: pint.UnitRegistry.Quantity,
    intervals_y: pint.UnitRegistry.Quantity,
    control_points_wall_x: pint.UnitRegistry.Quantity,
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
pc_extended_monthly_l = []

for scenario, sda in pc_extended.groupby("scenario"):
    scenario_monthly = interpolate_annual_mean_to_monthly(
        values=sda.data.m.squeeze(),
        values_units=pc_extended.data.u,
        years=pc_extended["year"].values,
        algorithm=LaiKaplanInterpolator(
            progress_bar=True,
            get_wall_control_points_y_from_interval_ys=get_wall_control_points,
        ),
        unit_registry=openscm_units.unit_registry,
    ).assign_coords(scenario=scenario)
    pc_extended_monthly_l.append(scenario_monthly)

pc_extended_monthly = xr.concat(pc_extended_monthly_l, dim="scenario")
pc_extended_monthly.name = "principal-components-monthly"
pc_extended_monthly.attrs["description"] = "principal component values on a monthly timestep"
pc_extended_monthly

# %%
months_per_year = 12
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
fig, axes = plt.subplots(ncols=2)

for years, ax in (
    (np.arange(last_hist_year + 1, 2500 + 1), axes[0]),
    (np.arange(last_hist_year + 1, last_hist_year + 15), axes[1]),
):
    pc_extended_monthly.sel(time=pc_extended_monthly["time"].dt.year.isin(years)).pint.to("dimensionless").plot(
        ax=ax, alpha=0.6, hue="scenario"
    )
    # convert_year_to_time(pc_extended.sel(year=years)).pint.to("dimensionless").plot.scatter(ax=ax, hue="scenario")

plt.tight_layout()

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ## Prepare output

# %%
if ghg != "c8f18":
    eofs = cmip7_lat_gradient_pieces_ds["eofs"]

else:
    # Calculate the EOF based on the output on ESGF.
    # In general this doesn't work,
    # but for c8f18 it's ok because seasonality is zero
    # and we want a constant latitudinal gradient in the future.
    cmip7_historical_gn_monthly_ds = load_file_from_glob(f"*{ghg}_*nz_175001-*.nc", historical_data_root_dir_p)
    cmip7_historical_gn_monthly = cmip7_historical_gn_monthly_ds[ghg]
    cmip7_historical_gn_monthly["lat"].attrs.pop("units")
    months_per_year = 12
    eofs = (
        (
            cmip7_historical_gn_monthly.sel(
                time=(cmip7_historical_gn_monthly["time"].dt.year == last_hist_year)
                & (cmip7_historical_gn_monthly["time"].dt.month == months_per_year)
            ).pint.quantify(unit_registry=ur)
            - cmip7_historical_gm_monthly.sel(
                time=(cmip7_historical_gm_monthly["time"].dt.year == last_hist_year)
                & (cmip7_historical_gm_monthly["time"].dt.month == months_per_year)
            ).pint.quantify(unit_registry=ur)
        )
        .isel(time=0)
        .assign_coords(eof=0)
        .expand_dims("eof")
    )
    eofs.name = "eofs"

out = xr.merge(
    [
        eofs,
        pc_extended_monthly.pint.to("dimensionless"),
    ],
    combine_attrs="drop_conflicts",
).pint.dequantify()

out

# %%
# Double check this can work
tmp = out["eofs"] * out["principal-components-monthly"]
# tmp

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ## Save

# %% editable=true slideshow={"slide_type": ""}
out_file_p.parent.mkdir(exist_ok=True, parents=True)
out.to_netcdf(out_file_p)
out_file_p

# %%
# Added for the CMIP7 ScenarioMIP GHG manuscript.
# The original run only saves the monthly PCs, we want the annual ones too.
import pandas as pd

manuscript_pcs_l = []
for manuscript_eof, manuscript_pc in (
    (0, pc_extended),
):
    manuscript_pc = manuscript_pc.pint.to("dimensionless").pint.dequantify()
    manuscript_pc = manuscript_pc.sel(year=manuscript_pc["year"] >= last_hist_year)
    manuscript_pc_df = (
        manuscript_pc.rename("value").to_dataframe()[["value"]].reset_index()
    )
    manuscript_pc_df["eof"] = manuscript_eof
    manuscript_pcs_l.append(manuscript_pc_df[["eof", "scenario", "year", "value"]])

manuscript_pcs_file = Path("manuscript-outputs/rerun-markers/latitudinal-gradient/hfc134a_pcs-annual.csv")
manuscript_pcs_file.parent.mkdir(exist_ok=True, parents=True)
pd.concat(manuscript_pcs_l).to_csv(manuscript_pcs_file, index=False)
manuscript_pcs_file
