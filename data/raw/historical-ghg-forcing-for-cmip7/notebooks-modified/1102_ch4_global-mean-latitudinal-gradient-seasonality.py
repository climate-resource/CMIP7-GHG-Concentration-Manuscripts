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
# # CH$_4$ - calculate global- annual-mean, latitudinal gradient and seasonality
#
# Calculate global-mean, latitudinal gradient and seasonality.

# %% [markdown]
# ## Imports

# %%
from pathlib import Path

import cf_xarray.units
import matplotlib.axes
import matplotlib.pyplot as plt
import pint_xarray
import xarray as xr
from pydoit_nb.config_handling import get_config_for_step_id

import local.binned_data_interpolation
import local.binning
import local.latitudinal_gradient
import local.raw_data_processing
import local.seasonality
import local.xarray_space
import local.xarray_time
from local.config import load_config_from_file

# %%
cf_xarray.units.units.define("ppm = 1 / 1000000")
cf_xarray.units.units.define("ppb = ppm / 1000")

pint_xarray.accessors.default_registry = pint_xarray.setup_registry(cf_xarray.units.units)

# %% [markdown]
# ## Define branch this notebook belongs to

# %% editable=true slideshow={"slide_type": ""}
step: str = "calculate_ch4_monthly_fifteen_degree_pieces"

# %% [markdown]
# ## Parameters

# %% editable=true slideshow={"slide_type": ""} tags=["parameters"]
config_file: str = "../../dev-config-absolute.yaml"  # config file
step_config_id: str = "only"  # config ID to select for this branch

# %% tags=["injected-parameters"]
# Parameters
config_file = "v1.0.0-config-raw.yaml"
step_config_id = "only"

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ## Load config

# %% editable=true slideshow={"slide_type": ""}
config = load_config_from_file(Path(config_file))
config_step = get_config_for_step_id(config=config, step=step, step_config_id=step_config_id)


# %% [markdown]
# ## Action

# %% [markdown]
# ### Load data

# %%
interpolated_spatial: xr.DataArray = xr.load_dataarray(  # type: ignore
    config_step.observational_network_interpolated_file
).pint.quantify()
interpolated_spatial

# %% [markdown]
# ### Drop out any years that have nan
#
# These break everything later on

# %%
interpolated_spatial_nan_free = local.xarray_time.convert_year_month_to_time(
    local.xarray_time.convert_time_to_year_month(interpolated_spatial).dropna("year")
)

# %% [markdown]
# ### Longitudinal-mean

# %%
lon_mean = interpolated_spatial_nan_free.mean(dim="lon")
lon_mean.plot(hue="lat")  # type: ignore

# %% [markdown]
# ### Global-, annual-mean

# %%
global_mean = local.xarray_space.calculate_global_mean_from_lon_mean(lon_mean)
global_mean.plot()  # type: ignore

# %%
global_annual_mean = global_mean.groupby("time.year").mean()
global_annual_mean.plot()  # type: ignore

# %% [markdown]
# ### Latitudinal gradient

# %%
(
    lat_residuals_annual_mean,
    full_eofs_pcs,
) = local.latitudinal_gradient.calculate_eofs_pcs(lon_mean, global_mean)
full_eofs_pcs

# %%
fig, ax = plt.subplots()
n_eofs_to_show = 3

for i in range(n_eofs_to_show):
    ax.plot(
        full_eofs_pcs["principal-components"].sel(eof=i).data.m,
        label=f"Principal component {i}",
    )

ax.legend(loc="center left", bbox_to_anchor=(1.05, 0.5))

# %%
fig, ax = plt.subplots()

for i in range(3):
    ax.plot(
        full_eofs_pcs["eofs"].sel(eof=i).data.m,
        full_eofs_pcs["lat"].data,
        label=f"EOF {i}",
        zorder=2 - i / 10,
    )

ax.legend(loc="center left", bbox_to_anchor=(1.05, 0.5))

# %%
fig, ax = plt.subplots()

for i in range(3):
    ax.plot(
        full_eofs_pcs["principal-components"].sel(eof=i).isel(year=1) @ full_eofs_pcs["eofs"].sel(eof=i),
        full_eofs_pcs["lat"],
        label=f"EOF {i}",
        zorder=2 - i / 10,
    )

ax.legend(loc="center left", bbox_to_anchor=(1.05, 0.5))

# %%
eofs_pcs = full_eofs_pcs.sel(eof=range(config_step.lat_gradient_n_eofs_to_use))
eofs_pcs

# %%
latitudinal_anomaly_from_eofs = eofs_pcs["principal-components"] @ eofs_pcs["eofs"]

for year in latitudinal_anomaly_from_eofs["year"]:
    if year % 5:
        continue

    fig, axes = plt.subplots(nrows=3, sharex=True, sharey=True)
    if isinstance(axes, matplotlib.axes.Axes):
        raise TypeError(type(axes))

    selected = lat_residuals_annual_mean.sel(year=year)
    axes[0].plot(selected.data.m, lat_residuals_annual_mean.lat)
    axes[0].set_title("Full anomaly")

    selected_eofs = latitudinal_anomaly_from_eofs.sel(year=year)
    axes[1].plot(selected_eofs.data.m, latitudinal_anomaly_from_eofs.lat)
    axes[1].set_title("Anomaly from EOFs")

    axes[2].plot(selected.data.m - selected_eofs.data.m, lat_residuals_annual_mean.lat)
    axes[2].set_title("Difference")

    axes[0].set_ylim([-90, 90])

    plt.suptitle(str(int(year)))
    plt.tight_layout()
    plt.show()

# %% [markdown]
# ### Seasonality

# %%
(
    seasonality,
    relative_seasonality,
    lon_mean_ym_monthly_anomalies,
) = local.seasonality.calculate_seasonality(
    lon_mean=lon_mean,
    global_mean=global_mean,
)

# %%
lon_mean_ym_monthly_anomalies.plot.line(hue="lat", col="year", col_wrap=3)

# %%
seasonality.plot.line(hue="lat")

# %%
relative_seasonality.plot.line(hue="lat")

# %% [markdown]
# ### Save

# %%
config_step.observational_network_global_annual_mean_file.parent.mkdir(exist_ok=True, parents=True)
global_annual_mean.pint.dequantify().to_netcdf(config_step.observational_network_global_annual_mean_file)
global_annual_mean

# %%
config_step.observational_network_latitudinal_gradient_eofs_file.parent.mkdir(exist_ok=True, parents=True)
eofs_pcs.pint.dequantify().to_netcdf(config_step.observational_network_latitudinal_gradient_eofs_file)
eofs_pcs

# %%
# Use relative seasonality for CH4
config_step.observational_network_seasonality_file.parent.mkdir(exist_ok=True, parents=True)
relative_seasonality.pint.dequantify().to_netcdf(config_step.observational_network_seasonality_file)
relative_seasonality

# %%
# Added for the CMIP7 GHG manuscript.
# The original run only ever saved the EOFs it keeps,
# but we want to show how much of the variance every EOF explains,
# so that keeping only the first few can be justified.
from pathlib import Path

import numpy as np
import pandas as pd

# The EOFs are the right singular vectors of the residuals, so they are
# orthonormal, which makes the principal components uncorrelated:
# the principal components' cross-product is `diag(D) ** 2`,
# i.e. the singular values squared with nothing off the diagonal.
# Bound to a name of our own, so nothing here can shadow the notebook's.
manuscript_eofs_pcs = full_eofs_pcs
pcs = manuscript_eofs_pcs["principal-components"].transpose("year", "eof").data.m
singular_values_squared = pcs.T @ pcs

off_diagonal = singular_values_squared - np.diag(np.diag(singular_values_squared))
if not np.allclose(off_diagonal, 0.0, atol=1e-10 * np.trace(singular_values_squared)):
    msg = "The principal components are not uncorrelated, so this is not an SVD"
    raise AssertionError(msg)

variance_explained = np.diag(singular_values_squared) / np.trace(
    singular_values_squared
)

full_eofs_pcs_out_file = Path("manuscript-outputs/ch4_lat-gradient-full-eofs-pcs.nc")
full_eofs_pcs_out_file.parent.mkdir(exist_ok=True, parents=True)
to_save = manuscript_eofs_pcs.pint.dequantify()
# A decomposition can be over a stacked dimension -- CO2's seasonality change
# is over latitude and month together -- and netCDF cannot hold a MultiIndex.
# Unstacking it into its levels keeps everything, in a form netCDF can write.
multi_indexes = [
    name
    for name, index in to_save.indexes.items()
    if isinstance(index, pd.MultiIndex)
]
if multi_indexes:
    to_save = to_save.reset_index(multi_indexes)

to_save.to_netcdf(full_eofs_pcs_out_file)

manuscript_out_file = Path("manuscript-outputs/ch4_lat-gradient-variance-explained.csv")
manuscript_out_file.parent.mkdir(exist_ok=True, parents=True)
pd.DataFrame(
    {
        "eof": manuscript_eofs_pcs["eof"].values,
        "variance_explained_fraction": variance_explained,
    }
).to_csv(manuscript_out_file, index=False)
manuscript_out_file
