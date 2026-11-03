"""
Loading of the published (input4MIPs) concentrations, for any gas

Everything comes from ESGF, via the same database and local cache
as the historical manuscript's comparison data
(see [local.historical_ghg_forcing_for_cmip7.comparison_data][]).
It is downloaded the first time and re-used from disk after that.
"""

from __future__ import annotations

from functools import cache

import numpy as np
import openscm_units
import pandas as pd
import xarray as xr

from local.data_loading import fetch_and_load_ghg_dataset
from local.esgf.search.search_query import KnownIndexNode
from local.historical_ghg_forcing_for_cmip7.comparison_data import (
    ESGF_LOCAL_DATA_ROOT_DIR,
    get_esgf_engine,
)
from local.scenario_ghg_forcing_for_cmip7.scenarios import (
    CMIP6_HISTORY_SOURCE_ID,
    CMIP7_HISTORY_SOURCE_ID,
    Scenario,
)

HEMISPHERES = ("Northern hemisphere", "Southern hemisphere")
"""Names of the hemispheres, as named in the figures"""


@cache
def load_esgf_data(
    gas: str, cmip_era: str, source_id: str, grid: str, time_sampling: str
) -> xr.DataArray:
    """
    Load one published dataset

    Parameters
    ----------
    gas
        Gas of interest (CMIP7 name, also for CMIP6 data)

    cmip_era
        CMIP era of the dataset

    source_id
        Source ID of the dataset

    grid
        Grid of the dataset ("gm" or "gnz", CMIP7 names, also for CMIP6 data)

    time_sampling
        Time sampling of the dataset ("yr" or "mon")

    Returns
    -------
    :
        The dataset's data, sorted in time
        (its files aren't necessarily loaded in time order)
    """
    ds = fetch_and_load_ghg_dataset(
        local_data_root_dir=ESGF_LOCAL_DATA_ROOT_DIR,
        index_node=KnownIndexNode.ORNL,
        ghg=gas,
        grid=grid,
        time_sampling=time_sampling,
        cmip_era=cmip_era,
        source_id=source_id,
        engine=get_esgf_engine(),
    )

    return ds[gas].sortby("time").compute()


def to_decimal_year(da: xr.DataArray, time_sampling: str) -> xr.DataArray:
    """
    Put data's time axis in decimal years, at the middle of each time step

    The same way the historical manuscript's figures place their data.

    Parameters
    ----------
    da
        Data, with a cftime `time` axis

    time_sampling
        Time sampling of `da` ("yr" or "mon")

    Returns
    -------
    :
        `da`, with `time` as decimal years
    """
    if time_sampling == "yr":
        time = np.array([t.year + 0.5 for t in da["time"].values])
    elif time_sampling == "mon":
        time = np.array([t.year + (t.month - 0.5) / 12.0 for t in da["time"].values])
    else:
        raise NotImplementedError(time_sampling)

    return da.assign_coords(time=time)


def convert_units(da: xr.DataArray, units: str) -> xr.DataArray:
    """
    Convert data to the given units

    Parameters
    ----------
    da
        Data, with a `units` attribute

    units
        Units to convert to

    Returns
    -------
    :
        `da`, in `units`
    """
    if da.attrs["units"] == units:
        return da

    factor = openscm_units.unit_registry.Quantity(1.0, da.attrs["units"]).to(units).m
    res = da * factor
    res.attrs = {**da.attrs, "units": units}

    return res


def load_concatenated(
    gas: str,
    cmip_era: str,
    source_ids: tuple[str, ...],
    grid: str,
    time_sampling: str,
) -> xr.DataArray:
    """
    Load the datasets which together make up one timeseries, joined in time

    Parameters
    ----------
    gas
        Gas of interest

    cmip_era
        CMIP era of the datasets

    source_ids
        Source IDs of the datasets, in time order

    grid
        Grid of the datasets

    time_sampling
        Time sampling of the datasets

    Returns
    -------
    :
        The datasets' data, joined in time,
        with `time` as decimal years and in the first dataset's units
    """
    pieces = [
        to_decimal_year(
            load_esgf_data(gas, cmip_era, source_id, grid, time_sampling),
            time_sampling,
        )
        for source_id in source_ids
    ]
    units = pieces[0].attrs["units"]
    res = xr.concat(
        [convert_units(piece, units) for piece in pieces],
        dim="time",
        coords="minimal",
        compat="override",
        combine_attrs="override",
    )
    if not (np.diff(res["time"].values) > 0).all():
        msg = f"{source_ids} overlap or are out of order"
        raise AssertionError(msg)

    return res


def get_scenario_yearly_global_mean(gas: str, scenario: Scenario) -> xr.DataArray:
    """
    Get a scenario's yearly global-mean

    Parameters
    ----------
    gas
        Gas of interest

    scenario
        Scenario of interest

    Returns
    -------
    :
        Yearly global-mean, with `time` as decimal years
    """
    return load_concatenated(gas, scenario.cmip_era, scenario.source_ids, "gm", "yr")


def get_history_yearly_global_mean(gas: str, cmip_era: str) -> xr.DataArray:
    """
    Get the history's yearly global-mean

    Parameters
    ----------
    gas
        Gas of interest

    cmip_era
        CMIP era of interest

    Returns
    -------
    :
        Yearly global-mean, with `time` as decimal years
    """
    source_id = {
        "CMIP7": CMIP7_HISTORY_SOURCE_ID,
        "CMIP6": CMIP6_HISTORY_SOURCE_ID,
    }[cmip_era]

    return load_concatenated(gas, cmip_era, (source_id,), "gm", "yr")


def get_monthly_native_resolution(
    gas: str, source_ids: tuple[str, ...]
) -> xr.DataArray:
    """
    Get CMIP7 monthly data at its native resolution (15 degree zonal means)

    Parameters
    ----------
    gas
        Gas of interest

    source_ids
        Source IDs of the datasets to join, in time order

    Returns
    -------
    :
        Monthly data with `time` (as decimal years) and `lat` dimensions
    """
    return load_concatenated(gas, "CMIP7", source_ids, "gnz", "mon").drop_vars(
        ["lat_bnds", "time_bnds"], errors="ignore"
    )


def get_spatial_means(native_resolution: xr.DataArray) -> xr.DataArray:
    """
    Get the global- and hemispheric-means of data at native resolution

    Weighted by the cosine of latitude,
    which reproduces the published global-means.

    Parameters
    ----------
    native_resolution
        Data with a `lat` dimension

    Returns
    -------
    :
        Data with a `region` dimension ("Global" and [HEMISPHERES][])
        in place of `lat`
    """
    weights = np.cos(np.deg2rad(native_resolution["lat"]))
    north = native_resolution["lat"] > 0.0
    res = xr.concat(
        [
            native_resolution.weighted(weights).mean("lat"),
            native_resolution.where(north).weighted(weights).mean("lat"),
            native_resolution.where(~north).weighted(weights).mean("lat"),
        ],
        dim=pd.Index(["Global", *HEMISPHERES], name="region"),
    )
    res.attrs = native_resolution.attrs

    return res
