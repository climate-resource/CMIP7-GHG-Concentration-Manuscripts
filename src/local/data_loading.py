"""
Data loading

High-level helpers to hide the challenges of grabbing data from different sources
"""

from __future__ import annotations

import xarray as xr


def fetch_and_load_ghg_file(ghg: str, grid: str, frequency: str) -> xr.Dataset:
    """
    Fetch (if needed) and load a greenhouse gas concentration file

    This also deals with unifying variable names, metadata, units etc.
    across the CMIP6 and CMIP7 data

    Parameters
    ----------
    ghg
        Greenhouse gas for which to fetch data

    grid
        Grid on which to fetch the data

    time_sampling
        Time sampling to fetch

    Returns
    -------
    :
        Loaded data
    """
    # import pdb
    #
    # pdb.set_trace()
