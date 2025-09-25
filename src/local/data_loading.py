"""
Data loading

High-level helpers to hide the challenges of grabbing data from different sources
"""

from __future__ import annotations

import xarray as xr

from local.esgf.search import SearchQuery


def fetch_and_load_ghg_file(  # noqa: PLR0913
    ghg: str,
    grid: str,
    time_sampling: str,
    cmip_era: str,
    source_id: str,
    index_node: str,
) -> xr.Dataset:
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

    cmip_era
        CMIP era from which the data should be retrieved

    source_id
        Source ID for which to retrieve data

    Returns
    -------
    :
        Loaded data
    """
    query = SearchQuery(
        project="input4MIPs",
        variable=ghg,
        grid=grid,
        time_sampling=time_sampling,
        cmip_era=cmip_era,
        source_id=source_id,
    )
    search_results = query.get_results(index_node=index_node)
    breakpoint()
    save_search_query_and_results_to_db(
        # Database into which to save search results
    )
    # Could also have:
    # - `compare_search_query_and_results_to_db`
