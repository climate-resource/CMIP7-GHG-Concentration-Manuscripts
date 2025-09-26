"""
Data loading

High-level helpers to hide the challenges of grabbing data from different sources
"""

from __future__ import annotations

import sqlalchemy.engine.base
import xarray as xr
from sqlmodel import Session, select

from local.esgf.models import ESGFDataset, ESGFDatasetDB
from local.esgf.search import SearchQuery


def fetch_and_load_ghg_dataset(  # noqa: PLR0913
    ghg: str,
    grid: str,
    time_sampling: str,
    cmip_era: str,
    source_id: str,
    index_node: str,
    engine: sqlalchemy.engine.base.Engine,
) -> xr.Dataset:
    """
    Fetch (if needed) and load a greenhouse gas concentration dataset

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

    index_node
        Index node to use for searching

    engine
        Database engine for saving results

    Returns
    -------
    :
        Loaded data
    """
    with Session(engine) as session:
        statement = select(ESGFDatasetDB).where(
            # I Wonder if we can do `where(query)` or something
            # to be less repetitive
            ESGFDatasetDB.variable == ghg,
            ESGFDatasetDB.grid == grid,
            ESGFDatasetDB.time_sampling == time_sampling,
            ESGFDatasetDB.cmip_era == cmip_era,
            ESGFDatasetDB.source_id == source_id,
        )
        result_exec = session.exec(statement)

        result = result_exec.one_or_none()

        if result is not None:
            esgf_dataset: ESGFDataset | None = ESGFDataset.model_validate(result)

        else:
            esgf_dataset = None

    if esgf_dataset is None:
        query = SearchQuery(
            project="input4MIPs",
            variable=ghg,
            grid=grid,
            time_sampling=time_sampling,
            cmip_era=cmip_era,
            source_id=source_id,
        )
        # Check if results already in the database
        # Optionally check for new results from ESGF
        # (don't implement for now)
        # If no results in database, grab them from ESGF
        # Then save into database
        esgf_datasets = query.get_results(index_node=index_node)
        if len(esgf_datasets) != 1:
            raise AssertionError(esgf_datasets)

        esgf_dataset = esgf_datasets[0]

        # TODO: add something that allows you to check whether new results
        # differ from existing if the user wants to check this

        with Session(engine) as session:
            session.add(esgf_dataset.to_db_model())
            session.commit()

    assert False, "Up to here"
    local_files = [
        esgff.to_local_file(local_root_dir=root_data_dir)
        for esgff in esgf_dataset.esgf_files
    ]
    to_download = [
        lf
        for lf in local_files
        if not lf.available_locally()
        # TODO: Logic to put in available_locally
        # if f.local_path is None or not f.local_path.exists()
    ]
    # (only downloading time slices needed would be smarter, but ok).
    if to_download:
        local_files_downloaded = download_files(
            to_download,
            # tqdm stuff here
            # can use logging in here too
        )
        # TODO: Update the LocalFile table to record downloaded status
        local_files = [
            *local_files_downloaded,
            *(lf for lf in local_files if lf.available_locally()),
        ]

    # Finally, load the local files and return the xr.Dataset
    res = xr.open_mfdataset([lf.local_path for lf in local_files])

    return res
