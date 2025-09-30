"""
Data loading

High-level helpers to hide the challenges of grabbing data from different sources
"""

from __future__ import annotations

from pathlib import Path

import sqlalchemy.engine.base
import xarray as xr
from sqlmodel import Session, select

from local.esgf.esgf_dataset_collection import ESGFDatasetCollection
from local.esgf.models import ESGFDataset, ESGFDatasetDB
from local.esgf.search import SearchQuery


def fetch_and_load_ghg_dataset(  # noqa: PLR0913
    local_data_root_dir: Path,
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
    local_data_root_dir
        Path to the root directory in which local data should be downloaded

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
    # Check for values already in the database
    with Session(engine) as session:
        statement = select(ESGFDatasetDB).where(
            # I Wonder if we can do `where(query)`
            # where `query` is an object or something
            # to be less repetitive
            ESGFDatasetDB.variable == ghg,
            ESGFDatasetDB.grid == grid,
            ESGFDatasetDB.time_sampling == time_sampling,
            ESGFDatasetDB.cmip_era == cmip_era,
            ESGFDatasetDB.source_id == source_id,
        )
        result_exec = session.exec(statement)

        results = result_exec.all()

        if results is not None:
            esgf_dataset_collection: ESGFDatasetCollection | None = (
                ESGFDatasetCollection(
                    esgf_datasets=tuple(ESGFDataset.model_validate(r) for r in results)
                )
            )

        else:
            esgf_dataset_collection = None

    # TODO: optionally allow for checking for new results from ESGF,
    # even if we already have search results in the DB

    if esgf_dataset_collection is None:
        # Query the DB and add the results to the DB
        query = SearchQuery(
            project="input4MIPs",
            variable=ghg,
            grid=grid,
            time_sampling=time_sampling,
            cmip_era=cmip_era,
            source_id=source_id,
        )
        esgf_dataset_collection = query.get_results(index_node=index_node)

        # TODO: add something that allows you to check whether new results
        # differ from existing if the user wants to check this

        esgf_dataset_collection.set_local_files_root_dir(local_data_root_dir)

        with Session(engine) as session:
            session.add_all(
                esgf_dataset.to_db_model()
                for esgf_dataset in esgf_dataset_collection.esgf_datasets
            )
            session.commit()

    # Download
    local_paths = esgf_dataset_collection.ensure_all_files_available_locally()

    # Finally, load the local files and return the xr.Dataset

    # TODO: abstract into load_xarray_from_db_dataset
    # so we can do pre- and post-processing
    # based on metadata in dataset and handle the year 0 issue for CMIP6

    # Make sure that we can load our paths into a single dataset
    if len(esgf_dataset_collection.esgf_datasets) != 1:
        raise AssertionError(esgf_dataset_collection)
    esgf_dataset = esgf_dataset_collection.esgf_datasets[0]
    # TODO: other loading methods e.g. load into xr.Datatree
    res = xr.open_mfdataset(
        local_paths, use_cftime=True, data_vars=None, compat="no_conflicts"
    )

    return res
