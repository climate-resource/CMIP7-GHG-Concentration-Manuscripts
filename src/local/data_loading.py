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

CMIP7_TO_CMIP6_VARIABLE_MAP = {
    "co2": "mole-fraction-of-carbon-dioxide-in-air",
    "ch4": "mole-fraction-of-methane-in-air",
    "n2o": "mole-fraction-of-nitrous-oxide-in-air",
    "c2f6": "mole-fraction-of-c2f6-in-air",
    "c3f8": "mole-fraction-of-c3f8-in-air",
    "c4f10": "mole-fraction-of-c4f10-in-air",
    "c5f12": "mole-fraction-of-c5f12-in-air",
    "c6f14": "mole-fraction-of-c6f14-in-air",
    "c7f16": "mole-fraction-of-c7f16-in-air",
    "c8f18": "mole-fraction-of-c8f18-in-air",
    "cc4f8": "mole-fraction-of-c-c4f8-in-air",
    "ccl4": "mole-fraction-of-carbon-tetrachloride-in-air",
    "cf4": "mole-fraction-of-cf4-in-air",
    "cfc11": "mole-fraction-of-cfc11-in-air",
    "cfc113": "mole-fraction-of-cfc113-in-air",
    "cfc114": "mole-fraction-of-cfc114-in-air",
    "cfc115": "mole-fraction-of-cfc115-in-air",
    "cfc12": "mole-fraction-of-cfc12-in-air",
    "ch2cl2": "mole-fraction-of-ch2cl2-in-air",
    "ch3br": "mole-fraction-of-methyl-bromide-in-air",
    "ch3ccl3": "mole-fraction-of-ch3ccl3-in-air",
    "ch3cl": "mole-fraction-of-methyl-chloride-in-air",
    "chcl3": "mole-fraction-of-chcl3-in-air",
    "halon1211": "mole-fraction-of-halon1211-in-air",
    "halon1301": "mole-fraction-of-halon1301-in-air",
    "halon2402": "mole-fraction-of-halon2402-in-air",
    "hcfc141b": "mole-fraction-of-hcfc141b-in-air",
    "hcfc142b": "mole-fraction-of-hcfc142b-in-air",
    "hcfc22": "mole-fraction-of-hcfc22-in-air",
    "hfc125": "mole-fraction-of-hfc125-in-air",
    "hfc134a": "mole-fraction-of-hfc134a-in-air",
    "hfc143a": "mole-fraction-of-hfc143a-in-air",
    "hfc152a": "mole-fraction-of-hfc152a-in-air",
    "hfc227ea": "mole-fraction-of-hfc227ea-in-air",
    "hfc23": "mole-fraction-of-hfc23-in-air",
    "hfc236fa": "mole-fraction-of-hfc236fa-in-air",
    "hfc245fa": "mole-fraction-of-hfc245fa-in-air",
    "hfc32": "mole-fraction-of-hfc32-in-air",
    "hfc365mfc": "mole-fraction-of-hfc365mfc-in-air",
    "hfc4310mee": "mole-fraction-of-hfc4310mee-in-air",
    "nf3": "mole-fraction-of-nf3-in-air",
    "sf6": "mole-fraction-of-sf6-in-air",
    "so2f2": "mole-fraction-of-so2f2-in-air",
    "cfc11eq": "mole-fraction-of-cfc11eq-in-air",
    "cfc12eq": "mole-fraction-of-cfc12eq-in-air",
    "hfc134aeq": "mole-fraction-of-hfc134aeq-in-air",
}


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

        if results:
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
        if cmip_era == "CMIP6":
            # Translate names
            variable = CMIP7_TO_CMIP6_VARIABLE_MAP[ghg]

        else:
            variable = ghg

        # Query the DB and add the results to the DB
        query = SearchQuery(
            project="input4MIPs",
            variable=variable,
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
