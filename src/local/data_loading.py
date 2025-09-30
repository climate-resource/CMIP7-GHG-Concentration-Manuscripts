"""
Data loading

High-level helpers to hide the challenges of grabbing data from different sources
"""

from __future__ import annotations

from pathlib import Path
from typing import Callable

import ncdata
import sqlalchemy.engine.base
import xarray as xr
from sqlmodel import Session, select

from local.esgf.esgf_dataset_collection import ESGFDatasetCollection
from local.esgf.models import ESGFDataset, ESGFDatasetDB
from local.esgf.search import SearchQuery
from local.xarray_loading import load_xarray_from_esgf_dataset

CMIP7_TO_CMIP6_VARIABLE_MAP = {
    "co2": "mole_fraction_of_carbon_dioxide_in_air",
    "ch4": "mole_fraction_of_methane_in_air",
    "n2o": "mole_fraction_of_nitrous_oxide_in_air",
    "c2f6": "mole_fraction_of_c2f6_in_air",
    "c3f8": "mole_fraction_of_c3f8_in_air",
    "c4f10": "mole_fraction_of_c4f10_in_air",
    "c5f12": "mole_fraction_of_c5f12_in_air",
    "c6f14": "mole_fraction_of_c6f14_in_air",
    "c7f16": "mole_fraction_of_c7f16_in_air",
    "c8f18": "mole_fraction_of_c8f18_in_air",
    "cc4f8": "mole_fraction_of_c_c4f8_in_air",
    "ccl4": "mole_fraction_of_carbon_tetrachloride_in_air",
    "cf4": "mole_fraction_of_cf4_in_air",
    "cfc11": "mole_fraction_of_cfc11_in_air",
    "cfc113": "mole_fraction_of_cfc113_in_air",
    "cfc114": "mole_fraction_of_cfc114_in_air",
    "cfc115": "mole_fraction_of_cfc115_in_air",
    "cfc12": "mole_fraction_of_cfc12_in_air",
    "ch2cl2": "mole_fraction_of_ch2cl2_in_air",
    "ch3br": "mole_fraction_of_methyl_bromide_in_air",
    "ch3ccl3": "mole_fraction_of_ch3ccl3_in_air",
    "ch3cl": "mole_fraction_of_methyl_chloride_in_air",
    "chcl3": "mole_fraction_of_chcl3_in_air",
    "halon1211": "mole_fraction_of_halon1211_in_air",
    "halon1301": "mole_fraction_of_halon1301_in_air",
    "halon2402": "mole_fraction_of_halon2402_in_air",
    "hcfc141b": "mole_fraction_of_hcfc141b_in_air",
    "hcfc142b": "mole_fraction_of_hcfc142b_in_air",
    "hcfc22": "mole_fraction_of_hcfc22_in_air",
    "hfc125": "mole_fraction_of_hfc125_in_air",
    "hfc134a": "mole_fraction_of_hfc134a_in_air",
    "hfc143a": "mole_fraction_of_hfc143a_in_air",
    "hfc152a": "mole_fraction_of_hfc152a_in_air",
    "hfc227ea": "mole_fraction_of_hfc227ea_in_air",
    "hfc23": "mole_fraction_of_hfc23_in_air",
    "hfc236fa": "mole_fraction_of_hfc236fa_in_air",
    "hfc245fa": "mole_fraction_of_hfc245fa_in_air",
    "hfc32": "mole_fraction_of_hfc32_in_air",
    "hfc365mfc": "mole_fraction_of_hfc365mfc_in_air",
    "hfc4310mee": "mole_fraction_of_hfc4310mee_in_air",
    "nf3": "mole_fraction_of_nf3_in_air",
    "sf6": "mole_fraction_of_sf6_in_air",
    "so2f2": "mole_fraction_of_so2f2_in_air",
    "cfc11eq": "mole_fraction_of_cfc11eq_in_air",
    "cfc12eq": "mole_fraction_of_cfc12eq_in_air",
    "hfc134aeq": "mole_fraction_of_hfc134aeq_in_air",
}


def fix_year_zero_nonsense(
    ncdatas: list[ncdata.NcData], esgf_dataset: ESGFDataset, file_paths: list[Path]
) -> list[ncdata.NcData]:
    for ncd in ncdatas:
        if ncd.variables["time"].attributes["units"].value == "days since 0-1-1" and (
            ncd.variables["time"].attributes["calendar"].value == "gregorian"
        ):
            # # Time units broken - fix here
            # Drop out year 0 values (year 0 doesn't exist)
            ncd.variables[
                str(ncd.attributes["variable_id"].value)
            ].data = ncd.variables[str(ncd.attributes["variable_id"].value)].data[1:, :]
            ncd.variables["time"].data = ncd.variables["time"].data[1:]
            ncd.variables["time_bnds"].data = ncd.variables["time_bnds"].data[1:, :]
            # Use units that exist
            ncd.variables["time"].attributes["units"].value = "days since 1-1-1"
            # Align values with new units
            ncd.variables["time"].data -= 365
            ncd.variables["time_bnds"].data -= 365

    return ncdatas


def load_cmip_ghg_ds(esgf_dataset: ESGFDataset) -> xr.Dataset:
    return load_xarray_from_esgf_dataset(
        esgf_dataset=esgf_dataset,
        pre_to_xarray=fix_year_zero_nonsense,
    )


def fetch_and_load_ghg_dataset(  # noqa: PLR0913
    local_data_root_dir: Path,
    ghg: str,
    grid: str,
    time_sampling: str,
    cmip_era: str,
    source_id: str,
    index_node: str,
    engine: sqlalchemy.engine.base.Engine,
    load_xr_dataset: Callable[[ESGFDataset], xr.Dataset] = load_cmip_ghg_ds,
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

    load_xr_dataset
        Function to use to load an [xr.Dataset][xarray.Dataset]

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

        if cmip_era == "CMIP6":
            # Make sure that next step will hold
            if len(esgf_dataset_collection.esgf_datasets) != 1:
                raise AssertionError(esgf_dataset_collection)
            # Overwrite the query result value with the one we want
            # TODO: better logic around this
            esgf_dataset_collection.esgf_datasets[0].variable = ghg

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
    esgf_dataset_collection.ensure_all_files_available_locally()

    # Finally, load the local files and return the xr.Dataset
    # TODO: other loading methods e.g. load into xr.Datatree
    if len(esgf_dataset_collection.esgf_datasets) != 1:
        raise AssertionError(esgf_dataset_collection)
    esgf_dataset = esgf_dataset_collection.esgf_datasets[0]

    res = load_xr_dataset(esgf_dataset)

    # TODO: add support for realigning grid
    # to match CMIP7 style grids as needed
    return res
