"""
Getting the original run's files which aren't in its Zenodo archive

The original run's Zenodo archive (see [local.cmip_ghg_generation][])
doesn't contain every intermediate file the run created.
We keep a copy of the ones we need in this repository,
in the directory given by [get_zenodo_missing_dir][].
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from local.paths import DATA_RAW_DIR

CH4_ICE_CORE_FILES = {
    "epica": "epica_with_location.csv",
    "law-dome": "law-dome_ch4_smoothed_median.csv",
    "neem": "neem_with_location.csv",
}
"""File which holds each ice core's CH4 data (with its location)"""


def get_zenodo_missing_dir(data_raw_dir: Path = DATA_RAW_DIR) -> Path:
    """
    Get the directory in which we keep the files missing from Zenodo

    Parameters
    ----------
    data_raw_dir
        Raw data directory

    Returns
    -------
    :
        Directory in which we keep the files missing from Zenodo
    """
    return data_raw_dir / "historical-ghg-forcing-for-cmip7" / "zenodo-missing"


def get_ch4_ice_core_file(source: str, data_raw_dir: Path = DATA_RAW_DIR) -> Path:
    """
    Get the file which holds an ice core's CH4 data

    Parameters
    ----------
    source
        Ice core of interest (a key of [CH4_ICE_CORE_FILES][])

    data_raw_dir
        Raw data directory

    Returns
    -------
    :
        File which holds `source`'s CH4 data
    """
    return get_zenodo_missing_dir(data_raw_dir) / CH4_ICE_CORE_FILES[source]


def load_ch4_ice_core(source: str, data_raw_dir: Path = DATA_RAW_DIR) -> pd.DataFrame:
    """
    Load an ice core's CH4 data

    Parameters
    ----------
    source
        Ice core of interest (a key of [CH4_ICE_CORE_FILES][])

    data_raw_dir
        Raw data directory

    Returns
    -------
    :
        `source`'s CH4 data, with its location
    """
    return pd.read_csv(get_ch4_ice_core_file(source, data_raw_dir))


def get_ice_core_latitude(ice_core: pd.DataFrame, source: str) -> float:
    """
    Get an ice core's latitude

    Parameters
    ----------
    ice_core
        The ice core's data (e.g. from [load_ch4_ice_core][])

    source
        Name of the ice core (only used in the error message)

    Returns
    -------
    :
        The ice core's latitude

    Raises
    ------
    AssertionError
        `ice_core` has more than one latitude
    """
    lat_l = ice_core["latitude"].unique()
    if len(lat_l) > 1:
        msg = f"Expected a single latitude for {source}, found {lat_l}"
        raise AssertionError(msg)

    return float(lat_l[0])
