"""
Data we compare our output against in the results figures

Two kinds of comparison data go into those figures.

- Datasets with spatial information: measurements at a site (an ice core, a station)
  or a network of them. Each value has a latitude,
  so it can be read against our output at the same latitude,
  and it is drawn coloured by that latitude,
  as the observation network is.
- Datasets which are 'just' global-means, e.g. the CMIP6 concentrations.
  These have no latitude, so they are read against our global-mean.

Both are carried in a [ComparisonTimeseries][], the difference being
whether its data has a latitude column.
"""

from __future__ import annotations

import dataclasses
import urllib.request
from dataclasses import dataclass, field
from functools import cache
from pathlib import Path

import numpy as np
import openscm_units
import pandas as pd
import pint
import sqlalchemy.engine.base
import xarray as xr
from loguru import logger

from local.data_loading import fetch_and_load_ghg_dataset, fix_broken_calendar_spec
from local.esgf.db_helpers import create_all_tables, get_sqlite_engine
from local.esgf.models import ESGFDataset
from local.esgf.search.search_query import KnownIndexNode
from local.historical_ghg_forcing_for_cmip7.plotting import OKABE_ITO
from local.paths import DATA_RAW_DIR, REPO_ROOT
from local.xarray_loading import load_xarray_from_esgf_dataset

Q = openscm_units.unit_registry.Quantity

LATITUDE_COLUMN = "latitude"
"""Column which holds the latitude of each value of a spatial comparison dataset"""

TIME_COLUMN = "time"
"""Column which holds the time of each value, as a decimal year"""

VALUE_COLUMN = "value"
"""Column which holds each value"""


@dataclass(frozen=True)
class ComparisonTimeseries:
    """
    A dataset to compare our output against
    """

    label: str
    """Name of the dataset, as the legend shows it"""

    data: pd.DataFrame
    """
    The dataset's values

    Must have a [TIME_COLUMN][] (decimal year) and a [VALUE_COLUMN][].
    If it also has a [LATITUDE_COLUMN][], the dataset has spatial information
    and is drawn as points coloured by latitude.
    Otherwise it is a global-mean and is drawn as a line.
    """

    units: str
    """Units of the values"""

    marker: str = "o"
    """Marker to draw a spatial dataset's points with"""

    spatial_as_line: bool = False
    """
    Whether to draw a spatial dataset as a line rather than as points

    The line is still coloured by latitude.
    For dense records (e.g. a station's monthly means),
    whose points would run into one another.
    Ignored for global-mean datasets, which are always drawn as lines.
    """

    colour: str | None = None
    """
    Colour to draw a global-mean dataset's line in

    If `None`, the figure picks one.
    Ignored for spatial datasets, which are coloured by latitude.
    """

    linestyle: str = "-"
    """Line style to draw a global-mean dataset in"""

    alpha: float = 1.0
    """
    How opaque to draw a spatial dataset's points

    For datasets with so many points that, drawn solid,
    they hide what is drawn underneath them.
    """

    region: str = "Global"
    """
    Region a spatial mean dataset covers

    Only matters for datasets without spatial information.
    Used to draw e.g. CMIP6's hemispheric means in the same colour
    as our hemispheric means.
    """

    legend_group: str = "Comparison data"
    """
    Group to list the dataset under in the figures' legends

    The CMIP forcings are listed apart from the other comparison data,
    because they are what we are updating, not an independent check on it.
    """

    notes: tuple[str, ...] = field(default_factory=tuple)
    """Anything worth knowing about the dataset which the figure doesn't say"""

    @property
    def is_spatial(self) -> bool:
        """Whether the dataset's values each have a latitude"""
        return LATITUDE_COLUMN in self.data.columns

    def to_units(
        self, units: str, ur: pint.UnitRegistry = openscm_units.unit_registry
    ) -> ComparisonTimeseries:
        """
        Get the dataset in different units

        Parameters
        ----------
        units
            Units to convert to

        ur
            Unit registry to use for the conversion

        Returns
        -------
        :
            The dataset, in `units`
        """
        if units == self.units:
            return self

        data = self.data.copy()
        data[VALUE_COLUMN] = (
            ur.Quantity(data[VALUE_COLUMN].to_numpy(), self.units).to(units).m
        )

        return dataclasses.replace(self, data=data, units=units)


RADIATIVE_EFFICIENCIES: dict[str, pint.Quantity] = {
    # From the historical user guide
    # (Table 7.SM.7 of IPCC AR6 WG1 Ch. 7 Supplementary Material)
    "co2": Q(1.33e-5, "W / m^2 / ppb"),
    "ch4": Q(3.88e-4, "W / m^2 / ppb"),
    "n2o": Q(3.2e-3, "W / m^2 / ppb"),
    # The rest mirror `RADIATIVE_EFFICIENCIES` in the original run's
    # `3501_calculate-full-equivalence` notebook,
    # which cites Table 7.SM.6 of
    # https://www.ipcc.ch/report/ar6/wg1/downloads/report/IPCC_AR6_WGI_Chapter07_SM.pdf
    # Chlorofluorocarbons
    "cfc11": Q(0.291, "W / m^2 / ppb"),
    "cfc11eq": Q(0.291, "W / m^2 / ppb"),
    "cfc12": Q(0.358, "W / m^2 / ppb"),
    "cfc12eq": Q(0.358, "W / m^2 / ppb"),
    "cfc113": Q(0.301, "W / m^2 / ppb"),
    "cfc114": Q(0.314, "W / m^2 / ppb"),
    "cfc115": Q(0.246, "W / m^2 / ppb"),
    # Hydrofluorochlorocarbons
    "hcfc22": Q(0.214, "W / m^2 / ppb"),
    "hcfc141b": Q(0.161, "W / m^2 / ppb"),
    "hcfc142b": Q(0.193, "W / m^2 / ppb"),
    # Hydrofluorocarbons
    "hfc23": Q(0.191, "W / m^2 / ppb"),
    "hfc32": Q(0.111, "W / m^2 / ppb"),
    "hfc125": Q(0.234, "W / m^2 / ppb"),
    "hfc134a": Q(0.167, "W / m^2 / ppb"),
    "hfc134aeq": Q(0.167, "W / m^2 / ppb"),
    "hfc143a": Q(0.168, "W / m^2 / ppb"),
    "hfc152a": Q(0.102, "W / m^2 / ppb"),
    "hfc227ea": Q(0.273, "W / m^2 / ppb"),
    "hfc236fa": Q(0.251, "W / m^2 / ppb"),
    "hfc245fa": Q(0.245, "W / m^2 / ppb"),
    "hfc365mfc": Q(0.228, "W / m^2 / ppb"),
    "hfc4310mee": Q(0.357, "W / m^2 / ppb"),
    # Chlorocarbons and Hydrochlorocarbons
    "ch3ccl3": Q(0.065, "W / m^2 / ppb"),
    "ccl4": Q(0.166, "W / m^2 / ppb"),
    "ch3cl": Q(0.005, "W / m^2 / ppb"),
    "ch2cl2": Q(0.029, "W / m^2 / ppb"),
    "chcl3": Q(0.074, "W / m^2 / ppb"),
    # Bromocarbons, Hydrobromocarbons and Halons
    "ch3br": Q(0.004, "W / m^2 / ppb"),
    "halon1211": Q(0.300, "W / m^2 / ppb"),
    "halon1301": Q(0.299, "W / m^2 / ppb"),
    "halon2402": Q(0.312, "W / m^2 / ppb"),
    # Fully Fluorinated Species
    "nf3": Q(0.204, "W / m^2 / ppb"),
    "sf6": Q(0.567, "W / m^2 / ppb"),
    "so2f2": Q(0.211, "W / m^2 / ppb"),
    "cf4": Q(0.099, "W / m^2 / ppb"),
    "c2f6": Q(0.261, "W / m^2 / ppb"),
    "c3f8": Q(0.270, "W / m^2 / ppb"),
    "cc4f8": Q(0.314, "W / m^2 / ppb"),
    "c4f10": Q(0.369, "W / m^2 / ppb"),
    "c5f12": Q(0.408, "W / m^2 / ppb"),
    "c6f14": Q(0.449, "W / m^2 / ppb"),
    "c7f16": Q(0.503, "W / m^2 / ppb"),
    "c8f18": Q(0.558, "W / m^2 / ppb"),
}
"""Radiative efficiency of each gas

Used to put an approximate radiative effect next to concentrations,
i.e. the concentration multiplied by the radiative efficiency.
This is not effective radiative forcing,
but it puts differences between gases on the same footing.
"""


def get_radiative_effect_per_unit(
    gas: str, units: str, ur: pint.UnitRegistry = openscm_units.unit_registry
) -> float:
    """
    Get the radiative effect of one unit of a gas' concentration

    Parameters
    ----------
    gas
        Gas of interest

    units
        Units the gas' concentration is in

    ur
        Unit registry to use for the conversion

    Returns
    -------
    :
        Radiative effect of one `units` of `gas`, in W / m^2
    """
    return (ur.Quantity(1.0, units) * RADIATIVE_EFFICIENCIES[gas]).to("W / m^2").m


CMIP6_SOURCE_ID = "UoM-CMIP-1-2-0"
"""Source ID of the CMIP6 historical concentrations"""

ESGF_LOCAL_DATA_ROOT_DIR = DATA_RAW_DIR / "esgf"
"""Where data downloaded from ESGF is kept

The same place the user guide keeps it, so they share a cache.
"""

ESGF_DATABASE_FILE = REPO_ROOT / "download-test-database.db"
"""Database of what we have already found on and downloaded from ESGF

The same one the user guide uses, so they share a cache.
"""

CMIP6_SECTORS = {
    0: "Global",
    1: "Northern hemisphere",
    2: "Southern hemisphere",
}
"""Region each of the CMIP6 global-mean files' sectors stands for

Named as our own spatial means are named in the figures,
so they can be drawn in the same colours.
"""


@cache
def get_esgf_engine(
    database_file: Path = ESGF_DATABASE_FILE,
) -> sqlalchemy.engine.base.Engine:
    """
    Get the engine for the database of what we have found on ESGF

    Parameters
    ----------
    database_file
        Database file

    Returns
    -------
    :
        Engine for the database
    """
    engine = get_sqlite_engine(database_file)
    create_all_tables(engine)

    return engine


def fix_cmip6_conventions_keep_sectors(
    ds: xr.Dataset, esgf_dataset: ESGFDataset, file_paths: list[Path]
) -> xr.Dataset:
    """
    Fix the CMIP6 data's conventions, keeping its hemispheric means

    [local.data_loading.fix_conventions_to_match_cmip7][] throws away
    the hemispheric means, because CMIP7 puts them in a separate file.
    We want to compare against them, so we keep them.

    Parameters
    ----------
    ds
        Dataset to fix

    esgf_dataset
        [ESGFDataset][] from which `ds` was loaded

    file_paths
        File paths from which `ds` was loaded

    Returns
    -------
    :
        `ds`, with its units fixed and its sectors labelled by region
    """
    ds = ds.drop_vars("sector_bnds", errors="ignore")
    ds = ds.assign_coords(sector=[CMIP6_SECTORS[int(v)] for v in ds["sector"].values])

    unit_map = {
        "1.e-6": "ppm",
        "1.e-9": "ppb",
        "1.e-12": "ppt",
    }
    variable = ds.attrs["variable_id"]
    ds[variable].attrs["units"] = unit_map[ds[variable].attrs["units"]]

    return ds


def load_cmip6_ghg_ds_keep_sectors(esgf_dataset: ESGFDataset) -> xr.Dataset:
    """
    Load a CMIP6 greenhouse gas dataset, keeping its hemispheric means

    Parameters
    ----------
    esgf_dataset
        [ESGFDataset][] from which to load

    Returns
    -------
    :
        Loaded data
    """
    return load_xarray_from_esgf_dataset(
        esgf_dataset=esgf_dataset,
        pre_to_xarray=fix_broken_calendar_spec,
        post_to_xarray=fix_cmip6_conventions_keep_sectors,
        add_attributes_from_metadata=("cmip_era", "source_id"),
    )


def get_cmip6_spatial_means(gas: str, time_sampling: str) -> xr.DataArray:
    """
    Get the CMIP6 global- and hemispheric-means of a gas

    Downloaded from ESGF the first time, re-used from disk after that.

    Parameters
    ----------
    gas
        Gas of interest

    time_sampling
        Time sampling to get, "yr" or "mon"

    Returns
    -------
    :
        The CMIP6 global- and hemispheric-means, with a `sector` dimension
        (labelled with the region names used in the figures)
        and a `time` dimension (as decimal years)
    """
    ds = fetch_and_load_ghg_dataset(
        local_data_root_dir=ESGF_LOCAL_DATA_ROOT_DIR,
        index_node=KnownIndexNode.ORNL,
        ghg=gas,
        grid="gm",
        time_sampling=time_sampling,
        cmip_era="CMIP6",
        source_id=CMIP6_SOURCE_ID,
        engine=get_esgf_engine(),
        load_xr_dataset=load_cmip6_ghg_ds_keep_sectors,
    )

    da = ds[gas].compute()
    # Mid-point of each year or month, as a decimal year,
    # the same way our own output is placed on the time axis
    if time_sampling == "yr":
        time = np.array([t.year + 0.5 for t in da["time"].values])
    elif time_sampling == "mon":
        time = np.array([t.year + (t.month - 0.5) / 12.0 for t in da["time"].values])
    else:
        raise NotImplementedError(time_sampling)

    return da.assign_coords(time=time)


def get_cmip6_comparisons(
    gas: str,
    time_sampling: str,
    regions: tuple[str, ...] = ("Global",),
    label: str = "CMIP6",
) -> tuple[ComparisonTimeseries, ...]:
    """
    Get the CMIP6 concentrations as comparison datasets

    Parameters
    ----------
    gas
        Gas of interest

    time_sampling
        Time sampling to get, "yr" or "mon"

    regions
        Regions to get

    label
        Label to give the dataset

        The region is added to it.

    Returns
    -------
    :
        One comparison dataset per region
    """
    da = get_cmip6_spatial_means(gas, time_sampling)

    res = []
    for region in regions:
        region_da = da.sel(sector=region)
        res.append(
            ComparisonTimeseries(
                label=(
                    f"{label} {'global-mean' if region == 'Global' else region.lower()}"
                ),
                data=pd.DataFrame(
                    {
                        TIME_COLUMN: region_da["time"].values,
                        VALUE_COLUMN: region_da.values,
                    }
                ),
                units=da.attrs["units"],
                linestyle="--",
                region=region,
                legend_group="CMIP forcings",
            )
        )

    return tuple(res)


CH4_ICE_CORE_SUPPLEMENT_URL = (
    "https://media.springernature.com/original/springer-static/esm/"
    "art%3A10.1038%2Fs41586-026-10938-1/MediaObjects/"
    "41586_2026_10938_MOESM2_ESM.xlsx"
)
"""Supplementary data of https://www.nature.com/articles/s41586-026-10938-1"""

CH4_ICE_CORE_SUPPLEMENT_FILE = (
    DATA_RAW_DIR / "comparison-data" / "41586_2026_10938_MOESM2_ESM.xlsx"
)
"""Where we keep our copy of the CH4 ice core supplementary data

Ignored by git, because we can always download it again.
"""

# TODO: check
CH4_ICE_CORE_SUPPLEMENT_SITES = {
    # Column prefix: (label, latitude, marker)
    # The new record the paper presents
    "Summit": ("Summit", 72.58, "D"),
    "WAIS": ("WAIS Divide", -79.47, "v"),
    "GISP2": ("GISP2", 72.60, "^"),
    # Same site and latitude as the Law Dome data we build our record from
    "LawDome": ("Law Dome", -66.73, "s"),
    # Same site and latitude as the NEEM data we optimise our record against
    "NEEM": ("NEEM", 77.45, "P"),
    "ML": ("Mauna Loa", 19.54, "X"),
}
"""The sites in the CH4 ice core supplementary data

The spreadsheet has no locations in it, so they are written out here.
Each site gets its own marker, because colour is taken by latitude.
"""

CH4_SUPPLEMENT_ATMOSPHERIC_SITES = ("ML",)
"""The sites in the CH4 ice core supplementary data which are not ice cores

These are a handful of recent atmospheric measurements,
so they are drawn solid rather than see-through like the ice cores,
see [ICE_CORE_ALPHA][].
"""


def ensure_file_downloaded(
    url: str, out_file: Path, headers: dict[str, str] | None = None
) -> Path:
    """
    Download a file, unless we already have it

    Parameters
    ----------
    url
        URL to download from

    out_file
        Where to save the file

    headers
        Headers to send with the request

    Returns
    -------
    :
        `out_file`
    """
    if out_file.exists():
        logger.info(f"Using existing {out_file}")
        return out_file

    out_file.parent.mkdir(exist_ok=True, parents=True)
    logger.info(f"Downloading {url} to {out_file}")
    request = urllib.request.Request(url, headers=headers or {})  # noqa: S310
    with urllib.request.urlopen(request) as response:  # noqa: S310
        out_file.write_bytes(response.read())

    return out_file


ICE_CORE_ALPHA = 0.3
"""How opaque to draw the CH4 ice core records' points

There are hundreds of them before the industrial era,
and drawn solid they cover our output (and CMIP6's) completely.
"""


def get_ch4_ice_core_comparisons(
    sites: tuple[str, ...] | None = None,
) -> tuple[ComparisonTimeseries, ...]:
    """
    Get the CH4 records from the supplementary data of the Summit ice core paper

    https://www.nature.com/articles/s41586-026-10938-1

    Parameters
    ----------
    sites
        Sites to get, by their column prefix in the spreadsheet

        If `None`, every site in [CH4_ICE_CORE_SUPPLEMENT_SITES][].

    Returns
    -------
    :
        One comparison dataset per site
    """
    if sites is None:
        sites = tuple(CH4_ICE_CORE_SUPPLEMENT_SITES)

    raw = pd.read_excel(
        ensure_file_downloaded(
            CH4_ICE_CORE_SUPPLEMENT_URL, CH4_ICE_CORE_SUPPLEMENT_FILE
        ),
        sheet_name="Figure 1ab",
    )

    res = []
    for site in sites:
        label, latitude, marker = CH4_ICE_CORE_SUPPLEMENT_SITES[site]
        site_df = (
            raw[[f"{site}_Year", f"{site}_CH4"]]
            .dropna()
            .rename(columns={f"{site}_Year": TIME_COLUMN, f"{site}_CH4": VALUE_COLUMN})
            .sort_values(TIME_COLUMN)
        )
        site_df[LATITUDE_COLUMN] = latitude
        res.append(
            ComparisonTimeseries(
                label=label,
                data=site_df.reset_index(drop=True),
                units="ppb",
                marker=marker,
                alpha=(
                    1.0 if site in CH4_SUPPLEMENT_ATMOSPHERIC_SITES else ICE_CORE_ALPHA
                ),
            )
        )

    return tuple(res)


NOAA_TRENDS_URL = "https://gml.noaa.gov/webdata/ccgg/trends"
"""Where NOAA GML's trends data lives

See https://gml.noaa.gov/ccgg/trends/
(CO2, including https://gml.noaa.gov/ccgg/trends/mlo.html)
and https://gml.noaa.gov/ccgg/trends_ch4/ (CH4, N2O and SF6).
"""

NOAA_TRENDS_DIR = DATA_RAW_DIR / "comparison-data" / "noaa-trends"
"""Where we keep our copy of NOAA GML's trends data

Tracked in git, because NOAA revise these files every month,
so downloading them again won't necessarily give the same numbers.
Delete them to pick up the latest data.
"""

NOAA_TRENDS_UNITS = {
    "co2": "ppm",
    "ch4": "ppb",
    "n2o": "ppb",
    "sf6": "ppt",
}
"""Units of each gas NOAA GML report trends for, as their files state them"""

MAUNA_LOA_LATITUDE = 19.54
"""Latitude of the Mauna Loa Observatory"""


def load_noaa_trends_file(gas: str, filename: str, value_column: str) -> pd.DataFrame:
    """
    Load one of NOAA GML's monthly trends files

    Parameters
    ----------
    gas
        Gas the file is for

    filename
        Name of the file, e.g. `co2_mm_gl.csv`

    value_column
        Column to take the values from

    Returns
    -------
    :
        Data with a [TIME_COLUMN][] and a [VALUE_COLUMN][]
    """
    raw = pd.read_csv(
        ensure_file_downloaded(
            f"{NOAA_TRENDS_URL}/{gas}/{filename}", NOAA_TRENDS_DIR / filename
        ),
        comment="#",
    )
    # Mid-point of each month, as a decimal year,
    # the same way our own output is placed on the time axis
    res = pd.DataFrame(
        {
            TIME_COLUMN: raw["year"] + (raw["month"] - 0.5) / 12.0,
            VALUE_COLUMN: raw[value_column],
        }
    )
    # Missing values are marked with negative numbers
    return res[res[VALUE_COLUMN] > 0.0].reset_index(drop=True)


def get_noaa_comparisons(
    gas: str, deseasonalised: bool
) -> tuple[ComparisonTimeseries, ...]:
    """
    Get NOAA GML's records of a gas

    For every gas, NOAA's global-mean.
    For CO2, also the Mauna Loa record.

    Parameters
    ----------
    gas
        Gas of interest, one of [NOAA_TRENDS_UNITS][]

    deseasonalised
        Whether to get the records with their seasonal cycle removed

        If `False`, the monthly means.

    Returns
    -------
    :
        NOAA's records of `gas` (for CO2, the Mauna Loa record first)
    """
    units = NOAA_TRENDS_UNITS[gas]
    notes = ("Seasonal cycle removed by NOAA",) if deseasonalised else ()

    res = []
    if gas == "co2":
        mauna_loa = load_noaa_trends_file(
            gas,
            "co2_mm_mlo.csv",
            "deseasonalized" if deseasonalised else "average",
        )
        mauna_loa[LATITUDE_COLUMN] = MAUNA_LOA_LATITUDE
        res.append(
            ComparisonTimeseries(
                label="NOAA Mauna Loa",
                data=mauna_loa,
                units=units,
                spatial_as_line=True,
                notes=notes,
            )
        )

    res.append(
        ComparisonTimeseries(
            label="NOAA global-mean",
            data=load_noaa_trends_file(
                gas,
                f"{gas}_mm_gl.csv",
                # NOAA call the global-mean with its seasonal cycle removed the trend
                "trend" if deseasonalised else "average",
            ),
            units=units,
            colour=OKABE_ITO["bluish green"],
            notes=notes,
        )
    )

    return tuple(res)


GCP_CH4_2024_OBJECT_ID = "-MqGCn38zUlEi4_aQm-1_w2I"
"""ICOS Carbon Portal ID of the Global Methane Budget 2000-2020 data supplement

Saunois et al. (2025), https://doi.org/10.5194/essd-17-1873-2025,
data at https://doi.org/10.18160/GKQ9-2RHT.
"""

GCP_CH4_2024_FILE = (
    DATA_RAW_DIR
    / "comparison-data"
    / "gcp-ch4-2024"
    / "Global_methane_budget_2023_2000_2020_v1.xlsx"
)
"""Where we keep our copy of the Global Methane Budget 2000-2020 data supplement

Tracked in git, because the Global Carbon Project may revise it in place,
so downloading it again won't necessarily give the same numbers.
"""


def get_uci_ch4_comparison(deseasonalised: bool) -> ComparisonTimeseries:
    """
    Get UCI's global-mean CH4 record

    UCI (University of California, Irvine) sample the remote Pacific
    (71N to 46S) every three months,
    so the record is quarterly.
    The record's reference is Simpson et al. (2012),
    https://doi.org/10.1038/nature11342.
    The archived version of it (https://doi.org/10.3334/CDIAC/ATG.002)
    stops in 2009, so we take the version compiled for
    the Global Methane Budget 2000-2020, which runs to the end of 2022.

    Parameters
    ----------
    deseasonalised
        Whether to get the record with its seasonal cycle removed

        If `False`, the quarterly means.

    Returns
    -------
    :
        UCI's global-mean CH4 record
    """
    raw = pd.read_excel(
        ensure_file_downloaded(
            f"https://data.icos-cp.eu/objects/{GCP_CH4_2024_OBJECT_ID}",
            GCP_CH4_2024_FILE,
            # The portal only serves the file once its licence has been accepted,
            # which the browser records in this cookie
            headers={"Cookie": f"CpLicenseAcceptedFor={GCP_CH4_2024_OBJECT_ID}"},
        ),
        sheet_name="CH4_observation -fig 1",
        header=None,
    )
    # The sheet holds one block of five columns per network, side by side,
    # (date, mixing ratio, deseasonalised mixing ratio, trend, growth rate).
    # UCI's is the fourth block,
    # even though its header says "CSIRO" (CSIRO's is the third block).
    # The sheet's notes list the networks as NOAA, AGAGE, CSIRO, UCI,
    # and this block's quarterly sampling and 1978 start are UCI's.
    uci_first_column = 15
    header_row = 13
    time_column = uci_first_column
    value_column = uci_first_column + (2 if deseasonalised else 1)
    if not str(raw.iloc[header_row, value_column]).startswith(
        "Deseasonalized" if deseasonalised else "Mixing ratio"
    ):
        raise AssertionError(raw.iloc[header_row, value_column])

    data = (
        raw.iloc[header_row + 1 :, [time_column, value_column]].dropna().astype(float)
    )
    data.columns = [TIME_COLUMN, VALUE_COLUMN]

    return ComparisonTimeseries(
        label="UCI global-mean",
        data=data.sort_values(TIME_COLUMN).reset_index(drop=True),
        units="ppb",
        colour=OKABE_ITO["reddish purple"],
        notes=("Seasonal cycle removed by the GCP",) if deseasonalised else (),
    )


# TODO: check IGCC processing
IGCC_RELEASE = "v6.4.0"
"""Release of the Indicators of Global Climate Change (IGCC) forcing timeseries

The snapshot used in IGCC 2025 (Forster et al., 2026,
https://doi.org/10.5194/essd-18-3889-2026), archived at
https://doi.org/10.5281/zenodo.20498594.
Later releases (up to v6.4.2 at the time of writing)
add forcing categories but leave the concentrations unchanged.
"""

IGCC_CONCENTRATIONS_URL = (
    "https://raw.githubusercontent.com/ClimateIndicator/forcing-timeseries/"
    f"{IGCC_RELEASE}/output/ghg_concentrations.csv"
)
"""Where IGCC's global-, annual-mean concentrations live"""

IGCC_CONCENTRATIONS_FILE = (
    DATA_RAW_DIR / "comparison-data" / "igcc" / IGCC_RELEASE / "ghg_concentrations.csv"
)
"""Where we keep our copy of IGCC's concentrations

Ignored by git: the release is pinned, so downloading it again
gives the same numbers.
"""

IGCC_COLUMNS = {
    "co2": "CO2",
    "ch4": "CH4",
    "n2o": "N2O",
    "c2f6": "C2F6",
    "c3f8": "C3F8",
    "ccl4": "CCl4",
    "cf4": "CF4",
    "cfc11": "CFC-11",
    "cfc113": "CFC-113",
    "cfc114": "CFC-114",
    "cfc115": "CFC-115",
    "cfc12": "CFC-12",
    "ch2cl2": "CH2Cl2",
    "ch3br": "CH3Br",
    "ch3ccl3": "CH3CCl3",
    "ch3cl": "CH3Cl",
    "chcl3": "CHCl3",
    "halon1211": "Halon-1211",
    "halon1301": "Halon-1301",
    "halon2402": "Halon-2402",
    "hcfc141b": "HCFC-141b",
    "hcfc142b": "HCFC-142b",
    "hcfc22": "HCFC-22",
    "hfc125": "HFC-125",
    "hfc134a": "HFC-134a",
    "hfc143a": "HFC-143a",
    "hfc152a": "HFC-152a",
    "hfc227ea": "HFC-227ea",
    "hfc23": "HFC-23",
    "hfc236fa": "HFC-236fa",
    "hfc245fa": "HFC-245fa",
    "hfc32": "HFC-32",
    "hfc365mfc": "HFC-365mfc",
    "hfc4310mee": "HFC-43-10mee",
    "nf3": "NF3",
    "sf6": "SF6",
    "so2f2": "SO2F2",
    "cc4f8": "c-C4F8",
    "c4f10": "n-C4F10",
    "c5f12": "n-C5F12",
    # IGCC split C6F14 into its isomers, we only have the one
    "c6f14": "n-C6F14",
    "c7f16": "C7F16",
    "c8f18": "C8F18",
    # These group different gases to our equivalent species, with different
    # radiative efficiencies, so they are not like-for-like comparisons.
    # They are shown anyway, to make clear that the two are not comparable.
    "cfc12eq": "CFC[CFC-12-eq]",
    "hfc134aeq": "HFC[HFC-134a-eq]",
}
"""Column of IGCC's concentrations file which holds each of our gases

For the equivalent species, the column holds IGCC's equivalent,
which is not the same thing as ours.
IGCC's CFC-12 equivalent includes eight more gases than ours
(CFC-13, CFC-112, CFC-112a, CFC-113a, CFC-114a, HCFC-133a, HCFC-31, HCFC-124),
its HFC-134a equivalent includes only the HFCs
(i.e. none of the PFCs, SF6, NF3 or SO2F2),
and both are calculated with the radiative efficiencies of Hodnebrog et al. (2020)
rather than those of AR6.
See `notebooks/01_trace-gas-global-mean.py` in IGCC's repository.
IGCC has no CFC-11 equivalent.
"""


def get_igcc_units(gas: str) -> str:
    """
    Get the units of a gas in IGCC's concentrations file

    The file doesn't say, so this follows IGCC's convention.

    Parameters
    ----------
    gas
        Gas of interest

    Returns
    -------
    :
        Units of `gas` in IGCC's concentrations file
    """
    if gas == "co2":
        return "ppm"

    if gas in ("ch4", "n2o"):
        return "ppb"

    return "ppt"


def get_igcc_comparison(gas: str) -> ComparisonTimeseries:
    """
    Get IGCC's global-, annual-mean record of a gas

    Forster et al. (2026), Indicators of Global Climate Change 2025,
    https://doi.org/10.5194/essd-18-3889-2026.
    The record is compiled from NOAA and AGAGE data
    (and, before those, the AR6 concentrations,
    which are themselves partly based on the CMIP6 concentrations),
    so it is not independent of our output or of CMIP6.

    Parameters
    ----------
    gas
        Gas of interest, one of [IGCC_COLUMNS][]

        For an equivalent species, this is IGCC's equivalent,
        which is defined differently to ours (see [IGCC_COLUMNS][]).

    Returns
    -------
    :
        IGCC's record of `gas`
    """
    raw = pd.read_csv(
        ensure_file_downloaded(IGCC_CONCENTRATIONS_URL, IGCC_CONCENTRATIONS_FILE),
        index_col="YYYY",
    )
    # The file holds 1750, then jumps to 1850.
    # 1750 is the AR6 pre-industrial reference value rather than part of the record,
    # and drawn as part of the line it would be joined to 1850 by a straight line
    # which isn't in the data.
    first_year = 1850
    values = raw.loc[raw.index >= first_year, IGCC_COLUMNS[gas]]

    return ComparisonTimeseries(
        label="IGCC global-mean",
        data=pd.DataFrame(
            {
                # Annual-means, so placed at the middle of the year,
                # the same way our own output is placed on the time axis
                TIME_COLUMN: values.index.to_numpy(dtype=float) + 0.5,
                VALUE_COLUMN: values.to_numpy(dtype=float),
            }
        ),
        units=get_igcc_units(gas),
        colour=OKABE_ITO["orange"],
        linestyle="-.",
        notes=("1750 value (AR6 pre-industrial reference) not shown",),
    )
