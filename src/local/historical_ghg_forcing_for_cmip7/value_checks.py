"""
Checks of the values behind the statements in the historical manuscript

Each check calculates the value behind a statement,
and is matched by its tag to the `% value-check: {...}` comments in the latex.
See [local.value_checks][] for how that works.

Differences are always our (CMIP7) dataset minus the other dataset.
"""

from __future__ import annotations

import json
from collections.abc import Callable
from functools import cache, partial
from pathlib import Path
from typing import Any, Literal

import numpy as np
import openscm_units
import pandas as pd
import pint
import xarray as xr
import yaml

from local.cmip_ghg_generation import BUNDLE_CONFIG_FILE, DEFAULT_BUNDLE_DIR
from local.historical_ghg_forcing_for_cmip7.c4f10_like_methods_figure import (
    C4F10_LIKE_GASES,
)
from local.historical_ghg_forcing_for_cmip7.cfc12_like_methods_figure import (
    CFC12_LIKE_GASES,
    get_cfc12_like_all_data_with_bins,
    get_global_mean_supplement,
    get_step_config,
    interim_dir,
    supplement_replaces_obs_network,
)
from local.historical_ghg_forcing_for_cmip7.comparison_data import (
    MAUNA_LOA_LATITUDE,
    RADIATIVE_EFFICIENCIES,
    TIME_COLUMN,
    VALUE_COLUMN,
    ComparisonTimeseries,
    get_cmip6_comparisons,
    get_igcc_comparison,
    get_noaa_comparisons,
    get_uci_ch4_comparison,
)
from local.historical_ghg_forcing_for_cmip7.equivalent_species import (
    EquivalenceDataset,
    check_reproduction,
    decompose_difference,
    load_cmip6,
    load_cmip7,
    load_igcc,
)
from local.historical_ghg_forcing_for_cmip7.results_figure import load_output
from local.historical_ghg_forcing_for_cmip7.zenodo_missing import (
    CH4_ICE_CORE_FILES,
    get_ice_core_latitude,
    load_ch4_ice_core,
)
from local.paths import DATA_RAW_DIR
from local.value_checks import CheckValue, ValueCheck

Q = openscm_units.unit_registry.Quantity

LAST_N_YEARS = 10
"""Number of years the statements about 'the last ten years' cover"""

CMIP6_LAST_YEAR = 2014
"""
Last year of CMIP6's historical dataset

Checked against the CMIP6 data by the `{gas}-cmip6-last-year` value checks.
"""

ENERGY_BALANCE_THRESHOLD = Q(0.01, "W / m^2")
"""
Change in energy balance below which a choice is not critical

As stated in the output requirements section of the manuscript.
"""

NORTHERN_HEMISPHERE = "Northern hemisphere"
"""Name of the northern hemisphere in our hemispheric-mean output"""

SOUTHERN_HEMISPHERE = "Southern hemisphere"
"""Name of the southern hemisphere in our hemispheric-mean output"""


@cache
def get_output(
    gas: str, *, bundle_dir: Path
) -> tuple[xr.DataArray, xr.DataArray, xr.DataArray, xr.DataArray]:
    """
    Get our output for a gas

    Parameters
    ----------
    gas
        Gas of interest

    bundle_dir
        Directory in which to keep the original run's bundle

    Returns
    -------
    :
        Our native resolution output, yearly global-mean,
        monthly global-mean and monthly hemispheric-means
        (see [local.historical_ghg_forcing_for_cmip7.results_figure.load_output][])
    """
    return load_output(gas, bundle_dir)


def get_units(gas: str, *, bundle_dir: Path) -> str:
    """
    Get the units of our output for a gas

    Parameters
    ----------
    gas
        Gas of interest

    bundle_dir
        Directory in which to keep the original run's bundle

    Returns
    -------
    :
        Units of our output for `gas`
    """
    return str(get_output(gas, bundle_dir=bundle_dir)[1].attrs["units"])


def get_last_n_years(gas: str, *, bundle_dir: Path) -> np.ndarray:
    """
    Get the last [LAST_N_YEARS][] years of our output for a gas

    Parameters
    ----------
    gas
        Gas of interest

    bundle_dir
        Directory in which to keep the original run's bundle

    Returns
    -------
    :
        The years
    """
    last_year = int(get_output(gas, bundle_dir=bundle_dir)[1]["year"].max())

    return np.arange(last_year - LAST_N_YEARS + 1, last_year + 1)


def get_last_year(gas: str, *, bundle_dir: Path) -> pint.Quantity:
    """
    Get the last year of our output for a gas

    Parameters
    ----------
    gas
        Gas of interest

    bundle_dir
        Directory in which to keep the original run's bundle

    Returns
    -------
    :
        Last year of our global-, annual-mean
    """
    return Q(int(get_output(gas, bundle_dir=bundle_dir)[1]["year"].max()), "yr")


def get_written_year(
    key: Literal["start_year", "end_year"], *, bundle_dir: Path
) -> pint.Quantity:
    """
    Get the first or last year for which the original run wrote its output

    Parameters
    ----------
    key
        Which year to get

    bundle_dir
        Directory in which to keep the original run's bundle

    Returns
    -------
    :
        The year, which is the same for every gas

    Raises
    ------
    AssertionError
        The year is not the same for every gas
    """
    config = yaml.safe_load((bundle_dir / BUNDLE_CONFIG_FILE).read_text())
    years = {step_config[key] for step_config in config["write_input4mips"]}
    if len(years) != 1:
        msg = f"{key} is not the same for every gas: {sorted(years)}"
        raise AssertionError(msg)

    return Q(int(years.pop()), "yr")


def get_global_annual_mean(
    gas: str, year: int | None = None, *, bundle_dir: Path
) -> pint.Quantity:
    """
    Get our global-, annual-mean for a gas in a given year

    Parameters
    ----------
    gas
        Gas of interest

    year
        Year of interest

        If `None`, the last year of our output.

    bundle_dir
        Directory in which to keep the original run's bundle

    Returns
    -------
    :
        Our global-, annual-mean in `year`
    """
    ours = get_output(gas, bundle_dir=bundle_dir)[1].to_series()
    if year is None:
        year = int(ours.index.max())

    return Q(float(ours.loc[year]), get_units(gas, bundle_dir=bundle_dir))


def get_decimal_month(
    year: pd.Index | np.ndarray, month: pd.Index | np.ndarray
) -> np.ndarray:
    """
    Get the middle of each month as a decimal year

    The same way our output and the comparison data are placed on the time axis.
    Rounded, so that the same month from different sources lines up exactly.

    Parameters
    ----------
    year
        Year of each month

    month
        Month (1 to 12)

    Returns
    -------
    :
        Middle of each month, as a decimal year
    """
    return np.round(np.asarray(year) + (np.asarray(month) - 0.5) / 12.0, 4)


def monthly_to_series(da: xr.DataArray) -> pd.Series[float]:
    """
    Get data with a year and month dimension as a series indexed by decimal month

    Parameters
    ----------
    da
        Data with (only) `year` and `month` dimensions

    Returns
    -------
    :
        `da`, indexed by the middle of each month as a decimal year
    """
    res = da.to_series()
    res.index = get_decimal_month(
        res.index.get_level_values("year"), res.index.get_level_values("month")
    )

    return res


def comparison_to_series(
    comparison: ComparisonTimeseries, units: str, yearly: bool
) -> pd.Series[float]:
    """
    Get a (global-mean) comparison dataset as a series

    Parameters
    ----------
    comparison
        Comparison dataset

    units
        Units to convert to

    yearly
        Whether the data is annual-means

        If `True`, the result is indexed by (integer) year.
        Otherwise, by the middle of each month as a decimal year.

    Returns
    -------
    :
        The comparison dataset's values
    """
    data = comparison.to_units(units).data
    if yearly:
        index = np.floor(data[TIME_COLUMN].to_numpy()).astype(int)
    else:
        index = np.round(data[TIME_COLUMN].to_numpy(), 4)

    return pd.Series(data[VALUE_COLUMN].to_numpy(), index=index)


def get_radiative_effect(gas: str, value: pint.Quantity) -> pint.Quantity:
    """
    Get the approximate radiative effect of (a difference in) concentration

    As described in the results, this is the concentration
    multiplied by the radiative efficiency.

    Parameters
    ----------
    gas
        Gas of interest

    value
        Concentration (difference)

    Returns
    -------
    :
        Approximate radiative effect of the magnitude of `value`
    """
    # Change this to not being abs
    return (abs(value) * RADIATIVE_EFFICIENCIES[gas]).to("W / m^2")


def max_abs(diff: pd.Series[float], units: str) -> pint.Quantity:
    """
    Get the maximum absolute value of a difference

    Parameters
    ----------
    diff
        Difference

    units
        Units of `diff`

    Returns
    -------
    :
        Maximum absolute value of `diff`
    """
    return Q(float(diff.abs().max()), units)


def get_last_n_years_trend(gas: str, *, bundle_dir: Path) -> pint.Quantity:
    """
    Get the linear trend in our global-, annual-mean over the last ten years

    Parameters
    ----------
    gas
        Gas of interest

    bundle_dir
        Directory in which to keep the original run's bundle

    Returns
    -------
    :
        Trend
    """
    years = get_last_n_years(gas, bundle_dir=bundle_dir)
    values = get_output(gas, bundle_dir=bundle_dir)[1].sel(year=years).values
    gradient = np.polyfit(years, values, 1)[0]

    return Q(gradient, f"{get_units(gas, bundle_dir=bundle_dir)} / yr")


def get_last_n_years_second_derivative(gas: str, *, bundle_dir: Path) -> pint.Quantity:
    """
    Get the second derivative of our global-, annual-mean over the last ten years

    From a quadratic fit.
    Negative means the trend is getting more negative
    (e.g. a decline is accelerating).

    Parameters
    ----------
    gas
        Gas of interest

    bundle_dir
        Directory in which to keep the original run's bundle

    Returns
    -------
    :
        Second derivative
    """
    years = get_last_n_years(gas, bundle_dir=bundle_dir)
    values = get_output(gas, bundle_dir=bundle_dir)[1].sel(year=years).values
    quadratic_coefficient = np.polyfit(years, values, 2)[0]

    return Q(
        2.0 * quadratic_coefficient, f"{get_units(gas, bundle_dir=bundle_dir)} / yr^2"
    )


def get_last_n_years_seasonal_cycle(
    gas: str, regions: tuple[str, ...], *, bundle_dir: Path
) -> pint.Quantity:
    """
    Get the size of the seasonal cycle over the last ten years

    The size is the peak-to-trough amplitude of the monthly hemispheric-mean
    in each year, averaged over the years.

    Parameters
    ----------
    gas
        Gas of interest

    regions
        Hemispheres to get it for

    bundle_dir
        Directory in which to keep the original run's bundle

    Returns
    -------
    :
        Size of the seasonal cycle in each of `regions`
    """
    hm_monthly = get_output(gas, bundle_dir=bundle_dir)[3].sel(
        year=get_last_n_years(gas, bundle_dir=bundle_dir), lat=list(regions)
    )
    amplitude = (hm_monthly.max("month") - hm_monthly.min("month")).mean("year")

    return Q(amplitude.values, get_units(gas, bundle_dir=bundle_dir))


def get_last_n_years_lat_gradient(gas: str, *, bundle_dir: Path) -> pint.Quantity:
    """
    Get the latitudinal gradient of our output over the last ten years

    The gradient of a linear fit to our native resolution output,
    averaged over the last ten years, against latitude.
    Positive means concentrations are higher in the north.

    Parameters
    ----------
    gas
        Gas of interest

    bundle_dir
        Directory in which to keep the original run's bundle

    Returns
    -------
    :
        Latitudinal gradient
    """
    mean = (
        get_output(gas, bundle_dir=bundle_dir)[0]
        .sel(year=get_last_n_years(gas, bundle_dir=bundle_dir))
        .mean(["year", "month"])
    )
    gradient = np.polyfit(mean["lat"].values, mean.values, 1)[0]

    return Q(gradient, f"{get_units(gas, bundle_dir=bundle_dir)} / degree")


@cache
def get_diff_from_cmip6(gas: str, *, bundle_dir: Path) -> pd.Series[float]:
    """
    Get the difference between our global-, annual-mean and CMIP6's

    Parameters
    ----------
    gas
        Gas of interest

    bundle_dir
        Directory in which to keep the original run's bundle

    Returns
    -------
    :
        Difference, by year, in our units
    """
    (cmip6,) = get_cmip6_comparisons(gas, "yr")
    ours = get_output(gas, bundle_dir=bundle_dir)[1].to_series()

    return (
        ours
        - comparison_to_series(
            cmip6, get_units(gas, bundle_dir=bundle_dir), yearly=True
        )
    ).dropna()


def get_max_abs_diff_from_cmip6(
    gas: str, start: int | None = None, *, bundle_dir: Path
) -> pint.Quantity:
    """
    Get the maximum absolute difference from CMIP6

    Parameters
    ----------
    gas
        Gas of interest

    start
        First year to consider

        If `None`, all years.

    bundle_dir
        Directory in which to keep the original run's bundle

    Returns
    -------
    :
        Maximum absolute difference from CMIP6
    """
    return max_abs(
        get_diff_from_cmip6(gas, bundle_dir=bundle_dir).loc[start:],
        get_units(gas, bundle_dir=bundle_dir),
    )


def get_year_of_max_abs_diff_from_cmip6(
    gas: str, start: int | None = None, *, bundle_dir: Path
) -> pint.Quantity:
    """
    Get the year in which the absolute difference from CMIP6 is greatest

    Parameters
    ----------
    gas
        Gas of interest

    start
        First year to consider

        If `None`, all years.

    bundle_dir
        Directory in which to keep the original run's bundle

    Returns
    -------
    :
        Year of maximum absolute difference
    """
    return Q(
        int(get_diff_from_cmip6(gas, bundle_dir=bundle_dir).loc[start:].abs().idxmax()),
        "yr",
    )


def get_cmip6_feature_size(
    gas: str, start: int, end: int, kind: Literal["dip", "spike"], *, bundle_dir: Path
) -> pint.Quantity:
    """
    Get the size of a dip or spike in CMIP6, relative to our dataset

    Parameters
    ----------
    gas
        Gas of interest

    start
        First year to look in

    end
        Last year to look in

    kind
        Whether CMIP6 dips below ours ("dip")
        or spikes above ours ("spike")

    bundle_dir
        Directory in which to keep the original run's bundle

    Returns
    -------
    :
        How far CMIP6 dips below (or spikes above) our dataset
        between `start` and `end`.
        Negative if CMIP6 doesn't dip below (or spike above) our dataset at all.
    """
    diff = get_diff_from_cmip6(gas, bundle_dir=bundle_dir).loc[start:end]
    size = diff.max() if kind == "dip" else -diff.min()

    return Q(float(size), get_units(gas, bundle_dir=bundle_dir))


def get_native_grid_lat_band_width(*, bundle_dir: Path) -> pint.Quantity:
    """
    Get the width of the latitudinal bands on our native grid

    Parameters
    ----------
    bundle_dir
        Directory in which to keep the original run's bundle

    Returns
    -------
    :
        Width of each latitudinal band

    Raises
    ------
    AssertionError
        The bands aren't all the same width, or don't cover the globe
    """
    # Doesn't matter which gas we get, they're all on the same grid
    lat = get_output("ch4", bundle_dir=bundle_dir)[0]["lat"].to_numpy()
    widths = np.unique(np.diff(lat))
    if widths.size != 1:
        msg = f"Latitudinal bands aren't all the same width: {widths=}"
        raise AssertionError(msg)

    width = float(widths[0])
    if not np.isclose(width * lat.size, 180.0):
        msg = f"Latitudinal bands don't cover the globe: {lat=}"
        raise AssertionError(msg)

    return Q(width, "degree")


def get_binning_grid_lon_band_width(*, bundle_dir: Path) -> pint.Quantity:
    """
    Get the width of the longitudinal bands we bin the observations into

    Read from every gas' interpolated observation network,
    which is on the grid the observations are binned into.

    Parameters
    ----------
    bundle_dir
        Directory in which to keep the original run's bundle

    Returns
    -------
    :
        Width of each longitudinal band

    Raises
    ------
    FileNotFoundError
        No interpolated observation network files were found in `bundle_dir`

    AssertionError
        The bands aren't all the same width, don't go round the globe
        or differ between gases
    """
    files = sorted(
        (bundle_dir / "data" / "interim").glob(
            "*/*_observational-network_interpolated.nc"
        )
    )
    if not files:
        msg = f"No interpolated observation network files found in {bundle_dir}"
        raise FileNotFoundError(msg)

    widths = set()
    for file in files:
        lon = xr.load_dataarray(file)["lon"].to_numpy()
        file_widths = np.unique(np.diff(lon))
        if file_widths.size != 1:
            msg = f"Longitudinal bands aren't all the same width in {file}: {lon=}"
            raise AssertionError(msg)

        width = float(file_widths[0])
        if not np.isclose(width * lon.size, 360.0):
            msg = f"Longitudinal bands don't go round the globe in {file}: {lon=}"
            raise AssertionError(msg)

        widths.add(width)

    if len(widths) != 1:
        msg = f"Longitudinal band widths differ between gases: {widths=}"
        raise AssertionError(msg)

    return Q(widths.pop(), "degree")


@cache
def get_obs_network_years(gas: str, *, bundle_dir: Path) -> tuple[int, int]:
    """
    Get the years covered by the observation network

    Parameters
    ----------
    gas
        Gas of interest

    bundle_dir
        Directory in which to keep the original run's bundle

    Returns
    -------
    :
        First and last year of the observation network's global-, annual-mean
    """
    years = xr.load_dataset(
        bundle_dir
        / "data"
        / "interim"
        / gas
        / f"{gas}_observational-network_global-annual-mean.nc"
    )["year"]

    return int(years.min()), int(years.max())


def get_seasonality_diff_from_observed_average(
    gas: str, *, bundle_dir: Path
) -> pint.Quantity:
    """
    Get how far our seasonality's average is from the observed average seasonality

    Both are averaged over the observation network period.
    The observed average seasonality is the relative seasonality
    multiplied by the observation network's average global-, annual-mean
    (which undoes the division used to calculate the relative seasonality).

    Parameters
    ----------
    gas
        Gas of interest

    bundle_dir
        Directory in which to keep the original run's bundle

    Returns
    -------
    :
        Maximum absolute difference,
        relative to the maximum absolute observed average seasonality
    """
    gas_dir = bundle_dir / "data" / "interim" / gas
    relative_seasonality = xr.load_dataarray(
        gas_dir / f"{gas}_observational-network_seasonality.nc"
    )
    obs_network_global_annual_mean = xr.load_dataarray(
        gas_dir / f"{gas}_observational-network_global-annual-mean.nc"
    )
    seasonality = xr.load_dataarray(
        gas_dir / f"{gas}_seasonality_fifteen-degree_allyears-monthly.nc"
    )

    observed = relative_seasonality * float(obs_network_global_annual_mean.mean())
    ours = seasonality.sel(year=obs_network_global_annual_mean["year"]).mean("year")

    return Q(
        float(np.abs(ours - observed).max() / np.abs(observed).max()), "dimensionless"
    )


def get_max_lat_gradient_eof_area_weighted_mean(*, bundle_dir: Path) -> pint.Quantity:
    """
    Get the largest area-weighted mean of any gas' latitudinal gradient EOFs

    Uses the EOFs our output is built from,
    i.e. the ones the original run kept, for every gas.
    The weights are cos(latitude), as in the original run's global-mean.
    On our grid of equal-width latitudinal bands,
    these are exactly proportional to each band's area
    (sin(c + h) - sin(c - h) = 2 cos(c) sin(h)).

    Parameters
    ----------
    bundle_dir
        Directory in which to keep the original run's bundle

    Returns
    -------
    :
        Largest absolute area-weighted mean of any EOF,
        relative to the largest absolute value of that gas' EOFs

    Raises
    ------
    FileNotFoundError
        No EOF files were found in `bundle_dir`
    """
    files = sorted(
        (bundle_dir / "data" / "interim").glob("*/*_allyears-lat-gradient-eofs-pcs.nc")
    )
    if not files:
        msg = f"No latitudinal gradient EOF files found in {bundle_dir}"
        raise FileNotFoundError(msg)

    res = 0.0
    for file in files:
        eofs = xr.load_dataset(file)["eofs"]
        weights = np.cos(np.deg2rad(eofs["lat"]))
        weighted_mean = (eofs * weights).sum("lat") / weights.sum()
        res = max(res, float(np.abs(weighted_mean).max() / np.abs(eofs).max()))

    return Q(res, "dimensionless")


def get_max_abs_diff_from_cmip6_obs_network(
    gas: str, *, bundle_dir: Path
) -> pint.Quantity:
    """
    Get the maximum absolute difference from CMIP6 over the observation network era

    Parameters
    ----------
    gas
        Gas of interest

    bundle_dir
        Directory in which to keep the original run's bundle

    Returns
    -------
    :
        Maximum absolute difference from CMIP6
        from the first year of the observation network onwards
    """
    return get_max_abs_diff_from_cmip6(
        gas,
        start=get_obs_network_years(gas, bundle_dir=bundle_dir)[0],
        bundle_dir=bundle_dir,
    )


def get_noaa_global_mean(gas: str) -> ComparisonTimeseries:
    """
    Get NOAA's global-mean monthly record (seasonal cycle included)

    Parameters
    ----------
    gas
        Gas of interest

    Returns
    -------
    :
        NOAA's global-mean monthly record
    """
    (res,) = (
        c for c in get_noaa_comparisons(gas, deseasonalised=False) if not c.is_spatial
    )

    return res


@cache
def get_diff_from_noaa_monthly(gas: str, *, bundle_dir: Path) -> pd.Series[float]:
    """
    Get the difference between our global-mean monthly and NOAA's

    Parameters
    ----------
    gas
        Gas of interest

    bundle_dir
        Directory in which to keep the original run's bundle

    Returns
    -------
    :
        Difference, for each month both cover, in our units
    """
    ours = monthly_to_series(get_output(gas, bundle_dir=bundle_dir)[2])
    noaa = comparison_to_series(
        get_noaa_global_mean(gas), get_units(gas, bundle_dir=bundle_dir), yearly=False
    )

    return (ours - noaa).dropna()


def get_diff_from_mauna_loa(*, bundle_dir: Path) -> pd.Series[float]:
    """
    Get the difference between our CO2 and NOAA's Mauna Loa record

    Our native resolution output is linearly interpolated to Mauna Loa's latitude.

    Parameters
    ----------
    bundle_dir
        Directory in which to keep the original run's bundle

    Returns
    -------
    :
        Difference, for each month both cover
    """
    gas = "co2"
    ours = monthly_to_series(
        get_output(gas, bundle_dir=bundle_dir)[0].interp(lat=MAUNA_LOA_LATITUDE)
    )
    (mauna_loa,) = (
        c for c in get_noaa_comparisons(gas, deseasonalised=False) if c.is_spatial
    )

    return (
        ours
        - comparison_to_series(
            mauna_loa, get_units(gas, bundle_dir=bundle_dir), yearly=False
        )
    ).dropna()


def get_diff_from_igcc(
    gas: str, start: int | None = None, *, bundle_dir: Path
) -> pd.Series[float]:
    """
    Get the difference between our global-, annual-mean and IGCC's

    Parameters
    ----------
    gas
        Gas of interest

    start
        First year to consider

        If `None`, all the years both cover.

    bundle_dir
        Directory in which to keep the original run's bundle

    Returns
    -------
    :
        Difference, by year, in our units
    """
    ours = get_output(gas, bundle_dir=bundle_dir)[1].to_series()
    igcc = comparison_to_series(
        get_igcc_comparison(gas), get_units(gas, bundle_dir=bundle_dir), yearly=True
    )

    return (ours - igcc).dropna().loc[start:]


def get_max_abs_diff_from_igcc_noaa_era(gas: str, *, bundle_dir: Path) -> pint.Quantity:
    """
    Get the max absolute difference from IGCC over the years NOAA's global-mean covers

    IGCC's record is based on NOAA's over these years
    (alone for CO2, averaged with AGAGE for CH4 and N2O from 2019 on),
    which is what the statements about IGCC are about.

    Parameters
    ----------
    gas
        Gas of interest

    bundle_dir
        Directory in which to keep the original run's bundle

    Returns
    -------
    :
        Maximum absolute difference from IGCC
    """
    noaa_start = int(np.floor(get_noaa_global_mean(gas).data[TIME_COLUMN].min()))

    return max_abs(
        get_diff_from_igcc(gas, start=noaa_start, bundle_dir=bundle_dir),
        get_units(gas, bundle_dir=bundle_dir),
    )


def get_mean_diff_from_igcc_noaa_era(gas: str, *, bundle_dir: Path) -> pint.Quantity:
    """
    Get the mean difference from IGCC over the years NOAA's global-mean covers

    Parameters
    ----------
    gas
        Gas of interest

    bundle_dir
        Directory in which to keep the original run's bundle

    Returns
    -------
    :
        Mean difference from IGCC
    """
    noaa_start = int(np.floor(get_noaa_global_mean(gas).data[TIME_COLUMN].min()))

    return Q(
        float(get_diff_from_igcc(gas, start=noaa_start, bundle_dir=bundle_dir).mean()),
        get_units(gas, bundle_dir=bundle_dir),
    )


def get_mean_diff_from_uci_ch4(*, bundle_dir: Path) -> pint.Quantity:
    """
    Get the mean difference between our global-mean CH4 and UCI's

    UCI's record is quarterly, so our global-mean monthly
    is linearly interpolated to UCI's times.

    Parameters
    ----------
    bundle_dir
        Directory in which to keep the original run's bundle

    Returns
    -------
    :
        Mean difference
    """
    gas = "ch4"
    ours = monthly_to_series(get_output(gas, bundle_dir=bundle_dir)[2])
    uci = (
        get_uci_ch4_comparison(deseasonalised=False)
        .to_units(get_units(gas, bundle_dir=bundle_dir))
        .data
    )
    in_range = uci[
        (uci[TIME_COLUMN] >= ours.index.min()) & (uci[TIME_COLUMN] <= ours.index.max())
    ]
    ours_at_uci_times = np.interp(in_range[TIME_COLUMN], ours.index, ours.to_numpy())

    return Q(
        float(np.mean(ours_at_uci_times - in_range[VALUE_COLUMN])),
        get_units(gas, bundle_dir=bundle_dir),
    )


@cache
def get_cfc12_like_obs_network_diffs(gas: str, *, bundle_dir: Path) -> pd.DataFrame:
    """
    Get the difference between our output and each observation in the network

    Our native resolution output is linearly interpolated
    to each observation's latitude.

    Parameters
    ----------
    gas
        Gas of interest

    bundle_dir
        Directory in which to keep the original run's bundle

    Returns
    -------
    :
        One row per observation, with its year, network and the difference
    """
    obs = get_cfc12_like_all_data_with_bins(gas, bundle_dir)
    (obs_units,) = obs["unit"].unique()
    if obs_units != get_units(gas, bundle_dir=bundle_dir):
        raise NotImplementedError(obs_units)

    native = get_output(gas, bundle_dir=bundle_dir)[0]
    res = []
    for latitude, latitude_obs in obs.groupby("latitude"):
        ours = native.interp(lat=latitude).to_series()
        index = pd.MultiIndex.from_frame(latitude_obs[["year", "month"]])
        res.append(
            pd.DataFrame(
                {
                    "year": latitude_obs["year"].to_numpy(),
                    "network": latitude_obs["network"].to_numpy(),
                    "diff": ours.reindex(index).to_numpy()
                    - latitude_obs["value"].to_numpy(),
                }
            )
        )

    return pd.concat(res).dropna()


def get_mean_diff_from_network(
    gas: str, network: str, *, bundle_dir: Path
) -> pint.Quantity:
    """
    Get the mean difference between our output and a network's observations

    Parameters
    ----------
    gas
        Gas of interest

    network
        Network of interest

    bundle_dir
        Directory in which to keep the original run's bundle

    Returns
    -------
    :
        Mean difference over all the network's observations
    """
    diffs = get_cfc12_like_obs_network_diffs(gas, bundle_dir=bundle_dir)

    return Q(
        float(diffs.loc[diffs["network"] == network, "diff"].mean()),
        get_units(gas, bundle_dir=bundle_dir),
    )


def is_global_mean_from_source(gas: str, source: str, *, bundle_dir: Path) -> bool:
    """
    Get whether a gas' global-, annual-mean comes from a source, rather than the network

    Parameters
    ----------
    gas
        Gas of interest

    source
        Label of the source

    bundle_dir
        Directory in which to keep the original run's bundle

    Returns
    -------
    :
        Whether `source` replaces the observation network's global-, annual-mean
    """
    supplement = get_global_mean_supplement(gas, bundle_dir)
    if supplement is None:
        return False

    label, supplement_data = supplement
    max_year_extended = int(
        xr.load_dataset(
            interim_dir(gas, bundle_dir) / f"{gas}_global-annual-mean_allyears.nc"
        )["year"].max()
    )

    return label == source and supplement_replaces_obs_network(
        supplement_data, max_year_extended
    )


TRUDINGER_GASES = ("cf4", "c2f6", "c3f8")
"""Gases whose global-, annual-mean includes Trudinger et al. (2016) data"""


HARMONISATION_TRANSITION_YEARS = 100
"""
Years over which a harmonisation offset declines to zero

Mirrors `n_transition_years=100`, which the original run uses for every harmonisation:
Menking et al. (2025) for N2O (`1004_n2o_extend-global-annual-mean`),
Law Dome for CH4 (`1104_ch4_extend-global-annual-mean`),
the Mauna Loa - Law Dome merged record and Menking et al. (2025) for CO2
(`1204_co2_extend-global-annual-mean`)
and Trudinger et al. (2016) for CF4, C2F6 and C3F8
(`1304_sf6-like_create-global-annual-mean`).
"""


def get_pre_industrial_years(source: str, *, bundle_dir: Path) -> pint.Quantity:
    """
    Get the pre-industrial year(s) of the gases processed like CFC-12 with a source

    Parameters
    ----------
    source
        Pre-industrial source, as the original run's config names it

    bundle_dir
        Directory in which to keep the original run's bundle

    Returns
    -------
    :
        Pre-industrial year of every gas with this source

    Raises
    ------
    AssertionError
        No gas has this source
    """
    years = [
        step_config["pre_industrial"]["year"]
        for step_config in (
            get_step_config(gas, bundle_dir) for gas in CFC12_LIKE_GASES
        )
        if step_config["pre_industrial"]["source"] == source
    ]
    if not years:
        msg = f"No gas has pre-industrial {source=}"
        raise AssertionError(msg)

    return Q(np.array(years), "yr")


def get_trudinger_years(*, bundle_dir: Path) -> pint.Quantity:
    """
    Get the first and last year of the Trudinger et al. (2016) data we use

    Parameters
    ----------
    bundle_dir
        Directory in which to keep the original run's bundle

    Returns
    -------
    :
        First and last year

    Raises
    ------
    AssertionError
        The years differ between [TRUDINGER_GASES][]
    """
    res = set()
    for gas in TRUDINGER_GASES:
        supplement = get_global_mean_supplement(gas, bundle_dir)
        if supplement is None:
            msg = f"No global-mean source for {gas=}"
            raise AssertionError(msg)

        years = supplement[1]["year"]
        res.add((int(years.min()), int(years.max())))

    if len(res) != 1:
        msg = f"Trudinger et al. (2016)'s years differ between gases: {res=}"
        raise AssertionError(msg)

    return Q(np.array(res.pop()), "yr")


def has_expected_non_zero_pre_industrial_gases(*, bundle_dir: Path) -> bool:
    """
    Get whether exactly the gases with natural sources have a non-zero pre-industrial

    Of the gases processed like CFC-12, these are CF4, CH2Cl2, CH3Br, CH3Cl and CHCl3
    (following M17).

    Parameters
    ----------
    bundle_dir
        Directory in which to keep the original run's bundle

    Returns
    -------
    :
        Whether exactly these gases have a non-zero pre-industrial value
    """
    non_zero = {
        gas
        for gas in CFC12_LIKE_GASES
        if get_step_config(gas, bundle_dir)["pre_industrial"]["value"][0] > 0.0
    }

    return non_zero == {"cf4", "ch2cl2", "ch3br", "ch3cl", "chcl3"}


VELDERS_SOURCE = "Velders et al., 2022"
"""Velders et al. (2022) pre-industrial source, as the original run's config names it"""


def get_velders_first_year_values(*, bundle_dir: Path) -> pd.Series[float]:
    """
    Get Velders et al. (2022)'s values in their first year

    For the gases which take their pre-industrial value from Velders et al. (2022)
    (without adjustment).

    Parameters
    ----------
    bundle_dir
        Directory in which to keep the original run's bundle

    Returns
    -------
    :
        Value for each gas in Velders et al. (2022)'s first year (in ppt),
        named by the first year
    """
    gases = [
        gas
        for gas in CFC12_LIKE_GASES
        if get_step_config(gas, bundle_dir)["pre_industrial"]["source"]
        == VELDERS_SOURCE
    ]
    velders = pd.read_csv(
        bundle_dir
        / "data"
        / "interim"
        / "velders-et-al-2022"
        / "velders_et_al_2022.csv"
    )
    (unit,) = velders["unit"].unique()
    if unit != "ppt":
        raise NotImplementedError(unit)

    first_year = int(velders["year"].min())
    res = velders[(velders["year"] == first_year) & velders["gas"].isin(gases)]
    if set(res["gas"]) != set(gases):
        msg = f"Velders et al. (2022) is missing some of {gases=}"
        raise AssertionError(msg)

    return res.set_index("gas")["value"].rename(first_year)


def is_velders_first_year_zero_except(exception: str, *, bundle_dir: Path) -> bool:
    """
    Get whether Velders et al. (2022)'s first-year value is zero for all but one gas

    Parameters
    ----------
    exception
        Gas whose value isn't zero

    bundle_dir
        Directory in which to keep the original run's bundle

    Returns
    -------
    :
        Whether the value is zero for every gas except `exception`
        (and non-zero for `exception`)
    """
    values = get_velders_first_year_values(bundle_dir=bundle_dir)

    return bool((values.drop(exception) == 0.0).all() and values[exception] > 0.0)


DROSTE_SITE_HEMISPHERES = {
    "cape-grim": -1,
    "tacolneston": 1,
}
"""Hemisphere of each site in Droste et al. (2020) (-1 south, 1 north)

The data identifies the sites by latitude only,
and the two are in different hemispheres.
"""


def load_droste(*, bundle_dir: Path) -> pd.DataFrame:
    """
    Load the Droste et al. (2020) data, as the original run processed it

    Parameters
    ----------
    bundle_dir
        Directory in which to keep the original run's bundle

    Returns
    -------
    :
        Droste et al. (2020) data, for all of [C4F10_LIKE_GASES][]
    """
    res = pd.read_csv(
        bundle_dir / "data" / "interim" / "droste-et-al-2020" / "droste_et_al_2020.csv"
    )
    if set(res["gas"]) != set(C4F10_LIKE_GASES):
        msg = f"Expected data for {C4F10_LIKE_GASES=}, found {set(res['gas'])}"
        raise AssertionError(msg)

    return res


def get_droste_years(*, bundle_dir: Path) -> tuple[int, int]:
    """
    Get the first and last year of the Droste et al. (2020) data

    Parameters
    ----------
    bundle_dir
        Directory in which to keep the original run's bundle

    Returns
    -------
    :
        First and last year

    Raises
    ------
    AssertionError
        The years differ between gases or sites
    """
    droste = load_droste(bundle_dir=bundle_dir)
    years = droste.groupby(["gas", "lat"])["year"].agg(["min", "max"])
    if len(years.drop_duplicates()) != 1:
        msg = f"Droste et al. (2020)'s years differ between gases or sites: {years}"
        raise AssertionError(msg)

    return int(years["min"].iloc[0]), int(years["max"].iloc[0])


def get_droste_max_first_year_value(*, bundle_dir: Path) -> pint.Quantity:
    """
    Get the largest Droste et al. (2020) value in its first year

    Over all gases and both sites.

    Parameters
    ----------
    bundle_dir
        Directory in which to keep the original run's bundle

    Returns
    -------
    :
        Largest value in the first year
    """
    droste = load_droste(bundle_dir=bundle_dir)
    (unit,) = droste["unit"].unique()
    first_year = get_droste_years(bundle_dir=bundle_dir)[0]

    return Q(float(droste.loc[droste["year"] == first_year, "value"].max()), unit)


def get_droste_site_lat(site: str, *, bundle_dir: Path) -> pint.Quantity:
    """
    Get the latitude of a site in Droste et al. (2020)

    Parameters
    ----------
    site
        Site of interest (a key of [DROSTE_SITE_HEMISPHERES][])

    bundle_dir
        Directory in which to keep the original run's bundle

    Returns
    -------
    :
        Latitude of `site`

    Raises
    ------
    AssertionError
        The data doesn't have exactly one latitude in `site`'s hemisphere
    """
    lats = load_droste(bundle_dir=bundle_dir)["lat"].unique()
    (lat,) = (v for v in lats if np.sign(v) == DROSTE_SITE_HEMISPHERES[site])

    return Q(float(lat), "degree")


def get_c4f10_like_last_year(*, bundle_dir: Path) -> pint.Quantity:
    """
    Get the last year of our output for the gases processed like C4F10

    Parameters
    ----------
    bundle_dir
        Directory in which to keep the original run's bundle

    Returns
    -------
    :
        Last year

    Raises
    ------
    AssertionError
        The last year differs between gases
    """
    last_years = {
        int(get_output(gas, bundle_dir=bundle_dir)[1]["year"].max())
        for gas in C4F10_LIKE_GASES
    }
    if len(last_years) != 1:
        msg = f"Last year differs between gases: {last_years=}"
        raise AssertionError(msg)

    return Q(last_years.pop(), "yr")


def get_c4f10_like_erf(year: int, *, bundle_dir: Path) -> pint.Quantity:
    """
    Get the approximate ERF of all the gases processed like C4F10 together

    Each gas' ERF is its concentration change since 1750
    multiplied by its radiative efficiency
    (see [get_radiative_effect][]).

    Parameters
    ----------
    year
        Year of interest

    bundle_dir
        Directory in which to keep the original run's bundle

    Returns
    -------
    :
        Approximate ERF of the gases processed like C4F10 in `year`
    """
    res = Q(0.0, "W / m^2")
    for gas in C4F10_LIKE_GASES:
        global_annual_mean = get_output(gas, bundle_dir=bundle_dir)[1]
        change = Q(
            float(
                global_annual_mean.sel(year=year) - global_annual_mean.sel(year=1750)
            ),
            get_units(gas, bundle_dir=bundle_dir),
        )
        res += (change * RADIATIVE_EFFICIENCIES[gas]).to("W / m^2")

    return res


MAX_LAT_GRADIENT_FRACTION = 0.5
"""
Largest magnitude of the latitudinal gradient's most negative value

As a fraction of the global-mean.

Mirrors the `0.5 * month_da` in the original run's
`1305_sf6-like_create-pieces-for-gridding`
and `1405_c4f10-like_create-pieces-for-gridding` notebooks.
"""

MAX_SEASONALITY_FRACTION = 0.35
"""
Largest the seasonality can be, as a fraction of the global-mean

Mirrors `max_reasonable_seasonality_frac = 0.35`
in the original run's `1305_sf6-like_create-pieces-for-gridding` notebook.
"""


def load_monthly_pieces(
    gas: str, *, bundle_dir: Path
) -> tuple[xr.DataArray, xr.DataArray, xr.DataArray]:
    """
    Load the monthly pieces our native resolution output is built from

    Parameters
    ----------
    gas
        Gas of interest

    bundle_dir
        Directory in which to keep the original run's bundle

    Returns
    -------
    :
        Global-, annual-mean (interpolated to monthly steps),
        latitudinal gradient and seasonality,
        over the years they all cover and which are in the published output
        (the pieces can run beyond the last year which is written)
    """
    gas_dir = interim_dir(gas, bundle_dir)
    last_year = int(get_written_year("end_year", bundle_dir=bundle_dir).m)

    pieces = xr.align(
        xr.load_dataarray(gas_dir / f"{gas}_global-annual-mean_allyears-monthly.nc"),
        xr.load_dataarray(
            gas_dir / f"{gas}_latitudinal-gradient_fifteen-degree_allyears-monthly.nc"
        ),
        xr.load_dataarray(
            gas_dir / f"{gas}_seasonality_fifteen-degree_allyears-monthly.nc"
        ),
        join="inner",
    )

    return tuple(piece.sel(year=slice(None, last_year)) for piece in pieces)  # type: ignore[return-value]


def get_lat_gradient_capped_years(gas: str, *, bundle_dir: Path) -> np.ndarray:
    """
    Get the years in which the latitudinal gradient was scaled down

    In these months, the latitudinal gradient's most negative value
    is exactly [MAX_LAT_GRADIENT_FRACTION][] of the global-mean.

    Parameters
    ----------
    gas
        Gas of interest

    bundle_dir
        Directory in which to keep the original run's bundle

    Returns
    -------
    :
        Years with at least one month in which the latitudinal gradient was scaled down
    """
    global_mean, lat_gradient, _ = load_monthly_pieces(gas, bundle_dir=bundle_dir)
    capped = (global_mean > 0.0) & (
        np.abs(lat_gradient.min("lat") + MAX_LAT_GRADIENT_FRACTION * global_mean)
        <= 1e-6 * np.abs(global_mean)
    )

    return capped["year"].to_numpy()[capped.any("month").to_numpy()]


def get_seasonality_capped_years(gas: str, *, bundle_dir: Path) -> np.ndarray:
    """
    Get the years in which the seasonality was scaled down

    In these years, the seasonality's largest magnitude
    is exactly [MAX_SEASONALITY_FRACTION][] of the global-mean.

    Parameters
    ----------
    gas
        Gas of interest

    bundle_dir
        Directory in which to keep the original run's bundle

    Returns
    -------
    :
        Years in which the seasonality was scaled down
    """
    global_mean, _, seasonality = load_monthly_pieces(gas, bundle_dir=bundle_dir)
    fraction = (np.abs(seasonality) / global_mean.where(global_mean > 0.0)).max(
        ["lat", "month"]
    )

    return fraction["year"].to_numpy()[
        np.isclose(fraction.to_numpy(), MAX_SEASONALITY_FRACTION, rtol=1e-6, atol=0.0)
    ]


def get_capped_years_extent(
    years_getter: Callable[..., np.ndarray],
    gases: tuple[str, ...],
    *,
    bundle_dir: Path,
) -> pint.Quantity:
    """
    Get the first and last year in which any of a group of gases was scaled down

    Parameters
    ----------
    years_getter
        Gets the scaled-down years for a gas
        (e.g. [get_lat_gradient_capped_years][])

    gases
        Gases of interest

    bundle_dir
        Directory in which to keep the original run's bundle

    Returns
    -------
    :
        First and last year

    Raises
    ------
    AssertionError
        None of the gases was scaled down in any year
    """
    years = np.concatenate([years_getter(gas, bundle_dir=bundle_dir) for gas in gases])
    if years.size == 0:
        msg = f"None of {gases=} was scaled down"
        raise AssertionError(msg)

    return Q(np.array([int(years.min()), int(years.max())]), "yr")


def get_max_lat_gradient_capped_years_after_pre_industrial(
    exclude: tuple[str, ...], *, bundle_dir: Path
) -> pint.Quantity:
    """
    Get how long after the pre-industrial year the latitudinal gradient is scaled down

    Over the gases processed like CFC-12, excluding some.

    Parameters
    ----------
    exclude
        Gases to exclude

    bundle_dir
        Directory in which to keep the original run's bundle

    Returns
    -------
    :
        Largest number of years between a gas' pre-industrial year
        and the last year in which its latitudinal gradient was scaled down
    """
    res = 0
    for gas in CFC12_LIKE_GASES:
        if gas in exclude:
            continue

        years = get_lat_gradient_capped_years(gas, bundle_dir=bundle_dir)
        if years.size == 0:
            continue

        pre_industrial_year = get_step_config(gas, bundle_dir)["pre_industrial"]["year"]
        res = max(res, int(years.max()) - pre_industrial_year)

    return Q(res, "yr")


def get_seasonality_capped_single_year_gases(*, bundle_dir: Path) -> bool:
    """
    Get whether the seasonality is scaled down in the gases and years the methods say

    I.e. in every year from its first for HFC-236fa,
    in a single year for HFC-32, HFC-152a and HFC-365mfc
    and never for any other gas processed like CFC-12.

    Parameters
    ----------
    bundle_dir
        Directory in which to keep the original run's bundle

    Returns
    -------
    :
        Whether this is the case
    """
    n_years = {
        gas: get_seasonality_capped_years(gas, bundle_dir=bundle_dir).size
        for gas in CFC12_LIKE_GASES
    }
    single = {gas for gas, n in n_years.items() if n == 1}
    many = {gas for gas, n in n_years.items() if n > 1}

    return single == {"hfc32", "hfc152a", "hfc365mfc"} and many == {"hfc236fa"}


def get_co2_seasonality_change_regression_start_year(
    *, bundle_dir: Path
) -> pint.Quantity:
    """
    Get the first year in which the CO2 seasonality change PC comes from the regression

    Before this year, the original run
    (`1205_co2_extend-seasonality-change-pcs`)
    keeps the regression's composite, and hence the PC, constant.
    So this is the last year in which the PC still has its year-1 value.

    Parameters
    ----------
    bundle_dir
        Directory in which to keep the original run's bundle

    Returns
    -------
    :
        First year of the regression
    """
    pc = xr.load_dataset(
        interim_dir("co2", bundle_dir) / "co2_allyears-seasonality-change-eofs-pcs.nc"
    )["principal-components"].sel(eof=0)
    values = pc.to_numpy()
    first_change = np.argmax(~np.isclose(values, values[0], rtol=0.0, atol=1e-12))

    return Q(int(pc["year"][first_change - 1]), "yr")


MIN_POINTS_FOR_SPATIAL_INTERPOLATION = 4
"""
Minimum number of binned values needed to spatially interpolate a month

Mirrors `MIN_POINTS_FOR_SPATIAL_INTERPOLATION` in the original run's
`*_interpolate-observational-network` notebooks (N2O, CH4, CO2 and SF6-like).
"""


CO2_MAUNA_LOA_START = 1959
"""
First year of the Mauna Loa - Law Dome merged record we use

Mirrors `mauna_loa_start` in the original run's `1204_co2_extend-global-annual-mean`
(the first full year of Mauna Loa data).
The Mauna Loa record isn't saved on its own, so this can't be read from the data.
"""


def get_variance_explained_fraction(
    decomposition: str, *, bundle_dir: Path
) -> pd.Series[float]:
    """
    Get the fraction of the variance each EOF of a decomposition explains

    These come from re-running the original run's notebooks
    (see [local.historical_ghg_forcing_for_cmip7.variance_explained][]),
    which the methods figures do.

    Parameters
    ----------
    decomposition
        Decomposition of interest, e.g. `"co2_lat-gradient"`
        or `"co2_seasonality-change"`

    bundle_dir
        Directory in which to keep the original run's bundle

    Returns
    -------
    :
        Fraction of the variance explained, indexed by EOF
    """
    return pd.read_csv(
        bundle_dir / "manuscript-outputs" / f"{decomposition}-variance-explained.csv",
        index_col="eof",
    )["variance_explained_fraction"]


def get_min_lat_gradient_first_eof_variance(
    exclude: tuple[str, ...], *, bundle_dir: Path
) -> pint.Quantity:
    """
    Get the smallest variance explained by the first EOF over the CFC12-like gases

    Parameters
    ----------
    exclude
        Gases to exclude

    bundle_dir
        Directory in which to keep the original run's bundle

    Returns
    -------
    :
        Smallest fraction of the variance explained by the first EOF
    """
    return Q(
        min(
            float(
                get_variance_explained_fraction(
                    f"{gas}_lat-gradient", bundle_dir=bundle_dir
                ).loc[0]
            )
            for gas in CFC12_LIKE_GASES
            if gas not in exclude
        ),
        "dimensionless",
    )


def get_annual_mean_at_lat(gas: str, lat: float, *, bundle_dir: Path) -> xr.DataArray:
    """
    Get the annual-mean of our native resolution output in a latitudinal bin

    Parameters
    ----------
    gas
        Gas of interest

    lat
        Latitude of interest (the bin whose centre is nearest is used)

    bundle_dir
        Directory in which to keep the original run's bundle

    Returns
    -------
    :
        Annual-mean in the bin, by year
    """
    native = get_output(gas, bundle_dir=bundle_dir)[0]

    return native.sel(lat=lat, method="nearest").mean("month")


def load_menking(gas: str, *, bundle_dir: Path) -> tuple[pd.Series[float], float]:
    """
    Load the Menking et al. (2025) data for a gas, as the original run processed it

    Parameters
    ----------
    gas
        Gas of interest

    bundle_dir
        Directory in which to keep the original run's bundle

    Returns
    -------
    :
        Values, by year, and the latitude of the data
    """
    menking = pd.read_csv(
        bundle_dir
        / "data"
        / "interim"
        / "menking-et-al-2025"
        / "menking_et_al_2025.csv"
    )
    menking = menking[menking["gas"] == gas]
    (lat,) = menking["latitude"].unique()

    return menking.set_index("year")["value"], float(lat)


def get_n2o_menking_offset(*, bundle_dir: Path) -> pint.Quantity:
    """
    Get the offset between our N2O global-, annual-mean and Menking et al. (2025)

    In the first year of the observation network period,
    which is where the original run harmonises them.

    Parameters
    ----------
    bundle_dir
        Directory in which to keep the original run's bundle

    Returns
    -------
    :
        Observation network-derived global-, annual-mean minus Menking et al. (2025)
    """
    menking, _ = load_menking("n2o", bundle_dir=bundle_dir)
    harmonisation_year = get_obs_network_years("n2o", bundle_dir=bundle_dir)[0]
    obs_network = xr.load_dataarray(
        interim_dir("n2o", bundle_dir)
        / "n2o_observational-network_global-annual-mean.nc"
    )

    return Q(
        float(obs_network.sel(year=harmonisation_year))
        - menking.loc[harmonisation_year],
        get_units("n2o", bundle_dir=bundle_dir),
    )


def get_co2_mauna_loa_offset(*, bundle_dir: Path) -> pint.Quantity:
    """
    Get the offset between our CO2 global-, annual-mean and the Mauna Loa merged record

    In the first year of the observation network period,
    which is where the original run harmonises them
    (using the merged record's mid-year value, as the original run does).

    Parameters
    ----------
    bundle_dir
        Directory in which to keep the original run's bundle

    Returns
    -------
    :
        Observation network-derived global-, annual-mean minus the merged record
    """
    merged = load_co2_mauna_loa_merged(bundle_dir=bundle_dir)
    harmonisation_year = get_obs_network_years("co2", bundle_dir=bundle_dir)[0]
    obs_network = xr.load_dataarray(
        interim_dir("co2", bundle_dir)
        / "co2_observational-network_global-annual-mean.nc"
    )

    return Q(
        float(obs_network.sel(year=harmonisation_year))
        - merged.loc[harmonisation_year + 0.5],
        get_units("co2", bundle_dir=bundle_dir),
    )


def load_co2_mauna_loa_merged(*, bundle_dir: Path) -> pd.Series[float]:
    """
    Load the Mauna Loa - Law Dome merged record, as the original run processed it

    Parameters
    ----------
    bundle_dir
        Directory in which to keep the original run's bundle

    Returns
    -------
    :
        Values, indexed by (decimal) time
    """
    return pd.read_csv(
        bundle_dir / "data" / "interim" / "mauna_loa" / "merged_ice_core.csv"
    ).set_index("time")["value"]


def get_co2_menking_offset_and_match(
    *, bundle_dir: Path
) -> tuple[pint.Quantity, pint.Quantity]:
    """
    Get the CO2 Menking et al. (2025) offset and how well our output then matches it

    The original run harmonises Menking et al. (2025)
    to our output in Menking et al. (2025)'s latitudinal bin
    in [CO2_MAUNA_LOA_START][],
    then sets our earlier output to match the harmonised data in that bin.

    Parameters
    ----------
    bundle_dir
        Directory in which to keep the original run's bundle

    Returns
    -------
    :
        Offset (our output minus Menking et al. (2025)) in the harmonisation year
        and the largest absolute difference between our output
        and the harmonised Menking et al. (2025) data before that year
    """
    menking, lat = load_menking("co2", bundle_dir=bundle_dir)
    ours = get_annual_mean_at_lat("co2", lat, bundle_dir=bundle_dir)
    harmonisation_year = CO2_MAUNA_LOA_START
    offset = float(ours.sel(year=harmonisation_year)) - menking.loc[harmonisation_year]

    years = np.arange(int(ours["year"].min()), harmonisation_year)
    harmonised = menking.loc[years].to_numpy() + offset * np.clip(
        (years - (harmonisation_year - HARMONISATION_TRANSITION_YEARS))
        / HARMONISATION_TRANSITION_YEARS,
        0.0,
        None,
    )
    units = get_units("co2", bundle_dir=bundle_dir)

    return (
        Q(offset, units),
        Q(float(np.abs(ours.sel(year=years).to_numpy() - harmonised).max()), units),
    )


def get_ch4_law_dome_offset(*, data_raw_dir: Path, bundle_dir: Path) -> pint.Quantity:
    """
    Get the size of the offset between our CH4 output and the smoothed Law Dome data

    In Law Dome's latitudinal bin,
    in the first year of the observation network period,
    which is where the original run harmonises them.

    Parameters
    ----------
    data_raw_dir
        Raw data directory

    bundle_dir
        Directory in which to keep the original run's bundle

    Returns
    -------
    :
        Absolute offset
    """
    law_dome = load_ch4_ice_core("law-dome", data_raw_dir)
    lat = get_ice_core_latitude(law_dome, "law-dome")
    harmonisation_year = get_obs_network_years("ch4", bundle_dir=bundle_dir)[0]
    ours = get_annual_mean_at_lat("ch4", lat, bundle_dir=bundle_dir)

    return Q(
        abs(
            float(ours.sel(year=harmonisation_year))
            - law_dome.set_index("year")["value"].loc[harmonisation_year]
        ),
        get_units("ch4", bundle_dir=bundle_dir),
    )


def get_ch4_max_neem_relative_diff(
    *, data_raw_dir: Path, bundle_dir: Path
) -> pint.Quantity:
    """
    Get how far our CH4 output is from the NEEM data, at most

    In NEEM's latitudinal bin, in the years of NEEM's observations
    (rounded to the nearest year, as the original run does).

    Parameters
    ----------
    data_raw_dir
        Raw data directory

    bundle_dir
        Directory in which to keep the original run's bundle

    Returns
    -------
    :
        Largest absolute difference, relative to NEEM
    """
    neem = load_ch4_ice_core("neem", data_raw_dir)
    lat = get_ice_core_latitude(neem, "neem")
    ours = get_annual_mean_at_lat("ch4", lat, bundle_dir=bundle_dir).sel(
        year=neem["year"].round(0).to_numpy()
    )
    neem_values = neem["value"].to_numpy()

    return Q(
        float(np.max(np.abs(ours.to_numpy() - neem_values) / neem_values)),
        "dimensionless",
    )


def get_ch4_pc0_optimised_years(*, bundle_dir: Path) -> tuple[int, int]:
    """
    Get the first and last year in which the CH4 first PC is optimised against ice cores

    These come from re-running the original run's `1103_ch4_extend-pcs` notebook,
    which the CH4 methods figure does.

    Parameters
    ----------
    bundle_dir
        Directory in which to keep the original run's bundle

    Returns
    -------
    :
        First and last year
    """
    years = json.loads(
        (bundle_dir / "manuscript-outputs" / "ch4_pc0-optimised-years.json").read_text()
    )

    return min(years), max(years)


def get_ch4_law_dome_smoothing_config(*, bundle_dir: Path) -> dict[str, Any]:
    """
    Get the original run's config for smoothing the CH4 Law Dome data

    Parameters
    ----------
    bundle_dir
        Directory in which to keep the original run's bundle

    Returns
    -------
    :
        The config, as loaded from the YAML

    Raises
    ------
    AssertionError
        The bundle's config has no CH4 Law Dome smoothing config
    """
    config = yaml.safe_load((bundle_dir / BUNDLE_CONFIG_FILE).read_text())
    for step_config in config["smooth_law_dome_data"]:
        if step_config["gas"] == "ch4":
            return step_config  # type: ignore[no-any-return]

    msg = "No CH4 Law Dome smoothing config"
    raise AssertionError(msg)


def get_input_lat_gradient_weakening(gas: str, *, bundle_dir: Path) -> pint.Quantity:
    """
    Get how fast the latitudinal gradient in the input data is weakening

    For each year, the latitudinal gradient is the gradient of a linear fit
    of the observations' annual-mean at each latitude against latitude.
    The rate of weakening is the trend in this gradient's magnitude
    over the last ten years of our output, with the sign flipped,
    i.e. positive if the gradient is weakening.

    Parameters
    ----------
    gas
        Gas of interest

    bundle_dir
        Directory in which to keep the original run's bundle

    Returns
    -------
    :
        Rate of weakening of the latitudinal gradient
    """
    obs = get_cfc12_like_all_data_with_bins(gas, bundle_dir)
    years = get_last_n_years(gas, bundle_dir=bundle_dir)
    annual_mean = (
        obs[obs["year"].isin(years)]
        .groupby(["year", "latitude"])["value"]
        .mean()
        .reset_index()
    )
    gradients = pd.Series(
        {
            year: np.polyfit(year_df["latitude"], year_df["value"], 1)[0]
            for year, year_df in annual_mean.groupby("year")
        }
    )
    trend = np.polyfit(gradients.index, gradients.abs().to_numpy(), 1)[0]

    return Q(-trend, f"{get_units(gas, bundle_dir=bundle_dir)} / degree / yr")


REPRODUCTION_TOLERANCES = {"CMIP7": 1e-9, "CMIP6": 1e-3, "IGCC": 1e-9}
"""
How closely each dataset's definition must reproduce its published equivalent species

The CMIP6 data is stored as 32-bit floats, hence the looser tolerance.
"""


def check_definition_reproduces(
    dataset: EquivalenceDataset, equivalent_species: str
) -> None:
    """
    Check that a dataset's equivalent species definition reproduces its published values

    If it does, we know we have copied the definition correctly
    (and, e.g., that no gas is double counted or missing).

    Parameters
    ----------
    dataset
        Dataset of interest

    equivalent_species
        Equivalent species of interest

    Raises
    ------
    AssertionError
        The definition doesn't reproduce the published values
    """
    residual = check_reproduction(equivalent_species, dataset).abs().max()
    tolerance = REPRODUCTION_TOLERANCES[dataset.label]
    if not residual < tolerance:
        msg = (
            f"{dataset.label}'s definition of {equivalent_species} "
            f"doesn't reproduce its published values: {residual=}, {tolerance=}"
        )
        raise AssertionError(msg)


@cache
def get_equivalent_species_decomposition(
    equivalent_species: str, other: Literal["CMIP6", "IGCC"], *, bundle_dir: Path
) -> pd.DataFrame:
    """
    Get the difference between our equivalent species and another dataset's, split up

    See
    [local.historical_ghg_forcing_for_cmip7.equivalent_species.decompose_difference][].

    Parameters
    ----------
    equivalent_species
        Equivalent species of interest

    other
        Dataset to compare against

    bundle_dir
        Directory in which to keep the original run's bundle

    Returns
    -------
    :
        Decomposition of the difference, by year
    """
    cmip7 = load_cmip7(bundle_dir)
    other_dataset = load_cmip6() if other == "CMIP6" else load_igcc()
    for dataset in (cmip7, other_dataset):
        check_definition_reproduces(dataset, equivalent_species)

    return decompose_difference(equivalent_species, cmip7, other_dataset)


def get_difference_share(
    equivalent_species: str,
    other: Literal["CMIP6", "IGCC"],
    part: str,
    year: Literal["max", "last"],
    *,
    bundle_dir: Path,
) -> pint.Quantity:
    """
    Get the share of the difference from another dataset due to one part

    Parameters
    ----------
    equivalent_species
        Equivalent species of interest

    other
        Dataset to compare against

    part
        Part of interest, one of
        [local.historical_ghg_forcing_for_cmip7.equivalent_species.DIFFERENCE_PARTS][]

    year
        Year to get the share in:
        the year of maximum absolute difference ("max")
        or the last year both datasets cover ("last")

    bundle_dir
        Directory in which to keep the original run's bundle

    Returns
    -------
    :
        Share of the difference due to `part`
    """
    decomposition = get_equivalent_species_decomposition(
        equivalent_species, other, bundle_dir=bundle_dir
    )
    if year == "max":
        selected_year = decomposition["total difference"].abs().idxmax()
    else:
        selected_year = decomposition.index.max()

    share = (
        decomposition.loc[selected_year, part]
        / decomposition.loc[selected_year, "total difference"]
    )

    return Q(float(share), "dimensionless")


def get_max_abs_diff_from_cmip6_recalculated(
    equivalent_species: str, *, bundle_dir: Path
) -> pint.Quantity:
    """
    Get the max abs difference from CMIP6 recalculated with AR6 radiative efficiencies

    I.e. the difference which is left
    once the difference in radiative efficiencies is removed.

    Parameters
    ----------
    equivalent_species
        Equivalent species of interest

    bundle_dir
        Directory in which to keep the original run's bundle

    Returns
    -------
    :
        Maximum absolute difference
    """
    decomposition = get_equivalent_species_decomposition(
        equivalent_species, "CMIP6", bundle_dir=bundle_dir
    )
    recalculated = (
        decomposition["total difference"]
        - decomposition["difference in radiative efficiency"]
    )

    return max_abs(recalculated, get_units(equivalent_species, bundle_dir=bundle_dir))


def get_diff_from_igcc_last_year(gas: str, *, bundle_dir: Path) -> pint.Quantity:
    """
    Get the difference from IGCC in the last year both datasets cover

    Parameters
    ----------
    gas
        Gas of interest

    bundle_dir
        Directory in which to keep the original run's bundle

    Returns
    -------
    :
        Difference from IGCC
    """
    return Q(
        float(get_diff_from_igcc(gas, bundle_dir=bundle_dir).iloc[-1]),
        get_units(gas, bundle_dir=bundle_dir),
    )


def radiative_effect_of(
    gas: str, calculate: Callable[[], pint.Quantity]
) -> pint.Quantity:
    """
    Get the approximate radiative effect of a calculated (difference in) concentration

    Parameters
    ----------
    gas
        Gas of interest

    calculate
        Calculates the (difference in) concentration

    Returns
    -------
    :
        Approximate radiative effect
    """
    return get_radiative_effect(gas, calculate())


def get_ch4_ice_core_lat(source: str, *, data_raw_dir: Path) -> pint.Quantity:
    """
    Get an ice core's latitude

    Parameters
    ----------
    source
        Ice core of interest (a key of [CH4_ICE_CORE_FILES][])

    data_raw_dir
        Raw data directory

    Returns
    -------
    :
        `source`'s latitude
    """
    lat = get_ice_core_latitude(load_ch4_ice_core(source, data_raw_dir), source)

    return Q(lat, "degree")


def get_ch4_law_dome_start_year(*, data_raw_dir: Path) -> pint.Quantity:
    """
    Get the first year of the smoothed Law Dome CH4 data

    The original run (`1104_ch4_extend-global-annual-mean`)
    uses EPICA for every year before this one.

    Parameters
    ----------
    data_raw_dir
        Raw data directory

    Returns
    -------
    :
        First year of the smoothed Law Dome CH4 data
    """
    return Q(int(load_ch4_ice_core("law-dome", data_raw_dir)["year"].min()), "yr")


def get_ch4_ice_core_lat_bin(
    source: str, *, bundle_dir: Path, data_raw_dir: Path
) -> pint.Quantity:
    """
    Get the native latitudinal bin which an ice core's data is matched in

    Mirrors the original run (`1104_ch4_extend-global-annual-mean`),
    which takes the bin whose centre is nearest the ice core's latitude.

    Parameters
    ----------
    source
        Ice core of interest (a key of [CH4_ICE_CORE_FILES][])

    bundle_dir
        Directory in which to keep the original run's bundle

    data_raw_dir
        Raw data directory

    Returns
    -------
    :
        Centre of the latitudinal bin
    """
    lat = get_ch4_ice_core_lat(source, data_raw_dir=data_raw_dir)

    # Doesn't matter which gas we get, they're all on the same grid
    native_output = get_output("ch4", bundle_dir=bundle_dir)[0]

    return Q(
        float(native_output["lat"].sel(lat=lat.to("degree").m, method="nearest")),
        "degree",
    )


def get_value_checks(  # noqa: PLR0912, PLR0915
    bundle_dir: Path = DEFAULT_BUNDLE_DIR, data_raw_dir: Path = DATA_RAW_DIR
) -> tuple[ValueCheck, ...]:
    """
    Get the value checks for the historical manuscript

    Parameters
    ----------
    bundle_dir
        Directory in which to keep the original run's bundle

    data_raw_dir
        Raw data directory

        Only used for the original run's files which aren't on Zenodo
        (see [local.historical_ghg_forcing_for_cmip7.zenodo_missing][]).
        The comparison data still comes from
        [local.historical_ghg_forcing_for_cmip7.comparison_data][]'s
        default locations.

    Returns
    -------
    :
        Value checks
    """
    res: list[ValueCheck] = []

    def add(tag: str, calculate: Callable[[], CheckValue], description: str) -> None:
        res.append(ValueCheck(tag=tag, calculate=calculate, description=description))

    def add_with_radiative_effect(
        tag: str, gas: str, calculate: Callable[[], pint.Quantity], description: str
    ) -> None:
        add(tag, calculate, description)
        add(
            f"{tag}-radiative-effect",
            partial(radiative_effect_of, gas, calculate),
            f"approx. radiative effect of the {description}",
        )

    seasonal_cycle_regions = {
        "co2": (NORTHERN_HEMISPHERE,),
        "ch4": (NORTHERN_HEMISPHERE, SOUTHERN_HEMISPHERE),
        "n2o": (NORTHERN_HEMISPHERE, SOUTHERN_HEMISPHERE),
        "cfc12": (NORTHERN_HEMISPHERE, SOUTHERN_HEMISPHERE),
    }
    for gas in ("co2", "ch4", "n2o", "cfc12", "cfc12eq", "hfc134aeq"):
        add(
            f"{gas}-last-10-years-trend",
            partial(get_last_n_years_trend, gas, bundle_dir=bundle_dir),
            "linear trend in the global-, annual-mean over the last ten years",
        )
        add(
            f"{gas}-last-10-years-second-derivative",
            partial(get_last_n_years_second_derivative, gas, bundle_dir=bundle_dir),
            "second derivative of the global-, annual-mean over the last ten years",
        )
        if gas in seasonal_cycle_regions:
            add(
                f"{gas}-last-10-years-lat-grad",
                partial(get_last_n_years_lat_gradient, gas, bundle_dir=bundle_dir),
                "latitudinal gradient (linear fit) over the last ten years",
            )
            regions = seasonal_cycle_regions[gas]
            add(
                f"{gas}-last-10-years-seasonal-cycle",
                partial(
                    get_last_n_years_seasonal_cycle, gas, regions, bundle_dir=bundle_dir
                ),
                "mean peak-to-trough seasonal cycle over the last ten years "
                f"({', '.join(r.lower() for r in regions)})",
            )

        add_with_radiative_effect(
            f"{gas}-diff-from-cmip6-all",
            gas,
            partial(get_max_abs_diff_from_cmip6, gas, bundle_dir=bundle_dir),
            "max abs difference from CMIP6 global-, annual-mean, all years",
        )
        add(
            f"{gas}-diff-from-cmip6-all-year",
            partial(get_year_of_max_abs_diff_from_cmip6, gas, bundle_dir=bundle_dir),
            "year of max abs difference from CMIP6 global-, annual-mean",
        )
        add_with_radiative_effect(
            f"{gas}-diff-from-cmip6-1850-on",
            gas,
            partial(get_max_abs_diff_from_cmip6, gas, 1850, bundle_dir=bundle_dir),
            "max abs difference from CMIP6 global-, annual-mean, 1850 onwards",
        )
        add(
            f"{gas}-diff-from-cmip6-1850-on-year",
            partial(
                get_year_of_max_abs_diff_from_cmip6, gas, 1850, bundle_dir=bundle_dir
            ),
            "year of max abs difference from CMIP6 global-, annual-mean, 1850 onwards",
        )

    add(
        "co2-for-energy-balance-threshold",
        lambda: (ENERGY_BALANCE_THRESHOLD / RADIATIVE_EFFICIENCIES["co2"]).to("ppm"),
        f"CO2 concentration change with an approx. radiative effect of "
        f"{ENERGY_BALANCE_THRESHOLD:~}",
    )
    add(
        "output-first-year",
        partial(get_written_year, "start_year", bundle_dir=bundle_dir),
        "first year for which the output is written (the same for every gas)",
    )
    add(
        "output-last-year",
        partial(get_written_year, "end_year", bundle_dir=bundle_dir),
        "last year for which the output is written (the same for every gas)",
    )
    for gas in ("co2", "ch4", "n2o", "cfc12"):
        add(
            f"{gas}-last-year",
            partial(get_last_year, gas, bundle_dir=bundle_dir),
            "last year of our global-, annual-mean",
        )
        add(
            f"{gas}-cmip6-last-year",
            lambda gas=gas: Q(
                int(get_diff_from_cmip6(gas, bundle_dir=bundle_dir).index.max()), "yr"
            ),
            "last year of CMIP6's global-, annual-mean",
        )
        add(
            f"{gas}-global-annual-mean-cmip6-last-year",
            partial(
                get_global_annual_mean, gas, CMIP6_LAST_YEAR, bundle_dir=bundle_dir
            ),
            f"our global-, annual-mean in {CMIP6_LAST_YEAR} "
            "(the last year of CMIP6's historical dataset)",
        )
        add(
            f"{gas}-global-annual-mean-last-year",
            partial(get_global_annual_mean, gas, bundle_dir=bundle_dir),
            "our global-, annual-mean in the last year of our output",
        )
        add(
            f"{gas}-obs-network-start",
            lambda gas=gas: Q(
                get_obs_network_years(gas, bundle_dir=bundle_dir)[0], "yr"
            ),
            "first year of the observation network's global-, annual-mean",
        )
        add(
            f"{gas}-obs-network-end",
            lambda gas=gas: Q(
                get_obs_network_years(gas, bundle_dir=bundle_dir)[1], "yr"
            ),
            "last year of the observation network's global-, annual-mean",
        )
        add(
            f"{gas}-diff-from-cmip6-obs-network",
            partial(
                get_max_abs_diff_from_cmip6_obs_network, gas, bundle_dir=bundle_dir
            ),
            "max abs difference from CMIP6 global-, annual-mean, "
            "from the start of the observation network",
        )

    add(
        "n2o-seasonality-reproduces-observed-average",
        partial(
            get_seasonality_diff_from_observed_average, "n2o", bundle_dir=bundle_dir
        ),
        "difference between our seasonality and the observed seasonality, "
        "both averaged over the observation network period "
        "(relative to the observed seasonality's magnitude)",
    )
    for i, (tag, description) in enumerate(
        (
            ("trudinger-start-year", "first"),
            ("trudinger-end-year", "last"),
        )
    ):
        add(
            tag,
            lambda i=i: Q(get_trudinger_years(bundle_dir=bundle_dir).m[i], "yr"),
            f"{description} year of the Trudinger et al. (2016) data",
        )
    add(
        "cfc12-like-non-zero-pre-industrial-gases",
        partial(has_expected_non_zero_pre_industrial_gases, bundle_dir=bundle_dir),
        "whether exactly CF4, CH2Cl2, CH3Br, CH3Cl and CHCl3 "
        "have a non-zero pre-industrial value",
    )
    add(
        "velders-first-year",
        lambda: Q(get_velders_first_year_values(bundle_dir=bundle_dir).name, "yr"),
        "first year of the Velders et al. (2022) data",
    )
    add(
        "velders-first-year-zero-except-hfc143a",
        partial(is_velders_first_year_zero_except, "hfc143a", bundle_dir=bundle_dir),
        "whether Velders et al. (2022)'s first-year value is zero "
        "for every HFC with a 1980 pre-industrial year except HFC-143a",
    )
    add(
        "velders-first-year-hfc143a",
        lambda: Q(
            float(get_velders_first_year_values(bundle_dir=bundle_dir)["hfc143a"]),
            "ppt",
        ),
        "Velders et al. (2022)'s first-year HFC-143a value",
    )
    for i, (tag, description) in enumerate(
        (("droste-first-year", "first"), ("droste-last-year", "last"))
    ):
        add(
            tag,
            lambda i=i: Q(get_droste_years(bundle_dir=bundle_dir)[i], "yr"),
            f"{description} year of the Droste et al. (2020) data",
        )
    add(
        "droste-max-first-year-value",
        partial(get_droste_max_first_year_value, bundle_dir=bundle_dir),
        "largest Droste et al. (2020) value in its first year (all gases and sites)",
    )
    for site in DROSTE_SITE_HEMISPHERES:
        add(
            f"droste-{site}-lat",
            partial(get_droste_site_lat, site, bundle_dir=bundle_dir),
            f"latitude of {site} in the Droste et al. (2020) data",
        )
    add(
        "c4f10-like-last-year",
        partial(get_c4f10_like_last_year, bundle_dir=bundle_dir),
        "last year of our output for the gases processed like C4F10",
    )
    add(
        "c4f10-like-erf-2022",
        partial(get_c4f10_like_erf, 2022, bundle_dir=bundle_dir),
        "approx. ERF of the gases processed like C4F10 in 2022",
    )
    add(
        "max-lat-gradient-fraction",
        lambda: Q(MAX_LAT_GRADIENT_FRACTION, "dimensionless"),
        "largest the latitudinal gradient's most negative value can be "
        "as a fraction of the global-mean",
    )
    add(
        "max-seasonality-fraction",
        lambda: Q(MAX_SEASONALITY_FRACTION, "dimensionless"),
        "largest the seasonality can be as a fraction of the global-mean",
    )
    for tag, getter, gases, description in (
        (
            "hfc152a-lat-gradient-capped-years",
            get_lat_gradient_capped_years,
            ("hfc152a",),
            "HFC-152a latitudinal gradient",
        ),
        (
            "hfc236fa-seasonality-capped-years",
            get_seasonality_capped_years,
            ("hfc236fa",),
            "HFC-236fa seasonality",
        ),
        (
            "c4f10-like-lat-gradient-capped-years",
            get_lat_gradient_capped_years,
            C4F10_LIKE_GASES,
            "latitudinal gradient of the gases processed like C4F10",
        ),
    ):
        for i, which in enumerate(("first", "last")):
            add(
                f"{tag}-{which}",
                lambda getter=getter, gases=gases, i=i: Q(
                    get_capped_years_extent(getter, gases, bundle_dir=bundle_dir).m[i],
                    "yr",
                ),
                f"{which} year in which the {description} was scaled down",
            )
    add(
        "cfc12-like-lat-gradient-capped-years-after-pre-industrial",
        partial(
            get_max_lat_gradient_capped_years_after_pre_industrial,
            ("hfc152a",),
            bundle_dir=bundle_dir,
        ),
        "largest number of years after the pre-industrial year in which "
        "the latitudinal gradient is scaled down (CFC12-like gases except HFC-152a)",
    )
    add(
        "cfc12-like-seasonality-capped-gases",
        partial(get_seasonality_capped_single_year_gases, bundle_dir=bundle_dir),
        "whether the seasonality is scaled down in every year for HFC-236fa, "
        "in one year for HFC-32, HFC-152a and HFC-365mfc and never for other gases",
    )
    add(
        "co2-seasonality-change-regression-start-year",
        partial(
            get_co2_seasonality_change_regression_start_year, bundle_dir=bundle_dir
        ),
        "first year in which the CO2 seasonality change PC comes from the regression",
    )
    add(
        "min-points-for-spatial-interpolation",
        lambda: Q(MIN_POINTS_FOR_SPATIAL_INTERPOLATION, "dimensionless"),
        "minimum number of binned values needed to spatially interpolate a month",
    )
    add(
        "harmonisation-transition-years",
        lambda: Q(HARMONISATION_TRANSITION_YEARS, "yr"),
        "number of years over which harmonisation offsets decline to zero",
    )

    def variance(decomposition: str, eofs: tuple[int, ...]) -> pint.Quantity:
        fractions = get_variance_explained_fraction(
            decomposition, bundle_dir=bundle_dir
        )
        return Q(float(fractions.loc[list(eofs)].sum()), "dimensionless")

    for gas in ("n2o", "co2", "ch4"):
        add(
            f"{gas}-lat-gradient-first-two-eofs-variance",
            partial(variance, f"{gas}_lat-gradient", (0, 1)),
            "fraction of the variance explained by the first two "
            "latitudinal gradient EOFs",
        )
    add(
        "co2-seasonality-change-first-eof-variance",
        partial(variance, "co2_seasonality-change", (0,)),
        "fraction of the variance explained by the first seasonality change EOF",
    )
    add(
        "co2-seasonality-change-max-other-eof-variance",
        lambda: Q(
            float(
                get_variance_explained_fraction(
                    "co2_seasonality-change", bundle_dir=bundle_dir
                )
                .iloc[1:]
                .max()
            ),
            "dimensionless",
        ),
        "largest fraction of the variance explained by any other "
        "seasonality change EOF",
    )
    for gas in ("cfc12", "hfc236fa"):
        add(
            f"{gas}-lat-gradient-first-eof-variance",
            partial(variance, f"{gas}_lat-gradient", (0,)),
            "fraction of the variance explained by the first latitudinal gradient EOF",
        )
    add(
        "cfc12-like-min-lat-gradient-first-eof-variance-excl-hfc236fa",
        partial(
            get_min_lat_gradient_first_eof_variance,
            ("cfc12", "hfc236fa"),
            bundle_dir=bundle_dir,
        ),
        "smallest fraction of the variance explained by the first latitudinal "
        "gradient EOF (CFC12-like gases except CFC12 and HFC-236fa)",
    )
    add(
        "n2o-menking-offset",
        partial(get_n2o_menking_offset, bundle_dir=bundle_dir),
        "N2O observation network global-, annual-mean minus Menking et al. (2025) "
        "in the harmonisation year",
    )
    add(
        "co2-mauna-loa-start-year",
        lambda: Q(CO2_MAUNA_LOA_START, "yr"),
        "first year of the Mauna Loa - Law Dome merged record we use",
    )
    add(
        "co2-mauna-loa-offset",
        partial(get_co2_mauna_loa_offset, bundle_dir=bundle_dir),
        "CO2 observation network global-, annual-mean minus the Mauna Loa merged "
        "record in the harmonisation year",
    )
    for i, (tag, description) in enumerate(
        (
            ("co2-menking-offset", "offset from Menking et al. (2025)"),
            (
                "co2-menking-output-match",
                "max abs difference from harmonised Menking et al. (2025)",
            ),
        )
    ):
        add(
            tag,
            lambda i=i: get_co2_menking_offset_and_match(bundle_dir=bundle_dir)[i],
            f"CO2 {description} in Menking et al. (2025)'s latitudinal bin",
        )
    for i, which in enumerate(("first", "last")):
        add(
            f"ch4-pc0-optimised-years-{which}",
            lambda i=i: Q(get_ch4_pc0_optimised_years(bundle_dir=bundle_dir)[i], "yr"),
            f"{which} year in which the CH4 first PC is optimised against ice cores",
        )
    for tag, getter, unit in (
        (
            "ch4-law-dome-noise-value-sd",
            lambda c: c["noise_adder"]["y_random_error"],
            None,
        ),
        (
            "ch4-law-dome-noise-time-relative-sd",
            lambda c: c["noise_adder"]["x_relative_random_error"],
            None,
        ),
        ("ch4-law-dome-noise-time-ref-year", lambda c: c["noise_adder"]["x_ref"], None),
        (
            "ch4-law-dome-min-points-either-side",
            lambda c: c["point_selector_settings"]["minimum_data_points_either_side"],
            "dimensionless",
        ),
        (
            "ch4-law-dome-max-points-either-side",
            lambda c: c["point_selector_settings"]["maximum_data_points_either_side"],
            "dimensionless",
        ),
        (
            "ch4-law-dome-window-width",
            lambda c: c["point_selector_settings"]["window_width"],
            None,
        ),
        ("ch4-law-dome-n-draws", lambda c: c["n_draws"], "dimensionless"),
    ):
        add(
            tag,
            lambda getter=getter, unit=unit: (
                Q(
                    getter(get_ch4_law_dome_smoothing_config(bundle_dir=bundle_dir)),
                    unit,
                )
                if unit is not None
                else Q(
                    *getter(get_ch4_law_dome_smoothing_config(bundle_dir=bundle_dir))
                )
            ),
            f"CH4 Law Dome smoothing setting ({tag})",
        )
    add(
        "ch4-law-dome-offset",
        partial(
            get_ch4_law_dome_offset, data_raw_dir=data_raw_dir, bundle_dir=bundle_dir
        ),
        "absolute offset between our CH4 output and smoothed Law Dome "
        "in the harmonisation year",
    )
    add(
        "ch4-neem-max-relative-diff",
        partial(
            get_ch4_max_neem_relative_diff,
            data_raw_dir=data_raw_dir,
            bundle_dir=bundle_dir,
        ),
        "largest difference between our CH4 output and NEEM, relative to NEEM",
    )
    add(
        "trudinger-harmonisation-transition-years",
        lambda: Q(HARMONISATION_TRANSITION_YEARS, "yr"),
        "number of years over which the Trudinger et al. (2016) offset declines",
    )
    for tag, source in (
        ("velders-pre-industrial-year", VELDERS_SOURCE),
        (
            "velders-adjusted-pre-industrial-year",
            "Velders et al., 2022 (with adjustments to support interpolation)",
        ),
    ):
        add(
            tag,
            partial(get_pre_industrial_years, source, bundle_dir=bundle_dir),
            f"pre-industrial year of the gases whose source is {source!r}",
        )
    add(
        "binning-grid-lon-band-width",
        partial(get_binning_grid_lon_band_width, bundle_dir=bundle_dir),
        "width of the longitudinal bands we bin the observations into",
    )
    add(
        "native-grid-lat-band-width",
        partial(get_native_grid_lat_band_width, bundle_dir=bundle_dir),
        "width of the latitudinal bands on our native grid",
    )
    add(
        "lat-gradient-eofs-area-weighted-mean",
        partial(get_max_lat_gradient_eof_area_weighted_mean, bundle_dir=bundle_dir),
        "largest area-weighted mean of any gas' latitudinal gradient EOFs "
        "(relative to the EOF's magnitude)",
    )

    add(
        "ch4-law-dome-start-year",
        partial(get_ch4_law_dome_start_year, data_raw_dir=data_raw_dir),
        "first year of the smoothed Law Dome CH4 data",
    )
    add(
        "ch4-first-year",
        lambda: Q(int(get_output("ch4", bundle_dir=bundle_dir)[1]["year"].min()), "yr"),
        "first year of our CH4 output",
    )

    for source in CH4_ICE_CORE_FILES:
        add(
            f"ch4-{source}-lat-bin",
            partial(
                get_ch4_ice_core_lat_bin,
                source,
                bundle_dir=bundle_dir,
                data_raw_dir=data_raw_dir,
            ),
            "centre of the native latitudinal bin "
            f"in which the {source} data is matched",
        )
        add(
            f"ch4-{source}-lat",
            partial(get_ch4_ice_core_lat, source, data_raw_dir=data_raw_dir),
            f"{source} latitude",
        )

    for gas in ("co2", "n2o"):
        add(
            f"{gas}-monthly-diff-from-noaa",
            lambda gas=gas: max_abs(
                get_diff_from_noaa_monthly(gas, bundle_dir=bundle_dir),
                get_units(gas, bundle_dir=bundle_dir),
            ),
            "max abs difference from NOAA global-mean monthly",
        )
        add(
            f"{gas}-diff-from-igcc",
            partial(get_max_abs_diff_from_igcc_noaa_era, gas, bundle_dir=bundle_dir),
            "max abs difference from IGCC global-, annual-mean, "
            "over the years NOAA's global-mean covers",
        )

    # CO2
    add(
        # Typo kept, to match the latex
        "co2-monthly-diff-from-maunoa-loa",
        lambda: max_abs(
            get_diff_from_mauna_loa(bundle_dir=bundle_dir),
            get_units("co2", bundle_dir=bundle_dir),
        ),
        "max abs difference from NOAA's Mauna Loa monthly record "
        "(our output interpolated to Mauna Loa's latitude)",
    )
    add(
        "co2-diff-from-cmip6-obs-network-1940-dip",
        partial(
            get_cmip6_feature_size, "co2", 1935, 1945, "dip", bundle_dir=bundle_dir
        ),
        "how far CMIP6 global-, annual-mean dips below ours, 1935 to 1945",
    )
    add(
        "co2-diff-from-cmip6-obs-network-1970-spike",
        partial(
            get_cmip6_feature_size, "co2", 1965, 1975, "spike", bundle_dir=bundle_dir
        ),
        "how far CMIP6 global-, annual-mean spikes above ours, 1965 to 1975",
    )

    # CH4
    add(
        "ch4-monthly-diff-from-noaa",
        lambda: Q(
            float(get_diff_from_noaa_monthly("ch4", bundle_dir=bundle_dir).mean()),
            get_units("ch4", bundle_dir=bundle_dir),
        ),
        "mean difference from NOAA global-mean monthly",
    )
    add(
        "ch4-diff-from-igcc",
        partial(get_mean_diff_from_igcc_noaa_era, "ch4", bundle_dir=bundle_dir),
        "mean difference from IGCC global-, annual-mean, "
        "over the years NOAA's global-mean covers",
    )
    add(
        "ch4-monthly-diff-from-uci",
        partial(get_mean_diff_from_uci_ch4, bundle_dir=bundle_dir),
        "mean difference from UCI global-mean quarterly",
    )

    # N2O
    add(
        "n2o-diff-from-cmip6-1750-1850-hump",
        partial(
            get_cmip6_feature_size, "n2o", 1750, 1850, "spike", bundle_dir=bundle_dir
        ),
        "how far CMIP6 global-, annual-mean rises above ours, 1750 to 1850",
    )

    # CFC-12
    for network in ("NOAA", "AGAGE"):
        add(
            f"cfc12-monthly-diff-from-{network.lower()}",
            partial(
                get_mean_diff_from_network, "cfc12", network, bundle_dir=bundle_dir
            ),
            f"mean difference from the {network} network's monthly observations "
            "(our output interpolated to each observation's latitude)",
        )

    add(
        "daniel-et-al-is-cfc12-global-override",
        partial(
            is_global_mean_from_source,
            "cfc12",
            "Daniel et al. (2022)",
            bundle_dir=bundle_dir,
        ),
        "whether Daniel et al. (2022) replaces the observation network's "
        "CFC-12 global-, annual-mean",
    )
    add(
        "cfc12-diff-from-igcc",
        lambda: max_abs(
            get_diff_from_igcc("cfc12", bundle_dir=bundle_dir),
            get_units("cfc12", bundle_dir=bundle_dir),
        ),
        "max abs difference from IGCC global-, annual-mean, all years both cover",
    )
    add(
        "cfc12-input-datasets-weakening-lat-gradient",
        partial(get_input_lat_gradient_weakening, "cfc12", bundle_dir=bundle_dir),
        "rate of weakening of the input data's latitudinal gradient "
        "over the last ten years",
    )

    # Equivalent species
    add_with_radiative_effect(
        "cfc12eq-diff-from-cmip6-1995-on",
        "cfc12eq",
        partial(get_max_abs_diff_from_cmip6, "cfc12eq", 1995, bundle_dir=bundle_dir),
        "max abs difference from CMIP6 global-, annual-mean, 1995 onwards",
    )
    add(
        "cfc12eq-diff-from-cmip6-radiative-efficiency-share",
        partial(
            get_difference_share,
            "cfc12eq",
            "CMIP6",
            "difference in radiative efficiency",
            "max",
            bundle_dir=bundle_dir,
        ),
        "share of the difference from CMIP6 due to radiative efficiencies, "
        "in the year of max abs difference",
    )
    add_with_radiative_effect(
        "cfc12eq-diff-from-cmip6-recalculated",
        "cfc12eq",
        partial(
            get_max_abs_diff_from_cmip6_recalculated, "cfc12eq", bundle_dir=bundle_dir
        ),
        "max abs difference from CMIP6 recalculated with AR6 radiative efficiencies",
    )
    add(
        "hfc134aeq-diff-from-cmip6-concentration-share",
        partial(
            get_difference_share,
            "hfc134aeq",
            "CMIP6",
            "other",
            "max",
            bundle_dir=bundle_dir,
        ),
        "share of the difference from CMIP6 due to concentrations, "
        "in the year of max abs difference",
    )
    for gas, part, share_tag in (
        ("cfc12eq", "difference in radiative efficiency", "radiative-efficiency"),
        ("hfc134aeq", "difference in included gases", "included-gases"),
    ):
        add_with_radiative_effect(
            f"{gas}-diff-from-igcc-last-year",
            gas,
            partial(get_diff_from_igcc_last_year, gas, bundle_dir=bundle_dir),
            "difference from IGCC global-, annual-mean in the last year both cover",
        )
        add(
            f"{gas}-diff-from-igcc-{share_tag}-share",
            partial(
                get_difference_share, gas, "IGCC", part, "last", bundle_dir=bundle_dir
            ),
            f"share of the difference from IGCC due to the {part}, "
            "in the last year both cover",
        )
        add(
            f"{gas}-diff-from-igcc-max-radiative-effect",
            partial(
                radiative_effect_of,
                gas,
                lambda gas=gas: max_abs(
                    get_diff_from_igcc(gas, bundle_dir=bundle_dir),
                    get_units(gas, bundle_dir=bundle_dir),
                ),
            ),
            "approx. radiative effect of the max abs difference from IGCC "
            "global-, annual-mean",
        )

    return tuple(res)
