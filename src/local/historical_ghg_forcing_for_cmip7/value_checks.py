"""
Checks of the values behind the statements in the historical manuscript

Each check calculates the value behind a statement,
and is matched by its tag to the `% value-check: {...}` comments in the latex.
See [local.value_checks][] for how that works.

Differences are always our (CMIP7) dataset minus the other dataset.
"""

from __future__ import annotations

from collections.abc import Callable
from functools import cache, partial
from typing import Literal

import numpy as np
import openscm_units
import pandas as pd
import pint
import xarray as xr

from local.cmip_ghg_generation import DEFAULT_BUNDLE_DIR
from local.historical_ghg_forcing_for_cmip7.cfc12_like_methods_figure import (
    get_cfc12_like_all_data_with_bins,
    get_global_mean_supplement,
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
from local.value_checks import CheckValue, ValueCheck

Q = openscm_units.unit_registry.Quantity

LAST_N_YEARS = 10
"""Number of years the statements about 'the last ten years' cover"""

NORTHERN_HEMISPHERE = "Northern hemisphere"
"""Name of the northern hemisphere in our hemispheric-mean output"""

SOUTHERN_HEMISPHERE = "Southern hemisphere"
"""Name of the southern hemisphere in our hemispheric-mean output"""


@cache
def get_output(
    gas: str,
) -> tuple[xr.DataArray, xr.DataArray, xr.DataArray, xr.DataArray]:
    """
    Get our output for a gas

    Parameters
    ----------
    gas
        Gas of interest

    Returns
    -------
    :
        Our native resolution output, yearly global-mean,
        monthly global-mean and monthly hemispheric-means
        (see [local.historical_ghg_forcing_for_cmip7.results_figure.load_output][])
    """
    return load_output(gas, DEFAULT_BUNDLE_DIR)


def get_units(gas: str) -> str:
    """
    Get the units of our output for a gas

    Parameters
    ----------
    gas
        Gas of interest

    Returns
    -------
    :
        Units of our output for `gas`
    """
    return str(get_output(gas)[1].attrs["units"])


def get_last_n_years(gas: str) -> np.ndarray:
    """
    Get the last [LAST_N_YEARS][] years of our output for a gas

    Parameters
    ----------
    gas
        Gas of interest

    Returns
    -------
    :
        The years
    """
    last_year = int(get_output(gas)[1]["year"].max())

    return np.arange(last_year - LAST_N_YEARS + 1, last_year + 1)


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


def get_last_n_years_trend(gas: str) -> pint.Quantity:
    """
    Get the linear trend in our global-, annual-mean over the last ten years

    Parameters
    ----------
    gas
        Gas of interest

    Returns
    -------
    :
        Trend
    """
    years = get_last_n_years(gas)
    values = get_output(gas)[1].sel(year=years).values
    gradient = np.polyfit(years, values, 1)[0]

    return Q(gradient, f"{get_units(gas)} / yr")


def get_last_n_years_second_derivative(gas: str) -> pint.Quantity:
    """
    Get the second derivative of our global-, annual-mean over the last ten years

    From a quadratic fit.
    Negative means the trend is getting more negative
    (e.g. a decline is accelerating).

    Parameters
    ----------
    gas
        Gas of interest

    Returns
    -------
    :
        Second derivative
    """
    years = get_last_n_years(gas)
    values = get_output(gas)[1].sel(year=years).values
    quadratic_coefficient = np.polyfit(years, values, 2)[0]

    return Q(2.0 * quadratic_coefficient, f"{get_units(gas)} / yr^2")


def get_last_n_years_seasonal_cycle(
    gas: str, regions: tuple[str, ...]
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

    Returns
    -------
    :
        Size of the seasonal cycle in each of `regions`
    """
    hm_monthly = get_output(gas)[3].sel(year=get_last_n_years(gas), lat=list(regions))
    amplitude = (hm_monthly.max("month") - hm_monthly.min("month")).mean("year")

    return Q(amplitude.values, get_units(gas))


def get_last_n_years_lat_gradient(gas: str) -> pint.Quantity:
    """
    Get the latitudinal gradient of our output over the last ten years

    The gradient of a linear fit to our native resolution output,
    averaged over the last ten years, against latitude.
    Positive means concentrations are higher in the north.

    Parameters
    ----------
    gas
        Gas of interest

    Returns
    -------
    :
        Latitudinal gradient
    """
    mean = get_output(gas)[0].sel(year=get_last_n_years(gas)).mean(["year", "month"])
    gradient = np.polyfit(mean["lat"].values, mean.values, 1)[0]

    return Q(gradient, f"{get_units(gas)} / degree")


@cache
def get_diff_from_cmip6(gas: str) -> pd.Series[float]:
    """
    Get the difference between our global-, annual-mean and CMIP6's

    Parameters
    ----------
    gas
        Gas of interest

    Returns
    -------
    :
        Difference, by year, in our units
    """
    (cmip6,) = get_cmip6_comparisons(gas, "yr")
    ours = get_output(gas)[1].to_series()

    return (ours - comparison_to_series(cmip6, get_units(gas), yearly=True)).dropna()


def get_max_abs_diff_from_cmip6(gas: str, start: int | None = None) -> pint.Quantity:
    """
    Get the maximum absolute difference from CMIP6

    Parameters
    ----------
    gas
        Gas of interest

    start
        First year to consider

        If `None`, all years.

    Returns
    -------
    :
        Maximum absolute difference from CMIP6
    """
    return max_abs(get_diff_from_cmip6(gas).loc[start:], get_units(gas))


def get_year_of_max_abs_diff_from_cmip6(gas: str) -> pint.Quantity:
    """
    Get the year in which the absolute difference from CMIP6 is greatest

    Parameters
    ----------
    gas
        Gas of interest

    Returns
    -------
    :
        Year of maximum absolute difference
    """
    return Q(int(get_diff_from_cmip6(gas).abs().idxmax()), "yr")


def get_cmip6_feature_size(
    gas: str, start: int, end: int, kind: Literal["dip", "spike"]
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

    Returns
    -------
    :
        How far CMIP6 dips below (or spikes above) our dataset
        between `start` and `end`.
        Negative if CMIP6 doesn't dip below (or spike above) our dataset at all.
    """
    diff = get_diff_from_cmip6(gas).loc[start:end]
    size = diff.max() if kind == "dip" else -diff.min()

    return Q(float(size), get_units(gas))


@cache
def get_obs_network_years(gas: str) -> tuple[int, int]:
    """
    Get the years covered by the observational network

    Parameters
    ----------
    gas
        Gas of interest

    Returns
    -------
    :
        First and last year of the observational network's global-, annual-mean
    """
    years = xr.load_dataset(
        DEFAULT_BUNDLE_DIR
        / "data"
        / "interim"
        / gas
        / f"{gas}_observational-network_global-annual-mean.nc"
    )["year"]

    return int(years.min()), int(years.max())


def get_max_abs_diff_from_cmip6_obs_network(gas: str) -> pint.Quantity:
    """
    Get the maximum absolute difference from CMIP6 over the observational network era

    Parameters
    ----------
    gas
        Gas of interest

    Returns
    -------
    :
        Maximum absolute difference from CMIP6
        from the first year of the observational network onwards
    """
    return get_max_abs_diff_from_cmip6(gas, start=get_obs_network_years(gas)[0])


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
def get_diff_from_noaa_monthly(gas: str) -> pd.Series[float]:
    """
    Get the difference between our global-mean monthly and NOAA's

    Parameters
    ----------
    gas
        Gas of interest

    Returns
    -------
    :
        Difference, for each month both cover, in our units
    """
    ours = monthly_to_series(get_output(gas)[2])
    noaa = comparison_to_series(get_noaa_global_mean(gas), get_units(gas), yearly=False)

    return (ours - noaa).dropna()


def get_diff_from_mauna_loa() -> pd.Series[float]:
    """
    Get the difference between our CO2 and NOAA's Mauna Loa record

    Our native resolution output is linearly interpolated to Mauna Loa's latitude.

    Returns
    -------
    :
        Difference, for each month both cover
    """
    gas = "co2"
    ours = monthly_to_series(get_output(gas)[0].interp(lat=MAUNA_LOA_LATITUDE))
    (mauna_loa,) = (
        c for c in get_noaa_comparisons(gas, deseasonalised=False) if c.is_spatial
    )

    return (
        ours - comparison_to_series(mauna_loa, get_units(gas), yearly=False)
    ).dropna()


def get_diff_from_igcc(gas: str, start: int | None = None) -> pd.Series[float]:
    """
    Get the difference between our global-, annual-mean and IGCC's

    Parameters
    ----------
    gas
        Gas of interest

    start
        First year to consider

        If `None`, all the years both cover.

    Returns
    -------
    :
        Difference, by year, in our units
    """
    ours = get_output(gas)[1].to_series()
    igcc = comparison_to_series(get_igcc_comparison(gas), get_units(gas), yearly=True)

    return (ours - igcc).dropna().loc[start:]


def get_max_abs_diff_from_igcc_noaa_era(gas: str) -> pint.Quantity:
    """
    Get the max absolute difference from IGCC over the years NOAA's global-mean covers

    IGCC's record is based on NOAA's over these years,
    which is what the statements about IGCC are about.

    Parameters
    ----------
    gas
        Gas of interest

    Returns
    -------
    :
        Maximum absolute difference from IGCC
    """
    noaa_start = int(np.floor(get_noaa_global_mean(gas).data[TIME_COLUMN].min()))

    return max_abs(get_diff_from_igcc(gas, start=noaa_start), get_units(gas))


def get_mean_diff_from_igcc_noaa_era(gas: str) -> pint.Quantity:
    """
    Get the mean difference from IGCC over the years NOAA's global-mean covers

    Parameters
    ----------
    gas
        Gas of interest

    Returns
    -------
    :
        Mean difference from IGCC
    """
    noaa_start = int(np.floor(get_noaa_global_mean(gas).data[TIME_COLUMN].min()))

    return Q(float(get_diff_from_igcc(gas, start=noaa_start).mean()), get_units(gas))


def get_mean_diff_from_uci_ch4() -> pint.Quantity:
    """
    Get the mean difference between our global-mean CH4 and UCI's

    UCI's record is quarterly, so our global-mean monthly
    is linearly interpolated to UCI's times.

    Returns
    -------
    :
        Mean difference
    """
    gas = "ch4"
    ours = monthly_to_series(get_output(gas)[2])
    uci = get_uci_ch4_comparison(deseasonalised=False).to_units(get_units(gas)).data
    in_range = uci[
        (uci[TIME_COLUMN] >= ours.index.min()) & (uci[TIME_COLUMN] <= ours.index.max())
    ]
    ours_at_uci_times = np.interp(in_range[TIME_COLUMN], ours.index, ours.to_numpy())

    return Q(float(np.mean(ours_at_uci_times - in_range[VALUE_COLUMN])), get_units(gas))


@cache
def get_cfc12_like_obs_network_diffs(gas: str) -> pd.DataFrame:
    """
    Get the difference between our output and each observation in the network

    Our native resolution output is linearly interpolated
    to each observation's latitude.

    Parameters
    ----------
    gas
        Gas of interest

    Returns
    -------
    :
        One row per observation, with its year, network and the difference
    """
    obs = get_cfc12_like_all_data_with_bins(gas, DEFAULT_BUNDLE_DIR)
    (obs_units,) = obs["unit"].unique()
    if obs_units != get_units(gas):
        raise NotImplementedError(obs_units)

    native = get_output(gas)[0]
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


def get_mean_diff_from_network(gas: str, network: str) -> pint.Quantity:
    """
    Get the mean difference between our output and a network's observations

    Parameters
    ----------
    gas
        Gas of interest

    network
        Network of interest

    Returns
    -------
    :
        Mean difference over all the network's observations
    """
    diffs = get_cfc12_like_obs_network_diffs(gas)

    return Q(
        float(diffs.loc[diffs["network"] == network, "diff"].mean()), get_units(gas)
    )


def is_global_mean_from_source(gas: str, source: str) -> bool:
    """
    Get whether a gas' global-, annual-mean comes from a source, rather than the network

    Parameters
    ----------
    gas
        Gas of interest

    source
        Label of the source

    Returns
    -------
    :
        Whether `source` replaces the observational network's global-, annual-mean
    """
    supplement = get_global_mean_supplement(gas, DEFAULT_BUNDLE_DIR)
    if supplement is None:
        return False

    label, supplement_data = supplement
    max_year_extended = int(
        xr.load_dataset(
            interim_dir(gas, DEFAULT_BUNDLE_DIR)
            / f"{gas}_global-annual-mean_allyears.nc"
        )["year"].max()
    )

    return label == source and supplement_replaces_obs_network(
        supplement_data, max_year_extended
    )


def get_input_lat_gradient_weakening(gas: str) -> pint.Quantity:
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

    Returns
    -------
    :
        Rate of weakening of the latitudinal gradient
    """
    obs = get_cfc12_like_all_data_with_bins(gas, DEFAULT_BUNDLE_DIR)
    years = get_last_n_years(gas)
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

    return Q(-trend, f"{get_units(gas)} / degree / yr")


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
    equivalent_species: str, other: Literal["CMIP6", "IGCC"]
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

    Returns
    -------
    :
        Decomposition of the difference, by year
    """
    cmip7 = load_cmip7(DEFAULT_BUNDLE_DIR)
    other_dataset = load_cmip6() if other == "CMIP6" else load_igcc()
    for dataset in (cmip7, other_dataset):
        check_definition_reproduces(dataset, equivalent_species)

    return decompose_difference(equivalent_species, cmip7, other_dataset)


def get_difference_share(
    equivalent_species: str,
    other: Literal["CMIP6", "IGCC"],
    part: str,
    year: Literal["max", "last"],
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

    Returns
    -------
    :
        Share of the difference due to `part`
    """
    decomposition = get_equivalent_species_decomposition(equivalent_species, other)
    if year == "max":
        selected_year = decomposition["total difference"].abs().idxmax()
    else:
        selected_year = decomposition.index.max()

    share = (
        decomposition.loc[selected_year, part]
        / decomposition.loc[selected_year, "total difference"]
    )

    return Q(float(share), "dimensionless")


def get_max_abs_diff_from_cmip6_recalculated(equivalent_species: str) -> pint.Quantity:
    """
    Get the max abs difference from CMIP6 recalculated with AR6 radiative efficiencies

    I.e. the difference which is left
    once the difference in radiative efficiencies is removed.

    Parameters
    ----------
    equivalent_species
        Equivalent species of interest

    Returns
    -------
    :
        Maximum absolute difference
    """
    decomposition = get_equivalent_species_decomposition(equivalent_species, "CMIP6")
    recalculated = (
        decomposition["total difference"]
        - decomposition["difference in radiative efficiency"]
    )

    return max_abs(recalculated, get_units(equivalent_species))


def get_diff_from_igcc_last_year(gas: str) -> pint.Quantity:
    """
    Get the difference from IGCC in the last year both datasets cover

    Parameters
    ----------
    gas
        Gas of interest

    Returns
    -------
    :
        Difference from IGCC
    """
    return Q(float(get_diff_from_igcc(gas).iloc[-1]), get_units(gas))


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


def get_value_checks() -> tuple[ValueCheck, ...]:
    """
    Get the value checks for the historical manuscript

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
            partial(get_last_n_years_trend, gas),
            "linear trend in the global-, annual-mean over the last ten years",
        )
        add(
            f"{gas}-last-10-years-second-derivative",
            partial(get_last_n_years_second_derivative, gas),
            "second derivative of the global-, annual-mean over the last ten years",
        )
        if gas in seasonal_cycle_regions:
            add(
                f"{gas}-last-10-years-lat-grad",
                partial(get_last_n_years_lat_gradient, gas),
                "latitudinal gradient (linear fit) over the last ten years",
            )
            regions = seasonal_cycle_regions[gas]
            add(
                f"{gas}-last-10-years-seasonal-cycle",
                partial(get_last_n_years_seasonal_cycle, gas, regions),
                "mean peak-to-trough seasonal cycle over the last ten years "
                f"({', '.join(r.lower() for r in regions)})",
            )

        add_with_radiative_effect(
            f"{gas}-diff-from-cmip6-all",
            gas,
            partial(get_max_abs_diff_from_cmip6, gas),
            "max abs difference from CMIP6 global-, annual-mean, all years",
        )
        add(
            f"{gas}-diff-from-cmip6-all-year",
            partial(get_year_of_max_abs_diff_from_cmip6, gas),
            "year of max abs difference from CMIP6 global-, annual-mean",
        )
        add_with_radiative_effect(
            f"{gas}-diff-from-cmip6-1850-on",
            gas,
            partial(get_max_abs_diff_from_cmip6, gas, 1850),
            "max abs difference from CMIP6 global-, annual-mean, 1850 onwards",
        )

    for gas in ("co2", "ch4", "n2o", "cfc12"):
        add(
            f"{gas}-obs-network-start",
            lambda gas=gas: Q(get_obs_network_years(gas)[0], "yr"),
            "first year of the observational network's global-, annual-mean",
        )
        add(
            f"{gas}-obs-network-end",
            lambda gas=gas: Q(get_obs_network_years(gas)[1], "yr"),
            "last year of the observational network's global-, annual-mean",
        )
        add(
            f"{gas}-diff-from-cmip6-obs-network",
            partial(get_max_abs_diff_from_cmip6_obs_network, gas),
            "max abs difference from CMIP6 global-, annual-mean, "
            "from the start of the observational network",
        )

    for gas in ("co2", "n2o"):
        add(
            f"{gas}-monthly-diff-from-noaa",
            lambda gas=gas: max_abs(get_diff_from_noaa_monthly(gas), get_units(gas)),
            "max abs difference from NOAA global-mean monthly",
        )
        add(
            f"{gas}-diff-from-igcc",
            partial(get_max_abs_diff_from_igcc_noaa_era, gas),
            "max abs difference from IGCC global-, annual-mean, "
            "over the years NOAA's global-mean covers",
        )

    # CO2
    add(
        # Typo kept, to match the latex
        "co2-monthly-diff-from-maunoa-loa",
        lambda: max_abs(get_diff_from_mauna_loa(), get_units("co2")),
        "max abs difference from NOAA's Mauna Loa monthly record "
        "(our output interpolated to Mauna Loa's latitude)",
    )
    add(
        "co2-diff-from-cmip6-obs-network-1940-dip",
        partial(get_cmip6_feature_size, "co2", 1935, 1945, "dip"),
        "how far CMIP6 global-, annual-mean dips below ours, 1935 to 1945",
    )
    add(
        "co2-diff-from-cmip6-obs-network-1970-spike",
        partial(get_cmip6_feature_size, "co2", 1965, 1975, "spike"),
        "how far CMIP6 global-, annual-mean spikes above ours, 1965 to 1975",
    )

    # CH4
    add(
        "ch4-monthly-diff-from-noaa",
        lambda: Q(float(get_diff_from_noaa_monthly("ch4").mean()), get_units("ch4")),
        "mean difference from NOAA global-mean monthly",
    )
    add(
        "ch4-diff-from-igcc",
        partial(get_mean_diff_from_igcc_noaa_era, "ch4"),
        "mean difference from IGCC global-, annual-mean, "
        "over the years NOAA's global-mean covers",
    )
    add(
        "ch4-monthly-diff-from-uci",
        get_mean_diff_from_uci_ch4,
        "mean difference from UCI global-mean quarterly",
    )

    # N2O
    add(
        "n2o-diff-from-cmip6-1750-1850-hump",
        partial(get_cmip6_feature_size, "n2o", 1750, 1850, "spike"),
        "how far CMIP6 global-, annual-mean rises above ours, 1750 to 1850",
    )

    # CFC-12
    for network in ("NOAA", "AGAGE"):
        add(
            f"cfc12-monthly-diff-from-{network.lower()}",
            partial(get_mean_diff_from_network, "cfc12", network),
            f"mean difference from the {network} network's monthly observations "
            "(our output interpolated to each observation's latitude)",
        )

    add(
        "daniel-et-al-is-cfc12-global-override",
        partial(is_global_mean_from_source, "cfc12", "Daniel et al. (2022)"),
        "whether Daniel et al. (2022) replaces the observational network's "
        "CFC-12 global-, annual-mean",
    )
    add(
        "cfc12-diff-from-igcc",
        lambda: max_abs(get_diff_from_igcc("cfc12"), get_units("cfc12")),
        "max abs difference from IGCC global-, annual-mean, all years both cover",
    )
    add(
        "cfc12-input-datasets-weakening-lat-gradient",
        partial(get_input_lat_gradient_weakening, "cfc12"),
        "rate of weakening of the input data's latitudinal gradient "
        "over the last ten years",
    )

    # Equivalent species
    add_with_radiative_effect(
        "cfc12eq-diff-from-cmip6-1995-on",
        "cfc12eq",
        partial(get_max_abs_diff_from_cmip6, "cfc12eq", 1995),
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
        ),
        "share of the difference from CMIP6 due to radiative efficiencies, "
        "in the year of max abs difference",
    )
    add_with_radiative_effect(
        "cfc12eq-diff-from-cmip6-recalculated",
        "cfc12eq",
        partial(get_max_abs_diff_from_cmip6_recalculated, "cfc12eq"),
        "max abs difference from CMIP6 recalculated with AR6 radiative efficiencies",
    )
    add(
        "hfc134aeq-diff-from-cmip6-concentration-share",
        partial(get_difference_share, "hfc134aeq", "CMIP6", "other", "max"),
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
            partial(get_diff_from_igcc_last_year, gas),
            "difference from IGCC global-, annual-mean in the last year both cover",
        )
        add(
            f"{gas}-diff-from-igcc-{share_tag}-share",
            partial(get_difference_share, gas, "IGCC", part, "last"),
            f"share of the difference from IGCC due to the {part}, "
            "in the last year both cover",
        )
        add(
            f"{gas}-diff-from-igcc-max-radiative-effect",
            partial(
                radiative_effect_of,
                gas,
                lambda gas=gas: max_abs(get_diff_from_igcc(gas), get_units(gas)),
            ),
            "approx. radiative effect of the max abs difference from IGCC "
            "global-, annual-mean",
        )

    return tuple(res)


VALUE_CHECKS = get_value_checks()
"""Value checks for the historical manuscript"""
