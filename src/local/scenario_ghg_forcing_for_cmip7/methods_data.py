"""
Loading of the data behind the scenario manuscript's methods figures

Everything is returned as yearly data, indexed by decimal year
at the middle of each year (the same way the published yearly data is placed,
see [local.scenario_ghg_forcing_for_cmip7.esgf_data.to_decimal_year][]),
so it can all be drawn on the same time axes.
"""

from __future__ import annotations

import json
from functools import cache
from pathlib import Path

import numpy as np
import openscm_units
import pandas as pd
import xarray as xr

from local.cmip_ghg_generation import DEFAULT_BUNDLE_DIR as HISTORICAL_BUNDLE_DIR
from local.historical_ghg_forcing_for_cmip7.comparison_data import (
    get_radiative_effect_per_unit,
)
from local.scenario_ghg_forcing_for_cmip7.bundle import BUNDLE_DIR
from local.scenario_ghg_forcing_for_cmip7.esgf_data import (
    convert_units,
    to_decimal_year,
)
from local.scenario_ghg_forcing_for_cmip7.magicc_output import (
    MAGICC_GAS_NAMES,
    MAGICC_VERSION,
    get_concentration_variable,
    get_magicc_quantile,
)
from local.scenario_ghg_forcing_for_cmip7.rerun_projections import (
    CO2_SEASONALITY_PC_FILE,
    CO2_SEASONALITY_REGRESSION_KEY,
    COMPLETE_EMISSIONS_DIR,
    MARKERS,
    ONE_EOF_GASES,
    SINGLE_PROJECTION_GASES,
    TWO_EOF_GASES,
    get_pcs_file,
    get_regression_file,
    get_regression_info_file,
)
from local.scenario_ghg_forcing_for_cmip7.scenarios import (
    CMIP7_HISTORY_SOURCE_ID,
    SCENARIOS,
)

LAST_HISTORICAL_YEAR = 2022
"""Last year of the historical dataset (the scenarios start from it)"""

MID_YEAR = 0.5
"""Offset from the start of a year to where its yearly value is drawn"""

HISTORICAL_INTERIM_DIR = HISTORICAL_BUNDLE_DIR / "data" / "interim"
"""The historical run's interim data"""

INTERIM_DIR = Path("data") / "interim"
"""The scenario run's interim data, relative to the bundle's root"""

RERUN_GASES = (*TWO_EOF_GASES, *ONE_EOF_GASES, *SINGLE_PROJECTION_GASES)
"""Gases whose latitudinal gradients we re-ran

The archived ones were overwritten by the VL-CF run.
"""


def to_decimal_years(values: pd.Series) -> pd.Series:
    """Index yearly values by decimal year, at the middle of each year"""
    res = values.copy()
    res.index = res.index.astype(float) + MID_YEAR
    res.index.name = "time"

    return res


def convert_series(values: pd.Series, units_from: str, units_to: str) -> pd.Series:
    """Convert values from one unit to another"""
    if units_from == units_to:
        return values

    factor = openscm_units.unit_registry.Quantity(1.0, units_from).to(units_to).m

    return values * factor


def da_to_series(da: xr.DataArray) -> pd.Series:
    """Get yearly data with a decimal year `time` as a series"""
    return pd.Series(da.to_numpy(), index=pd.Index(da["time"].to_numpy(), name="time"))


def get_history(gas: str, units: str) -> pd.Series:
    """
    Get the historical dataset's yearly global-mean

    Parameters
    ----------
    gas
        Gas of interest

    units
        Units to return the data in

    Returns
    -------
    :
        Yearly global-mean from 1750 up to and including [LAST_HISTORICAL_YEAR][]
    """
    da = convert_units(load_historical_yearly_global_mean(gas), units)

    return da_to_series(da.sel(time=da["time"] < LAST_HISTORICAL_YEAR + 1))


HISTORICAL_ESGF_READY_DIR = (
    HISTORICAL_BUNDLE_DIR
    / "data"
    / "processed"
    / "esgf-ready"
    / "input4MIPs"
    / "CMIP7"
    / "CMIP"
    / "CR"
    / CMIP7_HISTORY_SOURCE_ID
)
"""The published historical files, as archived in the historical run's bundle"""


@cache
def load_historical_yearly_global_mean(gas: str) -> xr.DataArray:
    """
    Load the historical dataset's published yearly global-mean (1750 to 2022)

    From the historical run's bundle's copy of the published files
    (Zenodo record 14892947, `data--processed.tar.gz`),
    i.e. the same files as are on ESGF, but without the hours of downloading.

    Parameters
    ----------
    gas
        Gas of interest

    Returns
    -------
    :
        Yearly global-mean, with `time` as decimal years
    """
    (file,) = (HISTORICAL_ESGF_READY_DIR / "atmos" / "yr" / gas / "gm").glob(
        "*/*_gm_1750-2022.nc"
    )
    da = xr.load_dataset(file, decode_times=xr.coders.CFDatetimeCoder(use_cftime=True))[
        gas
    ]

    return to_decimal_year(da, "yr")


ESGF_READY_DIR = (
    Path("data")
    / "processed"
    / "esgf-ready"
    / "input4MIPs"
    / "CMIP7"
    / "ScenarioMIP"
    / "CR"
)
"""The published files, as archived in the bundle, relative to the bundle's root"""


@cache
def load_published_yearly_global_mean(
    gas: str, scenario: str, bundle_dir: Path = BUNDLE_DIR
) -> xr.DataArray:
    """
    Load a scenario's published yearly global-mean (2022 to 2500)

    From the bundle's copy of the published files,
    i.e. the same files as are on ESGF
    (see [local.scenario_ghg_forcing_for_cmip7.bundle.BUNDLE_TARBALLS][] for why).

    Parameters
    ----------
    gas
        Gas of interest

    scenario
        Scenario of interest

    bundle_dir
        Directory which holds the bundle

    Returns
    -------
    :
        Yearly global-mean, with `time` as decimal years
    """
    pieces = []
    for source_id in SCENARIOS[scenario].source_ids:
        files = sorted(
            (
                bundle_dir / ESGF_READY_DIR / source_id / "atmos" / "yr" / gas / "gm"
            ).glob("*/*.nc")
        )
        if not files:
            msg = f"No published files for {gas=} {source_id=} in {bundle_dir}"
            raise FileNotFoundError(msg)

        for file in files:
            da = xr.load_dataset(
                file, decode_times=xr.coders.CFDatetimeCoder(use_cftime=True)
            )[gas]
            pieces.append(to_decimal_year(da, "yr"))

    res = xr.concat(pieces, dim="time").sortby("time")
    if not (np.diff(res["time"].values) > 0).all():
        msg = f"Published files for {gas=} {scenario=} overlap"
        raise AssertionError(msg)

    return res


def get_scenario(
    gas: str, scenario: str, units: str, bundle_dir: Path = BUNDLE_DIR
) -> pd.Series:
    """
    Get a scenario's published yearly global-mean (2022 to 2500)

    Parameters
    ----------
    gas
        Gas of interest

    scenario
        Scenario of interest

    units
        Units to return the data in

    bundle_dir
        Directory which holds the bundle

    Returns
    -------
    :
        Yearly global-mean
    """
    return da_to_series(
        convert_units(
            load_published_yearly_global_mean(gas, scenario, bundle_dir), units
        )
    )


def get_units(gas: str, bundle_dir: Path = BUNDLE_DIR) -> str:
    """Get the units the gas' published data is in"""
    return load_published_yearly_global_mean(gas, "l", bundle_dir).attrs["units"]


@cache
def load_magicc_summary(summary_file: Path) -> pd.DataFrame:
    """Load the MAGICC summary (cached, as it is used by every figure)"""
    return pd.read_feather(summary_file)


def get_magicc(  # noqa: PLR0913
    summary_file: Path,
    variable: str,
    scenario: str,
    quantile: float,
    units: str | None = None,
    climate_model: str = MAGICC_VERSION,
) -> pd.Series:
    """
    Get one quantile of one MAGICC variable in one scenario

    Parameters
    ----------
    summary_file
        File which holds the MAGICC summary

    variable
        MAGICC variable

    scenario
        Scenario of interest

    quantile
        Quantile of interest

    units
        Units to return the data in

        If `None`, MAGICC's units are kept.

    climate_model
        MAGICC version

    Returns
    -------
    :
        Values, indexed by decimal year, with the units as the series' name
    """
    res = get_magicc_quantile(
        load_magicc_summary(summary_file),
        variable,
        scenario,
        quantile,
        climate_model=climate_model,
    )
    magicc_units = str(res.name)
    if units is not None:
        res = convert_series(res, magicc_units, units)
        res.name = units

    return to_decimal_years(res)


def get_magicc_concentrations(  # noqa: PLR0913
    summary_file: Path,
    gas: str,
    scenario: str,
    quantile: float,
    units: str,
    climate_model: str = MAGICC_VERSION,
) -> pd.Series:
    """Get one quantile of MAGICC's concentrations of a gas in one scenario"""
    return get_magicc(
        summary_file,
        get_concentration_variable(gas),
        scenario,
        quantile,
        units=units,
        climate_model=climate_model,
    )


def get_historical_pcs(
    gas: str, kind: str = "lat-gradient", bundle_dir: Path = BUNDLE_DIR
) -> pd.DataFrame | None:
    """
    Get the historical dataset's PCs

    Parameters
    ----------
    gas
        Gas of interest

    kind
        "lat-gradient" or "seasonality-change"

    bundle_dir
        Directory which holds the bundle (only needed for C8F18)

    Returns
    -------
    :
        PCs, indexed by decimal year, one column per EOF,
        up to [LAST_HISTORICAL_YEAR][].
        `None` if the historical dataset has no PCs for the gas.
    """
    if gas == C8F18 and kind == "lat-gradient":
        return get_c8f18_historical_pc(bundle_dir)

    file = HISTORICAL_INTERIM_DIR / gas / f"{gas}_allyears-{kind}-eofs-pcs.nc"
    if not file.exists():
        return None

    pcs = xr.load_dataset(file)["principal-components"]
    pcs = pcs.sel(year=pcs["year"] <= LAST_HISTORICAL_YEAR)
    res = pcs.to_pandas().T
    res.columns = res.columns.astype(int)

    return to_decimal_years(res)


def get_projected_pcs(gas: str, bundle_dir: Path = BUNDLE_DIR) -> pd.DataFrame:
    """
    Get the projected latitudinal gradient PCs

    Parameters
    ----------
    gas
        Gas of interest

    bundle_dir
        Directory which holds the bundle (including our re-runs)

    Returns
    -------
    :
        Annual PCs from [LAST_HISTORICAL_YEAR][], indexed by decimal year,
        with columns `(eof, scenario)`.
        Gases whose projection doesn't depend on the scenario have one scenario, "all".
    """
    if gas in RERUN_GASES:
        pcs = pd.read_csv(bundle_dir / get_pcs_file(gas))
        res = pcs.pivot_table(index="year", columns=["eof", "scenario"], values="value")

    else:
        # Monthly PCs, which we take the annual-mean of
        # (the monthly values are a mean-preserving interpolation of annual values)
        info = xr.open_dataset(
            bundle_dir
            / INTERIM_DIR
            / "latitudinal-gradient"
            / f"{gas}_latitudinal-gradient-info.nc",
            decode_times=xr.coders.CFDatetimeCoder(use_cftime=True),
        )
        monthly = info["principal-components-monthly"]
        if "eof" not in monthly.dims:
            monthly = monthly.expand_dims(eof=info["eof"].values)

        annual = monthly.groupby("time.year").mean()
        annual = annual.sel(year=annual["year"] >= LAST_HISTORICAL_YEAR)
        res = annual.to_dataframe(name="value")["value"].unstack(["eof", "scenario"])

    res.index = res.index.astype(int)

    return to_decimal_years(res.sort_index(axis="columns"))


C8F18 = "c8f18"
"""C8F18, whose latitudinal gradient is handled differently to every other gas"""


def get_hemispheric_mean_difference(lat_gradient: xr.DataArray) -> xr.DataArray:
    """
    Get the difference between the northern and southern hemispheric means

    Parameters
    ----------
    lat_gradient
        Latitudinal gradient, with a `lat` dimension (15 degree bands)

    Returns
    -------
    :
        Northern minus southern hemispheric (area-weighted) mean
    """
    weights = np.cos(np.deg2rad(lat_gradient["lat"]))
    north = lat_gradient["lat"] > 0.0
    northern = lat_gradient.where(north).weighted(weights).mean("lat")
    southern = lat_gradient.where(~north).weighted(weights).mean("lat")

    return northern - southern


def get_c8f18_historical_pc(bundle_dir: Path = BUNDLE_DIR) -> pd.DataFrame:
    """
    Get the historical dataset's C8F18 latitudinal gradient PC

    In the historical dataset, C8F18's latitudinal gradient has an assumed shape
    (linear in latitude), scaled to reproduce the CMIP6 hemispheric means,
    so it is fully described by the difference between its hemispheric means.
    The original run never wrote out a PC for it.

    We express it relative to the EOF the scenario run uses,
    i.e. the December 2022 latitudinal gradient (notebook `1030`),
    so the PC is dimensionless (the units sit on the EOF)
    and lines up with the projected PC (which is one, see [get_projected_pcs][]).
    Both gradients have the same (linear) shape,
    so the PC is the ratio of the differences between their hemispheric means.

    Parameters
    ----------
    bundle_dir
        Directory which holds the bundle

    Returns
    -------
    :
        PC, indexed by decimal year, one column (EOF 0), up to [LAST_HISTORICAL_YEAR][]
    """
    ds = xr.load_dataset(
        HISTORICAL_INTERIM_DIR / C8F18 / f"{C8F18}_hemispheric-mean_annual-mean.nc"
    )
    (hemispheric_mean,) = ds.data_vars.values()
    difference = hemispheric_mean.sel(lat=45.0) - hemispheric_mean.sel(lat=-45.0)
    difference = difference.sel(year=difference["year"] <= LAST_HISTORICAL_YEAR)
    difference_units = hemispheric_mean.attrs["units"]

    info = xr.open_dataset(
        bundle_dir
        / INTERIM_DIR
        / "latitudinal-gradient"
        / f"{C8F18}_latitudinal-gradient-info.nc",
        decode_times=xr.coders.CFDatetimeCoder(use_cftime=True),
    )
    # The EOF is in the units of the published data
    eof_difference = float(
        get_hemispheric_mean_difference(info["eofs"].sel(eof=0)).values
    )
    eof_difference = convert_series(
        pd.Series([eof_difference]), get_units(C8F18), difference_units
    ).iloc[0]

    res = pd.DataFrame(
        {0: difference.to_numpy() / eof_difference},
        index=difference["year"].to_numpy().astype(int),
    )

    return to_decimal_years(res)


def get_projected_co2_seasonality_change_pc(
    bundle_dir: Path = BUNDLE_DIR,
) -> pd.DataFrame:
    """
    Get the projected CO2 seasonality change PC

    Parameters
    ----------
    bundle_dir
        Directory which holds the bundle (including our re-runs)

    Returns
    -------
    :
        Annual PC from [LAST_HISTORICAL_YEAR][], indexed by decimal year,
        one column per scenario
    """
    pcs = pd.read_csv(bundle_dir / CO2_SEASONALITY_PC_FILE)
    res = pcs.pivot_table(index="year", columns="scenario", values="value")
    res.index = res.index.astype(int)

    return to_decimal_years(res)


def get_regression(
    key: str, bundle_dir: Path = BUNDLE_DIR
) -> tuple[pd.DataFrame, dict[str, str | float]]:
    """
    Get the data behind, and result of, a regression we re-derived

    Parameters
    ----------
    key
        Gas (for the latitudinal gradient regressions)
        or [CO2_SEASONALITY_REGRESSION_KEY][]

    bundle_dir
        Directory which holds the bundle (including our re-runs)

    Returns
    -------
    :
        Data (columns `year`, `x` and `y`)
        and result (units, slope and intercept) of the regression
    """
    data = pd.read_csv(bundle_dir / get_regression_file(key))
    info = json.loads((bundle_dir / get_regression_info_file(key)).read_text())

    return data, info


CO2_SEASONALITY_REGRESSION = CO2_SEASONALITY_REGRESSION_KEY
"""Key of the CO2 seasonality change PC's regression, see [get_regression][]"""


def get_single_projection_inverse_emissions(
    gas: str, units: str, bundle_dir: Path = BUNDLE_DIR
) -> pd.Series:
    """
    Get the inverse emissions of a gas which uses a single concentration projection

    Parameters
    ----------
    gas
        Gas of interest

    units
        Units to return the emissions in (must be mass per time)

    bundle_dir
        Directory which holds the bundle

    Returns
    -------
    :
        Emissions, indexed by decimal year
    """
    raw = pd.read_feather(
        bundle_dir
        / INTERIM_DIR
        / "inverse-emissions"
        / f"single-concentration-projection_{gas}_inverse-emissions.feather"
    )
    mass = raw[~raw.index.get_level_values("unit").str.contains("ppt")]
    if mass.shape[0] != 1:
        raise AssertionError(mass)

    res = convert_series(
        mass.iloc[0], mass.index.get_level_values("unit")[0], units
    ).astype(float)
    res.index = res.index.astype(int)

    return to_decimal_years(res)


def get_official_emissions(
    gas: str, scenario: str, units: str, bundle_dir: Path = BUNDLE_DIR
) -> pd.Series:
    """
    Get a gas' emissions in a scenario, as used to run MAGICC

    These are the official emissions, including the official extensions beyond 2100.

    Parameters
    ----------
    gas
        Gas of interest

    scenario
        Scenario of interest

    units
        Units to return the emissions in

    bundle_dir
        Directory which holds the bundle

    Returns
    -------
    :
        Emissions, indexed by decimal year
    """
    ((model, iam_scenario),) = ((m, s) for m, s, name in MARKERS if name == scenario)
    stem = f"{iam_scenario}_{model}".replace(" ", "_").replace(".", "-")
    emissions = load_complete_emissions(
        bundle_dir / COMPLETE_EMISSIONS_DIR / f"{stem}.feather"
    )
    variable = f"Emissions|{MAGICC_GAS_NAMES[gas]}"
    selected = emissions[emissions.index.get_level_values("variable") == variable]
    if selected.shape[0] != 1:
        raise AssertionError(selected)

    res = convert_series(
        selected.iloc[0], selected.index.get_level_values("unit")[0], units
    ).astype(float)
    res.index = res.index.astype(int)

    return to_decimal_years(res)


@cache
def load_complete_emissions(file: Path) -> pd.DataFrame:
    """Load a scenario's complete emissions (cached, they are used repeatedly)"""
    return pd.read_feather(file)


def get_western_et_al_2024(
    gas: str, bundle_dir: Path = BUNDLE_DIR
) -> tuple[pd.Series, pd.Series]:
    """
    Get the Western et al. (2024) projection of an HCFC, raw and as we extended it

    Parameters
    ----------
    gas
        Gas of interest

    bundle_dir
        Directory which holds the bundle

    Returns
    -------
    :
        Raw projection (to 2100) and extended projection (to 2500),
        both indexed by decimal year, in ppt
    """
    western_dir = bundle_dir / INTERIM_DIR / "western-et-al-2024"
    raw = pd.read_feather(western_dir / "cleaned-mixing-ratios.feather")
    raw = raw.xs(gas, level="ghg")
    extended = pd.read_feather(
        western_dir / f"western-et-al-2024_{gas}_extended.feather"
    )
    extended = extended.xs(gas, level="ghg")
    for df in (raw, extended):
        if df.index.tolist() != ["ppt"]:
            raise AssertionError(df.index)

    return tuple(
        to_decimal_years(pd.Series(df.iloc[0].to_numpy(), index=df.columns.astype(int)))
        for df in (raw, extended)
    )


def get_wmo_2022(gas: str, bundle_dir: Path = BUNDLE_DIR) -> pd.Series:
    """
    Get the WMO (2022) projection of a gas, as cleaned by the original run

    Parameters
    ----------
    gas
        Gas of interest

    bundle_dir
        Directory which holds the bundle

    Returns
    -------
    :
        Projection, indexed by decimal year, in ppt
    """
    wmo = pd.read_feather(
        bundle_dir / INTERIM_DIR / "wmo-2022" / "cleaned-mixing-ratios.feather"
    )
    wmo = wmo.xs(gas, level="ghg")
    if wmo.index.tolist() != ["ppt"]:
        raise AssertionError(wmo.index)

    return to_decimal_years(
        pd.Series(wmo.iloc[0].to_numpy(), index=wmo.columns.astype(int))
    )


def get_max_radiative_effect_difference(
    summary_file: Path, gas: str, scenario: str
) -> float:
    """
    Get the largest radiative effect difference between MAGICC's and our concentrations

    For a single-projection gas.
    The largest over time of the absolute difference in approximate radiative effect
    between MAGICC's median and our concentrations.

    Parameters
    ----------
    summary_file
        File which holds the MAGICC summary

    gas
        Gas of interest

    scenario
        Scenario of interest

    Returns
    -------
    :
        Largest absolute difference, from 2022 to 2500, in W / m^2
    """
    units = get_units(gas)
    ours = get_scenario(gas, scenario, units)
    magicc = get_magicc_concentrations(summary_file, gas, scenario, 0.5, units)
    common = ours.index.intersection(magicc.index)
    diff = (magicc[common] - ours[common]) * get_radiative_effect_per_unit(gas, units)

    return float(np.abs(diff).max())
