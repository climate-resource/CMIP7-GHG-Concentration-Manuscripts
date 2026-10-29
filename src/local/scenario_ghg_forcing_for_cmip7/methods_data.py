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
from local.scenario_ghg_forcing_for_cmip7.rerun_vl_cf import (
    CO2_SEASONALITY_VL_CF_PC_FILE,
    VL_CF_MARKER,
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


ESGF_READY_DIRS = tuple(
    Path("data") / "processed" / "esgf-ready" / "input4MIPs" / "CMIP7" / activity / "CR"
    for activity in ("ScenarioMIP", "PolMIP")
)
"""The published files, as archived in the bundle, relative to the bundle's root

The markers are in ScenarioMIP, VL-CF in PolMIP.
"""


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
            file
            for esgf_ready_dir in ESGF_READY_DIRS
            for file in (
                bundle_dir / esgf_ready_dir / source_id / "atmos" / "yr" / gas / "gm"
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


def get_projected_pcs(
    gas: str, bundle_dir: Path = BUNDLE_DIR, include_vl_cf: bool = False
) -> pd.DataFrame:
    """
    Get the projected latitudinal gradient PCs

    Parameters
    ----------
    gas
        Gas of interest

    bundle_dir
        Directory which holds the bundle (including our re-runs)

    include_vl_cf
        Whether to include VL-CF
        (ignored for the gases whose projection doesn't depend on the scenario,
        as "all" already covers it)

    Returns
    -------
    :
        Annual PCs from [LAST_HISTORICAL_YEAR][], indexed by decimal year,
        with columns `(eof, scenario)`.
        Gases whose projection doesn't depend on the scenario have one scenario, "all".
        VL-CF's PCs start where it leaves history (2015) for the gases which use MAGICC.
    """
    if gas in RERUN_GASES:
        pcs = pd.read_csv(bundle_dir / get_pcs_file(gas))
        res = pcs.pivot_table(index="year", columns=["eof", "scenario"], values="value")
        if include_vl_cf and gas not in SINGLE_PROJECTION_GASES:
            res = res.join(get_vl_cf_pcs(gas, bundle_dir), how="outer")

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
        if include_vl_cf:
            # The VL-CF run didn't write these gases' latitudinal gradients,
            # because, for these gases, VL-CF is VL from 2022 onwards
            # (checked against the published native resolution data, 2022-2500).
            vl_cf = res.loc[:, pd.IndexSlice[:, "vl"]].rename(
                columns={"vl": VL_CF_MARKER[2]}, level="scenario"
            )
            res = res.join(vl_cf)

    res.index = res.index.astype(int)

    return to_decimal_years(res.sort_index(axis="columns"))


def get_vl_cf_pcs(gas: str, bundle_dir: Path = BUNDLE_DIR) -> pd.DataFrame:
    """
    Get VL-CF's projected latitudinal gradient PCs

    The original run's interim files for the gases in `TWO_EOF_GASES`
    and `ONE_EOF_GASES` only hold VL-CF (see `rerun_projections`),
    so we take them from there.

    Parameters
    ----------
    gas
        Gas of interest, one of `TWO_EOF_GASES` or `ONE_EOF_GASES`

    bundle_dir
        Directory which holds the bundle

    Returns
    -------
    :
        Annual PCs from the start of VL-CF's projection (2015), indexed by year,
        with columns `(eof, scenario)`
    """
    if gas not in (*TWO_EOF_GASES, *ONE_EOF_GASES):
        raise NotImplementedError(gas)

    info = xr.open_dataset(
        bundle_dir
        / INTERIM_DIR
        / "latitudinal-gradient"
        / f"{gas}_latitudinal-gradient-info.nc",
        decode_times=xr.coders.CFDatetimeCoder(use_cftime=True),
    )
    monthly = info["principal-components-monthly"]
    if list(monthly["scenario"].values) != [VL_CF_MARKER[2]]:
        msg = (
            f"Expected only VL-CF in the interim file, got {monthly['scenario'].values}"
        )
        raise AssertionError(msg)

    if "eof" not in monthly.dims:
        monthly = monthly.expand_dims(eof=info["eof"].values)

    # The monthly values are a mean-preserving interpolation of annual values
    annual = monthly.groupby("time.year").mean()

    return annual.to_dataframe(name="value")["value"].unstack(["eof", "scenario"])


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
    pcs = pd.concat(
        [
            pd.read_csv(bundle_dir / CO2_SEASONALITY_PC_FILE),
            pd.read_csv(bundle_dir / CO2_SEASONALITY_VL_CF_PC_FILE),
        ]
    )
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
    emissions = get_complete_emissions(scenario, bundle_dir)
    variable = f"Emissions|{MAGICC_GAS_NAMES[gas]}"
    selected = emissions[emissions.index.get_level_values("variable") == variable]
    if selected.shape[0] != 1:
        raise AssertionError(selected)

    res = convert_series(
        selected.iloc[0], selected.index.get_level_values("unit")[0], units
    ).astype(float)
    res.index = res.index.astype(int)

    return to_decimal_years(res)


def get_complete_emissions(
    scenario: str, bundle_dir: Path = BUNDLE_DIR
) -> pd.DataFrame:
    """
    Get a scenario's complete emissions, as used to run MAGICC

    Parameters
    ----------
    scenario
        Scenario of interest

    bundle_dir
        Directory which holds the bundle

    Returns
    -------
    :
        Every species' emissions, from 1750 to 2500
        (history up to the harmonisation year, then the scenario
        and its extension beyond 2100)
    """
    ((model, iam_scenario),) = (
        (m, s) for m, s, name in (*MARKERS, VL_CF_MARKER) if name == scenario
    )
    stem = f"{iam_scenario}_{model}".replace(" ", "_").replace(".", "-")

    return load_complete_emissions(
        bundle_dir / COMPLETE_EMISSIONS_DIR / f"{stem}.feather"
    )


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


WESTERN_ET_AL_2024_LIFETIMES = {"hcfc22": 11.6, "hcfc141b": 8.81, "hcfc142b": 17.1}
"""Lifetimes (yr) with which the original run extended Western et al. (2024)

Copied from notebook `0001_extend-western-et-al-2024`.
"""


def get_ppt_to_mass_factor(
    gas: str, units: str, bundle_dir: Path = BUNDLE_DIR
) -> float:
    """
    Get the factor which converts a gas' emissions from ppt / yr to mass per time

    Taken from the original run's inverse emissions (notebook `1020`),
    which are saved in both ppt / yr and mass per time,
    so the conversion is exactly the one the original run used.

    Parameters
    ----------
    gas
        Gas of interest

    units
        Mass per time units to convert to

    bundle_dir
        Directory which holds the bundle

    Returns
    -------
    :
        Factor by which to multiply emissions in ppt / yr
    """
    raw = pd.read_feather(
        bundle_dir
        / INTERIM_DIR
        / "inverse-emissions"
        / f"single-concentration-projection_{gas}_inverse-emissions.feather"
    )
    is_ppt = raw.index.get_level_values("unit").str.contains("ppt")
    ppt, mass = raw[is_ppt], raw[~is_ppt]
    if ppt.shape[0] != 1 or mass.shape[0] != 1:
        raise AssertionError(raw)

    mass_converted = convert_series(
        mass.iloc[0], mass.index.get_level_values("unit")[0], units
    ).astype(float)
    non_zero = ppt.iloc[0] != 0.0
    factors = mass_converted[non_zero] / ppt.iloc[0][non_zero].astype(float)
    np.testing.assert_allclose(factors, factors.iloc[0], rtol=1e-10)

    return float(factors.iloc[0])


def get_western_et_al_2024_emissions(
    gas: str, units: str, bundle_dir: Path = BUNDLE_DIR
) -> tuple[pd.Series, pd.Series]:
    """
    Get the emissions behind our extension of the Western et al. (2024) projection

    The original run (notebook `0001`) assumes the gas follows a one-box model,
    infers the emissions in 2100 from the concentration and its rate of change,
    then holds them constant and runs the one-box model forward to 2500.
    We apply the same inversion to every year of the projection,
    so the emissions which lead up to 2100 can be shown too,
    and check that running forward from 2100 reproduces our extension.

    Parameters
    ----------
    gas
        Gas of interest, a key of [WESTERN_ET_AL_2024_LIFETIMES][]

    units
        Mass per time units in which to return the emissions

    bundle_dir
        Directory which holds the bundle

    Returns
    -------
    :
        Emissions implied by the Western et al. (2024) projection (to 2100)
        and the constant emissions which drive our extension (2100 onwards),
        both indexed by decimal year
    """
    raw, extended = get_western_et_al_2024(gas, bundle_dir)
    tau = WESTERN_ET_AL_2024_LIFETIMES[gas]
    # Emissions in concentration units, alpha E,
    # with annual timesteps, as in notebook `0001`
    alpha_e = (raw.diff() + raw / tau).dropna()

    alpha_e_constant = alpha_e.iloc[-1]
    rerun = [raw.iloc[-1]]
    extension_years = extended.index[extended.index > raw.index[-1]]
    for _ in extension_years:
        rerun.append(rerun[-1] + alpha_e_constant - rerun[-1] / tau)

    np.testing.assert_allclose(rerun[1:], extended.loc[extension_years], rtol=1e-10)

    factor = get_ppt_to_mass_factor(gas, units, bundle_dir)
    constant = pd.Series(
        alpha_e_constant * factor,
        index=extended.index[extended.index >= raw.index[-1]],
    )

    return alpha_e * factor, constant


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
