# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.17.2
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Run MAGICC with and without the CMIP7 history
#
# Written by the CMIP7 ScenarioMIP GHG manuscript's build
# (`src/local/scenario_ghg_forcing_for_cmip7/magicc_jump.py`), not part of the original run.
#
# Built from the original run's `0021_run-magicc` (which produced the published MAGICC output)
# and `1200_run-magicc-concentration-driven` (which writes CMIP7 concentrations as MAGICC input files).
#
# The published runs use the CMIP6 concentrations (and, for the species MAGICC runs emissions-driven
# throughout, the CMIP6 emissions) up to 2014, then the CMIP7 emissions from 2015.
# This notebook runs one experiment in one of two set-ups:
#
# - `published`: exactly as the original run, to check we can reproduce the published output
# - `cmip7-history`: the same, except that MAGICC gets the CMIP7 historical concentrations
#   (by writing `*_CONC.IN` files) and the CMIP7 historical emissions
#   (by starting the scenario, i.e. the SCEN7 file, in 1750 rather than 2015),
#   so there is no jump from CMIP6 to CMIP7 inputs in 2014/2015
#
# The difference between the two is the largest effect we could expect from the jump.

# %% [markdown]
# ## Imports

# %%
import json
import os
import re
from functools import partial
from pathlib import Path
from typing import Any

import pandas as pd
import pandas_indexing as pix
import pandas_openscm
import tqdm.auto
import xarray as xr
from gcages.renaming import SupportedNamingConventions, convert_variable_name
from gcages.scm_running import convert_openscm_runner_output_names_to_magicc_output_names, run_scms
from pymagicc.definitions import convert_magicc7_to_openscm_variables
from pymagicc.io import MAGICCData

# %% [markdown]
# ## Parameters

# %% tags=["parameters"]
model: str = "MESSAGEix-GLOBIOM-GAINS 2.1-M-R12"
scenario: str = "SSP2 - Low Emissions"
complete_file: str = (
    "data/interim/complete-emissions/SSP2_-_Low_Emissions_MESSAGEix-GLOBIOM-GAINS_2-1-M-R12.feather"
)
setup: str = "published"
n_runs: int = 0
magicc_exe: str = "magicc/magicc-v7.6.0a3/bin/magicc-darwin-arm64"
magicc_prob_distribution: str = "magicc/magicc-v7.6.0a3/configs/magicc-ar7-fast-track-drawnset-v0-3-0.json"
n_magicc_workers: int = 4
historical_concentrations_root: str = "../../cmip-ghg-concentration-generation/v1.0.0/data/processed/esgf-ready"
out_file: str = "manuscript-outputs/magicc-jump/l_published.feather"

# %% [markdown]
# ## Parse parameters

# %%
SETUPS = ("published", "cmip7-history")
if setup not in SETUPS:
    raise ValueError(setup)

complete_file_p = Path(complete_file)
magicc_exe_p = Path(magicc_exe)
magicc_prob_distribution_p = Path(magicc_prob_distribution)
historical_concentrations_root_p = Path(historical_concentrations_root)
out_file_p = Path(out_file)
magicc_run_dir = magicc_exe_p.parents[1] / "run"

# %% [markdown]
# ## Set up

# %%
pandas_openscm.register_pandas_accessor()

# %% [markdown]
# ## Load the complete emissions (1750 to 2500)

# %%
complete = pd.read_feather(complete_file_p)
complete_openscm_runner = complete.openscm.update_index_levels(
    {
        "variable": partial(
            convert_variable_name,
            from_convention=SupportedNamingConventions.GCAGES,
            to_convention=SupportedNamingConventions.OPENSCM_RUNNER,
        )
    },
)

# %% [markdown]
# ## Configure
#
# The published runs pass MAGICC the scenario from 2015 only,
# so MAGICC uses its own historical (CMIP6, RCMIP v4.0.0) emissions files up to 2014.
# Passing the scenario from 1750 instead makes MAGICC use the CMIP7 historical emissions.
# The emissions outputs below let us check that it does.

# %%
MAGICC_START_SCENARIO_YEAR = {"published": 2015, "cmip7-history": 1750}[setup]
complete_openscm_runner_for_magicc = complete_openscm_runner.loc[:, MAGICC_START_SCENARIO_YEAR:]

# %%
output_variables = (
    # As in the published run
    "Surface Air Temperature Change",
    "Effective Radiative Forcing",
    "Effective Radiative Forcing|Anthropogenic",
    "Effective Radiative Forcing|Aerosols",
    "Effective Radiative Forcing|Greenhouse Gases",
    "Effective Radiative Forcing|CO2",
    "Effective Radiative Forcing|CH4",
    "Effective Radiative Forcing|N2O",
    "Atmospheric Concentrations|CO2",
    "Atmospheric Concentrations|CH4",
    "Atmospheric Concentrations|N2O",
    # To check which historical emissions MAGICC used
    "Emissions|BC|MAGICC Fossil and Industrial",
    "Emissions|NOx|MAGICC Fossil and Industrial",
    "Emissions|SOx|MAGICC Fossil and Industrial",
    "Emissions|CO2|MAGICC Fossil and Industrial",
)

# %% [markdown]
# ## Write the CMIP7 historical concentrations as MAGICC input files
#
# Only for the `cmip7-history` set-up.
# Every gas keeps switching from concentrations to emissions in 2015 (MAGICC's default),
# so only the concentrations up to 2014 are used.
# Gases which aren't in the CMIP7 history (Halon-1202) keep MAGICC's default file.

# %%
CONC_MAGICC_FLAG_MAP = {
    "Atmospheric Concentrations|CO2": "file_co2_conc",
    "Atmospheric Concentrations|CH4": "file_ch4_conc",
    "Atmospheric Concentrations|N2O": "file_n2o_conc",
    "Atmospheric Concentrations|CF4": "fgas_files_conc__0",
    "Atmospheric Concentrations|C2F6": "fgas_files_conc__1",
    "Atmospheric Concentrations|C3F8": "fgas_files_conc__2",
    "Atmospheric Concentrations|C4F10": "fgas_files_conc__3",
    "Atmospheric Concentrations|C5F12": "fgas_files_conc__4",
    "Atmospheric Concentrations|C6F14": "fgas_files_conc__5",
    "Atmospheric Concentrations|C7F16": "fgas_files_conc__6",
    "Atmospheric Concentrations|C8F18": "fgas_files_conc__7",
    "Atmospheric Concentrations|cC4F8": "fgas_files_conc__8",
    "Atmospheric Concentrations|HFC23": "fgas_files_conc__9",
    "Atmospheric Concentrations|HFC32": "fgas_files_conc__10",
    "Atmospheric Concentrations|HFC4310": "fgas_files_conc__11",
    "Atmospheric Concentrations|HFC125": "fgas_files_conc__12",
    "Atmospheric Concentrations|HFC134a": "fgas_files_conc__13",
    "Atmospheric Concentrations|HFC143a": "fgas_files_conc__14",
    "Atmospheric Concentrations|HFC152a": "fgas_files_conc__15",
    "Atmospheric Concentrations|HFC227ea": "fgas_files_conc__16",
    "Atmospheric Concentrations|HFC236fa": "fgas_files_conc__17",
    "Atmospheric Concentrations|HFC245fa": "fgas_files_conc__18",
    "Atmospheric Concentrations|HFC365mfc": "fgas_files_conc__19",
    "Atmospheric Concentrations|NF3": "fgas_files_conc__20",
    "Atmospheric Concentrations|SF6": "fgas_files_conc__21",
    "Atmospheric Concentrations|SO2F2": "fgas_files_conc__22",
    "Atmospheric Concentrations|CFC11": "mhalo_files_conc__0",
    "Atmospheric Concentrations|CFC12": "mhalo_files_conc__1",
    "Atmospheric Concentrations|CFC113": "mhalo_files_conc__2",
    "Atmospheric Concentrations|CFC114": "mhalo_files_conc__3",
    "Atmospheric Concentrations|CFC115": "mhalo_files_conc__4",
    "Atmospheric Concentrations|HCFC22": "mhalo_files_conc__5",
    "Atmospheric Concentrations|HCFC141b": "mhalo_files_conc__6",
    "Atmospheric Concentrations|HCFC142b": "mhalo_files_conc__7",
    "Atmospheric Concentrations|CH3CCl3": "mhalo_files_conc__8",
    "Atmospheric Concentrations|CCl4": "mhalo_files_conc__9",
    "Atmospheric Concentrations|CH3Cl": "mhalo_files_conc__10",
    "Atmospheric Concentrations|CH2Cl2": "mhalo_files_conc__11",
    "Atmospheric Concentrations|CHCl3": "mhalo_files_conc__12",
    "Atmospheric Concentrations|CH3Br": "mhalo_files_conc__13",
    "Atmospheric Concentrations|Halon1211": "mhalo_files_conc__14",
    "Atmospheric Concentrations|Halon1301": "mhalo_files_conc__15",
    "Atmospheric Concentrations|Halon2402": "mhalo_files_conc__16",
    "Atmospheric Concentrations|Halon1202": "mhalo_files_conc__17",
}


# %%
def read_default_file_list(key: str) -> list[str]:
    """
    Read a list of files from MAGICC's default configuration

    Parameters
    ----------
    key
        Key of the list, e.g. `fgas_files_conc`

    Returns
    -------
    :
        The default list
    """
    cfg = (magicc_run_dir / "MAGCFG_DEFAULTALL.CFG").read_text()
    (line,) = (line for line in cfg.splitlines() if line.strip().upper().startswith(f"{key.upper()} "))

    return re.findall(r"'([^']*)'", line)


# %%
def to_df(da: xr.DataArray, variable: str) -> pd.DataFrame:
    """
    Convert yearly global-means to the form MAGICCData expects
    """
    return (
        da.groupby("time.year")
        .mean()
        .to_dataframe()
        .T.rename_axis("variable")
        .pix.assign(variable=variable, scenario="historical", model="CMIP7", unit=da.attrs["units"], region="World")
    )


# %%
magicc_concentration_cfg: dict[str, Any] = {}
if setup == "cmip7-history":
    magicc_concentration_cfg = {
        "file_co2_conc": None,
        "file_ch4_conc": None,
        "file_n2o_conc": None,
        "fgas_files_conc": read_default_file_list("fgas_files_conc"),
        "mhalo_files_conc": read_default_file_list("mhalo_files_conc"),
    }
    magicc_file_dir = magicc_run_dir / "cmip7-history"
    magicc_file_dir.mkdir(exist_ok=True, parents=True)

    # The yearly files are split into chunks of time (e.g. 0001-0999, 1000-1749, 1750-2022)
    files_by_gas: dict[str, list[Path]] = {}
    for fp in sorted(historical_concentrations_root_p.rglob("**/yr/*/gm/**/*.nc")):
        files_by_gas.setdefault(fp.name.split("_")[0], []).append(fp)

    for ghg, fps in tqdm.auto.tqdm(sorted(files_by_gas.items())):
        if ghg.endswith("eq"):
            # Not needed by MAGICC
            continue

        da = xr.concat(
            [xr.open_dataset(fp, use_cftime=True)[ghg] for fp in fps],
            dim="time",
        ).sortby("time")
        da = da.sel(time=da["time"].dt.year >= 1700)

        ghg_magicc = ghg.replace("hfc4310mee", "hfc4310")
        openscm_runner_variable = convert_magicc7_to_openscm_variables(f"{ghg_magicc}_conc".upper())

        writer = MAGICCData(to_df(da, variable=openscm_runner_variable))
        writer["todo"] = "SET"
        writer.metadata = {"header": f"CMIP7 historical {ghg}", "source": fp.name}

        fn = magicc_file_dir / f"CMIP7_HISTORICAL_{ghg_magicc}_CONC.IN".upper()
        writer.write(str(fn), magicc_version=7)

        magicc_flag = CONC_MAGICC_FLAG_MAP[openscm_runner_variable]
        relative = str(fn.relative_to(magicc_run_dir))
        if "__" in magicc_flag:
            key, index = magicc_flag.split("__")
            magicc_concentration_cfg[key][int(index)] = relative
        else:
            magicc_concentration_cfg[magicc_flag] = relative

    missing = [k for k, v in magicc_concentration_cfg.items() if v is None]
    if missing:
        raise AssertionError(missing)

magicc_concentration_cfg


# %%
def load_magicc_cfgs(
    prob_distribution_path: Path,
    output_variables: tuple[str, ...],
    magicc_concentration_cfg: dict[str, Any],
    startyear: int = 1750,
) -> dict[str, list[dict[str, Any]]]:
    """
    Load MAGICC's configuration

    As in the original run's `0021_run-magicc` and `1200_run-magicc-concentration-driven`.
    """
    with open(prob_distribution_path) as fh:
        cfgs_raw = json.load(fh)

    cfgs_physical = [
        {
            "run_id": c["paraset_id"],
            **{k.lower(): v for k, v in c["nml_allcfgs"].items()},
        }
        for c in cfgs_raw["configurations"]
    ]

    common_cfg = {
        "startyear": startyear,
        "out_dynamic_vars": convert_openscm_runner_output_names_to_magicc_output_names(output_variables),
        "out_ascii_binary": "BINARY",
        "out_binary_format": 2,
    }

    run_config = [{**common_cfg, **magicc_concentration_cfg, **physical_cfg} for physical_cfg in cfgs_physical]

    return {"MAGICC7": run_config}


# %%
climate_models_cfgs = load_magicc_cfgs(
    prob_distribution_path=magicc_prob_distribution_p,
    output_variables=output_variables,
    magicc_concentration_cfg=magicc_concentration_cfg,
    startyear=1750,
)
if n_runs > 0:
    climate_models_cfgs["MAGICC7"] = climate_models_cfgs["MAGICC7"][:n_runs]

len(climate_models_cfgs["MAGICC7"])

# %%
badly_converted = [v for v in climate_models_cfgs["MAGICC7"][0]["out_dynamic_vars"] if v.upper() != v]
if badly_converted:
    raise AssertionError(badly_converted)

# %% [markdown]
# ## Run

# %%
os.environ["MAGICC_EXECUTABLE_7"] = str(magicc_exe_p)

# %%
# papermill_description=run-magicc
res = run_scms(
    scenarios=complete_openscm_runner_for_magicc,
    climate_models_cfgs=climate_models_cfgs,
    output_variables=output_variables,
    scenario_group_levels=["model", "scenario"],
    n_processes=n_magicc_workers,
    db=None,
    verbose=True,
    progress=True,
)
res = res.pix.assign(setup=setup)

# %% [markdown]
# ## Save

# %%
out_file_p.parent.mkdir(exist_ok=True, parents=True)
to_save = res.reset_index()
to_save.columns = to_save.columns.astype(str)
to_save.to_feather(out_file_p)
out_file_p
