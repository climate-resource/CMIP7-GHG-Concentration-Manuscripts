"""
Re-running the original run's seasonality and latitudinal gradient projections

In the archived bundle (v1.1.0), the interim files for the gases
which use gradient-aware harmonisation (CO2, CH4, N2O and some HFCs)
only hold VL-CF: the later VL-CF run overwrote them.
This includes the emissions which drive the latitudinal gradients,
the latitudinal gradient PCs and the CO2 seasonality.
The published (ESGF) files are fine,
it is only the interim files which were overwritten.

To get these interim results for the other scenarios,
we re-run the notebooks which make them, for all the marker scenarios,
with the same inputs as the original run.
We re-run with all the markers, not just the scenario we show,
because the monthly interpolation of the PCs fixes the value at the end of 2022
to the mean over all the scenarios which are run together.
Re-running with all the markers reproduces the published values.

We also save the regressions and the annual PCs,
which the original run only plots.

Everything is written below [RERUN_DIR][] in the bundle,
so the original run's data is never overwritten.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import pandas as pd

from local.cmip_ghg_generation import (
    DEFAULT_BUNDLE_DIR as HISTORICAL_BUNDLE_DIR,
)
from local.cmip_ghg_generation import (
    ensure_bundle_available,
    ensure_bundle_environment,
    run_in_bundle,
    run_notebook_from_bundle_dir,
    write_modified_notebook,
)
from local.scenario_ghg_forcing_for_cmip7.bundle import (
    BUNDLE_DIR,
    EXECUTED_NOTEBOOKS_DIR,
    MODIFIED_NOTEBOOKS_DIR,
    RERUN_NOTEBOOKS_DIR,
)
from local.scenario_ghg_forcing_for_cmip7.scenarios import CMIP7_HISTORY_SOURCE_ID

MARKERS = (
    ("REMIND-MAgPIE 3.5-4.11", "SSP1 - Very Low Emissions", "vl"),
    ("AIM 3.0", "SSP2 - Low Overshoot_a", "ln"),
    ("MESSAGEix-GLOBIOM-GAINS 2.1-M-R12", "SSP2 - Low Emissions", "l"),
    ("COFFEE 1.6", "SSP2 - Medium-Low Emissions", "ml"),
    ("IMAGE 3.4", "SSP2 - Medium Emissions", "m"),
    ("WITCH 6.0", "SSP5 - Medium-Low Emissions_a", "hl"),
    ("GCAM 8s", "SSP3 - High Emissions", "h"),
)
"""Marker scenarios as (model, scenario, CMIP scenario name)

From `markers` in the original run's `scripts/generate-concentration-files.py`.
"""

HARMONISATION_YEAR = 2023
"""Year in which the scenarios' emissions are harmonised to history

The original run's default (`scripts/generate-concentration-files.py`).
"""

EMISSIONS_BATCH_ID = "202603251220_202512071232_202511040855_202511040855"
"""ID of the emissions the original run used

From the name of the bundle's `data/interim/input-emissions` directory.
"""

RAW_FOSSIL_BIO_SPLIT_FILE = (
    Path("data")
    / "raw"
    / "input-scenarios"
    / f"{EMISSIONS_BATCH_ID}_harmonised-emissions-fossil-biosphere-aggregation.csv"
)
"""Emissions split into fossil and biosphere sources, relative to the bundle's root"""

COMPLETE_EMISSIONS_DIR = Path("data") / "interim" / "complete-emissions"
"""Complete emissions of each scenario, relative to the bundle's root"""

MAGICC_OUTPUT_DB_DIR = Path("data") / "interim" / "magicc-output" / "db"
"""MAGICC output database, relative to the bundle's root"""

RERUN_DIR = Path("manuscript-outputs") / "rerun-markers"
"""Where the re-runs write everything, relative to the bundle's root"""

FOSSIL_BIO_SPLIT_FILE = RERUN_DIR / "fossil-bio-split-markers.feather"
"""Fossil/biosphere split for the markers, relative to the bundle's root"""

SINGLE_VARIABLE_DIR = RERUN_DIR / "single-variable-files"
"""Emissions which drive the latitudinal gradients, relative to the bundle's root"""

HISTORICAL_CONCENTRATIONS_ROOT = (
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
"""Root of the historical concentrations, as published

The historical run's bundle's copy of the published files
(Zenodo record 14892947, `data--processed.tar.gz`),
i.e. the same files as are on ESGF, without querying ESGF.
The notebooks find the files they need in here with globs.
"""

HISTORICAL_INTERIM_ROOT = HISTORICAL_BUNDLE_DIR / "data" / "interim"
"""Root of the historical run's interim data (seasonality and latitudinal gradient info)

The same data the original run used (Zenodo record 14892947, `data--interim.tar.gz`).
"""

EOF_SCALING_VARIABLES = {
    "co2": "emissions|co2|fossil",
    "ch4": "emissions|ch4|fossil",
    "n2o": "emissions|n2o",
}
"""Emissions which drive the first latitudinal gradient PC of the gases with two EOFs"""

TWO_EOF_GASES = tuple(EOF_SCALING_VARIABLES)
"""Gases whose latitudinal gradient is projected with notebook `1031`"""

ONE_EOF_GASES = ("hfc134a", "hfc152a", "hfc23", "hfc245fa", "hfc32", "hfc365mfc")
"""Gases processed like N2O whose interim files were overwritten (notebook `1030`)"""

SINGLE_PROJECTION_GASES = (
    "ccl4",
    "cfc11",
    "cfc113",
    "cfc114",
    "cfc115",
    "cfc12",
    "ch3br",
    "ch3ccl3",
    "ch3cl",
    "halon1211",
    "halon1301",
    "halon2402",
    "hcfc141b",
    "hcfc142b",
    "hcfc22",
)
"""Gases which use a single concentration projection (in our output)

Their latitudinal gradient files were overwritten by the VL-CF run too
(which started them in 2015, rather than 2022),
so we re-run them with notebook `1030`.
As in the original run, they are driven by their inverse emissions
(see `single_concentration_projection_flow.py` in the original run),
which were not overwritten.
"""

SEASONALITY_NOTEBOOK = EXECUTED_NOTEBOOKS_DIR / (
    "1011_scale-seasonality-based-on-magicc-npp"
    "_modelling-based-projection_co2_seasonality-all-time.ipynb"
)
"""Notebook which projects the CO2 seasonality, as executed (for VL-CF)

Only its code matters, as we replace its parameters.
"""

TWO_EOF_LAT_GRADIENT_NOTEBOOK = EXECUTED_NOTEBOOKS_DIR / (
    "1031_scale-latitudinal-gradient-eofs_co2_latitudinal-gradient-info.ipynb"
)
"""Notebook which projects the two-EOF latitudinal gradients, as executed (for VL-CF)"""

ONE_EOF_LAT_GRADIENT_NOTEBOOK = EXECUTED_NOTEBOOKS_DIR / (
    "1030_scale-latitudinal-gradient-based-on-emissions"
    "_hfc125_latitudinal-gradient-info.ipynb"
)
"""Notebook which projects the one-EOF latitudinal gradients, as executed

As executed for HFC-125 in the main run (March 2026), which wasn't overwritten.
We start from it,
rather than one of the notebooks executed in the VL-CF run (July 2026),
because the VL-CF run changed this notebook:
it added a `harmonisation_year` parameter (default 2015)
and starts the projection from the earlier of that and the end of the history.
The main run's version always starts the projection from the end of the history (2022),
which is what the published data comes from.
(Notebooks `1011` and `1031` were not changed in a way which matters:
re-running the VL-CF versions with the main run's parameters
reproduces the published data.)
"""


def get_eof_scaling_emissions_file(ghg: str) -> Path:
    """
    Get the file of emissions which drive a gas' first latitudinal gradient PC

    Parameters
    ----------
    ghg
        Gas of interest

    Returns
    -------
    :
        File, relative to the bundle's root (named as in the original run)
    """
    if ghg in TWO_EOF_GASES:
        return SINGLE_VARIABLE_DIR / f"{ghg}_eof-one-scaling.feather"

    if ghg in SINGLE_PROJECTION_GASES:
        return (
            Path("data")
            / "interim"
            / "inverse-emissions"
            / f"single-concentration-projection_{ghg}_inverse-emissions.feather"
        )

    return SINGLE_VARIABLE_DIR / f"{ghg}_total.feather"


def get_lat_gradient_info_file(ghg: str) -> Path:
    """Get the re-run latitudinal gradient info file, relative to the bundle's root"""
    return RERUN_DIR / "latitudinal-gradient" / f"{ghg}_latitudinal-gradient-info.nc"


def get_pcs_file(ghg: str) -> Path:
    """Get the re-run annual latitudinal gradient PCs, relative to the bundle's root"""
    return RERUN_DIR / "latitudinal-gradient" / f"{ghg}_pcs-annual.csv"


def get_regression_file(ghg: str) -> Path:
    """Get the data behind a re-derived regression, relative to the bundle's root"""
    return RERUN_DIR / "regressions" / f"{ghg}_regression.csv"


def get_regression_info_file(ghg: str) -> Path:
    """Get the result of a re-derived regression, relative to the bundle's root"""
    return RERUN_DIR / "regressions" / f"{ghg}_regression.json"


CO2_SEASONALITY_FILE = (
    RERUN_DIR / "seasonality" / "modelling-based-projection_co2_seasonality-all-time.nc"
)
"""The re-run CO2 seasonality, relative to the bundle's root"""

CO2_SEASONALITY_PC_FILE = RERUN_DIR / "seasonality" / "co2_seasonality-change-pc.csv"
"""The re-run annual CO2 seasonality change PC, relative to the bundle's root"""

CO2_SEASONALITY_REGRESSION_KEY = "co2-seasonality-change"
"""Key of the CO2 seasonality change PC's regression (against NPP) in the file names"""


EXTRACT_EMISSIONS_SCRIPT = '''
"""
Extract the emissions which drive the latitudinal gradients, for all the markers

Written by the CMIP7 ScenarioMIP GHG manuscript's build.
Uses the original run's own functions, as the original run's main flow
(`src/cmip7_scenariomip_ghg_generation/main_flow.py`) does.
"""

import json
import sys
from pathlib import Path

from cmip7_scenariomip_ghg_generation.prefect_tasks import (
    extract_fossil_biosphere_timeseries,
    extract_specific_variable_from_collection,
)
from cmip7_scenariomip_ghg_generation.scenario_info import ScenarioInfo

config = json.loads(Path(sys.argv[1]).read_text())
markers = tuple(
    ScenarioInfo(cmip_scenario_name=cmip_name, model=model, scenario=scenario)
    for model, scenario, cmip_name in config["markers"]
)
complete_emissions_files = tuple(
    Path(config["complete_emissions_dir"]) / f"{si.to_file_stem()}.feather"
    for si in markers
)

# `.fn` so we call the functions, rather than going through prefect
fossil_bio_split = extract_fossil_biosphere_timeseries.fn(
    (Path(config["raw_fossil_bio_split_file"]),),
    scenario_infos=markers,
    out_file=Path(config["fossil_bio_split_file"]),
)

for ghg, (variable_lower, out_file) in config["to_extract"].items():
    if variable_lower.endswith("|fossil"):
        extract_from = (fossil_bio_split,)
    else:
        extract_from = complete_emissions_files

    extract_specific_variable_from_collection.fn(
        extract_from=extract_from,
        scenario_infos=markers,
        variable_lower=variable_lower,
        out_file=Path(out_file),
    )
    print(f"Wrote {out_file}")
'''
"""Script which extracts the emissions, run in the bundle's environment"""


def ensure_rerun_inputs(bundle_dir: Path = BUNDLE_DIR) -> None:
    """
    Ensure the re-runs' inputs which aren't in the scenario bundle are available

    These are the historical run's interim data
    and the historical concentrations (both from the historical run's bundle).

    Parameters
    ----------
    bundle_dir
        Directory which holds the scenario bundle
    """
    ensure_bundle_available(
        files_to_get=(),
        files_to_get_tarred=("data--interim.tar.gz", "data--processed.tar.gz"),
        bundle_dir=HISTORICAL_BUNDLE_DIR,
    )
    ensure_bundle_environment(bundle_dir)


def extract_emissions(bundle_dir: Path = BUNDLE_DIR) -> dict[str, Path]:
    """
    Extract the emissions which drive the latitudinal gradients, for all the markers

    Parameters
    ----------
    bundle_dir
        Directory which holds the bundle

    Returns
    -------
    :
        File of emissions for each gas, relative to the bundle's root
    """
    to_extract = {
        ghg: (
            EOF_SCALING_VARIABLES.get(ghg, f"emissions|{ghg}"),
            get_eof_scaling_emissions_file(ghg),
        )
        for ghg in (*TWO_EOF_GASES, *ONE_EOF_GASES)
    }
    config = {
        "markers": MARKERS,
        "complete_emissions_dir": COMPLETE_EMISSIONS_DIR.as_posix(),
        "raw_fossil_bio_split_file": RAW_FOSSIL_BIO_SPLIT_FILE.as_posix(),
        "fossil_bio_split_file": FOSSIL_BIO_SPLIT_FILE.as_posix(),
        "to_extract": {
            ghg: (variable, out_file.as_posix())
            for ghg, (variable, out_file) in to_extract.items()
        },
    }

    scripts_dir = bundle_dir / RERUN_DIR / "scripts"
    scripts_dir.mkdir(parents=True, exist_ok=True)
    script = scripts_dir / "extract_emissions.py"
    script.write_text(EXTRACT_EMISSIONS_SCRIPT)
    config_file = scripts_dir / "extract_emissions_config.json"
    config_file.write_text(json.dumps(config, indent=2))

    run_in_bundle(
        ["pixi", "run", "--frozen", "python", str(script), str(config_file)],
        bundle_dir=bundle_dir,
    )

    return {ghg: out_file for ghg, (_, out_file) in to_extract.items()}


def get_markers_parameter() -> str:
    """Get the markers in the form the notebooks take them"""
    return ";;".join(f"{model};{scenario};{name}" for model, scenario, name in MARKERS)


def get_save_regression_cell(ghg_key: str) -> str:
    """
    Get a cell which saves the regression a notebook derives

    The notebooks call the regression's inputs `x` and `y`,
    the years they use `regression_years` and the result `m` and `c`.
    """
    return f"""
# Added for the CMIP7 ScenarioMIP GHG manuscript.
# The original run only plots the regression, we need its data.
import json

import pandas as pd

manuscript_regression_file = Path("{get_regression_file(ghg_key).as_posix()}")
manuscript_regression_file.parent.mkdir(exist_ok=True, parents=True)
pd.DataFrame({{"year": regression_years, "x": x.m, "y": y.m}}).to_csv(
    manuscript_regression_file, index=False
)
Path("{get_regression_info_file(ghg_key).as_posix()}").write_text(
    json.dumps(
        {{
            "x_unit": str(x.units),
            "y_unit": str(y.units),
            "slope": float(m.m),
            "slope_unit": str(m.units),
            "intercept": float(c),
        }},
        indent=2,
    )
)
manuscript_regression_file
"""


def get_save_pcs_cell(out_file: Path, pcs: dict[int, str]) -> str:
    """
    Get a cell which saves annual PCs from the year the projection starts

    Parameters
    ----------
    out_file
        File to write, relative to the bundle's root

    pcs
        Variable which holds each PC, by EOF number

        Each must have `scenario` and `year` dimensions.
    """
    pieces = "\n".join(f"    ({eof}, {variable})," for eof, variable in pcs.items())
    return f"""
# Added for the CMIP7 ScenarioMIP GHG manuscript.
# The original run only saves the monthly PCs, we want the annual ones too.
import pandas as pd

manuscript_pcs_l = []
for manuscript_eof, manuscript_pc in (
{pieces}
):
    manuscript_pc = manuscript_pc.pint.to("dimensionless").pint.dequantify()
    manuscript_pc = manuscript_pc.sel(year=manuscript_pc["year"] >= last_hist_year)
    manuscript_pc_df = (
        manuscript_pc.rename("value").to_dataframe()[["value"]].reset_index()
    )
    manuscript_pc_df["eof"] = manuscript_eof
    manuscript_pcs_l.append(manuscript_pc_df[["eof", "scenario", "year", "value"]])

manuscript_pcs_file = Path("{out_file.as_posix()}")
manuscript_pcs_file.parent.mkdir(exist_ok=True, parents=True)
pd.concat(manuscript_pcs_l).to_csv(manuscript_pcs_file, index=False)
manuscript_pcs_file
"""


def path_from_bundle_root(path: Path, bundle_dir: Path) -> str:
    """
    Get a root bath from bundle

    Parameters
    ----------
    path
        Path to convert

    bundle_dir
        Directory which holds the bundle

    Returns
    -------
    :
        `path`, relative to `bundle_dir`
    """
    return Path(os.path.relpath(path.resolve(), bundle_dir.resolve())).as_posix()


def rerun(
    executed_notebook: Path,
    rerun_name: str,
    injected_parameters: list[str],
    extra_cells: list[str],
    bundle_dir: Path = BUNDLE_DIR,
) -> Path:
    """
    Re-run a notebook with new parameters and extra cells

    Parameters
    ----------
    executed_notebook
        The notebook, as executed in the original run

    rerun_name
        Name of the re-run, used for the modified and executed notebooks' names

    injected_parameters
        Parameters to run with

    extra_cells
        Cells to append

    bundle_dir
        Directory which holds the bundle

    Returns
    -------
    :
        Path to the executed notebook
    """
    ipynb_to_run = bundle_dir / RERUN_NOTEBOOKS_DIR / f"{rerun_name}.ipynb"
    to_run = write_modified_notebook(
        start_from=executed_notebook,
        out_py=MODIFIED_NOTEBOOKS_DIR / f"{rerun_name}.py",
        out_ipynb=ipynb_to_run,
        extra_cells=extra_cells,
        injected_parameters=injected_parameters,
    )

    return run_notebook_from_bundle_dir(to_run, ipynb_to_run, bundle_dir=bundle_dir)


def rerun_co2_seasonality(bundle_dir: Path = BUNDLE_DIR) -> Path:
    """
    Re-run the CO2 seasonality projection (notebook `1011`) for all the markers

    Parameters
    ----------
    bundle_dir
        Directory which holds the bundle

    Returns
    -------
    :
        Path to the executed notebook
    """
    return rerun(
        SEASONALITY_NOTEBOOK,
        rerun_name="1011_scale-seasonality-based-on-magicc-npp_co2_markers",
        injected_parameters=[
            'ghg = "co2"',
            f'scenario_info_markers = "{get_markers_parameter()}"',
            f"harmonisation_year = {HARMONISATION_YEAR}",
            f'magicc_output_db_dir = "{MAGICC_OUTPUT_DB_DIR.as_posix()}"',
            'magicc_db_backend_str = "feather"',
            "historical_data_seasonality_lat_gradient_info_root"
            f' = "{path_from_bundle_root(HISTORICAL_INTERIM_ROOT, bundle_dir)}"',
            f'out_file = "{CO2_SEASONALITY_FILE.as_posix()}"',
        ],
        extra_cells=[
            get_save_regression_cell(CO2_SEASONALITY_REGRESSION_KEY),
            get_save_pcs_cell(CO2_SEASONALITY_PC_FILE, {0: "pc0_extended"}),
        ],
        bundle_dir=bundle_dir,
    )


def rerun_two_eof_lat_gradient(ghg: str, bundle_dir: Path = BUNDLE_DIR) -> Path:
    """
    Re-run a two-EOF latitudinal gradient projection (notebook `1031`)

    For all the markers.

    Parameters
    ----------
    ghg
        Gas of interest, one of [TWO_EOF_GASES][]

    bundle_dir
        Directory which holds the bundle

    Returns
    -------
    :
        Path to the executed notebook
    """
    return rerun(
        TWO_EOF_LAT_GRADIENT_NOTEBOOK,
        rerun_name=f"1031_scale-latitudinal-gradient-eofs_{ghg}_markers",
        injected_parameters=[
            f'ghg = "{ghg}"',
            "annual_mean_emissions_file"
            f' = "{get_eof_scaling_emissions_file(ghg).as_posix()}"',
            f"harmonisation_year = {HARMONISATION_YEAR}",
            "historical_data_root_dir"
            f' = "{path_from_bundle_root(HISTORICAL_CONCENTRATIONS_ROOT, bundle_dir)}"',
            "historical_data_seasonality_lat_gradient_info_root"
            f' = "{path_from_bundle_root(HISTORICAL_INTERIM_ROOT, bundle_dir)}"',
            f'out_file = "{get_lat_gradient_info_file(ghg).as_posix()}"',
        ],
        extra_cells=[
            get_save_regression_cell(ghg),
            get_save_pcs_cell(
                get_pcs_file(ghg), {0: "pc0_extended", 1: "pc1_extended_all_scenarios"}
            ),
        ],
        bundle_dir=bundle_dir,
    )


def rerun_one_eof_lat_gradient(ghg: str, bundle_dir: Path = BUNDLE_DIR) -> Path:
    """
    Re-run a one-EOF latitudinal gradient projection (notebook `1030`)

    For all the markers.

    Parameters
    ----------
    ghg
        Gas of interest, one of [ONE_EOF_GASES][] or [SINGLE_PROJECTION_GASES][]

    bundle_dir
        Directory which holds the bundle

    Returns
    -------
    :
        Path to the executed notebook
    """
    return rerun(
        ONE_EOF_LAT_GRADIENT_NOTEBOOK,
        rerun_name=f"1030_scale-latitudinal-gradient-based-on-emissions_{ghg}_markers",
        injected_parameters=[
            f'ghg = "{ghg}"',
            "annual_mean_emissions_file"
            f' = "{get_eof_scaling_emissions_file(ghg).as_posix()}"',
            "historical_data_root_dir"
            f' = "{path_from_bundle_root(HISTORICAL_CONCENTRATIONS_ROOT, bundle_dir)}"',
            "historical_data_seasonality_lat_gradient_info_root"
            f' = "{path_from_bundle_root(HISTORICAL_INTERIM_ROOT, bundle_dir)}"',
            f'out_file = "{get_lat_gradient_info_file(ghg).as_posix()}"',
        ],
        extra_cells=[get_save_pcs_cell(get_pcs_file(ghg), {0: "pc_extended"})],
        bundle_dir=bundle_dir,
    )


def load_rerun_pcs(ghg: str, bundle_dir: Path = BUNDLE_DIR) -> pd.DataFrame:
    """
    Load the re-run annual PCs of a gas

    Parameters
    ----------
    ghg
        Gas of interest

    bundle_dir
        Directory which holds the bundle

    Returns
    -------
    :
        PCs, with columns `eof`, `scenario`, `year` and `value`
    """
    return pd.read_csv(bundle_dir / get_pcs_file(ghg))
