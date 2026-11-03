"""
Re-running the original run's CO2 seasonality projection for VL-CF

VL-CF was run after the markers (July 2026) and its CO2 seasonality change PC
isn't in the archived bundle
(the run only wrote the seasonality itself, which wasn't archived either).
We re-run the VL-CF execution of notebook `1011` with its original parameters
(only the paths change), saving the annual PC and the regression.

This lives in its own module (rather than in `rerun_projections`)
so that adding it doesn't make pytask re-run all the other notebooks.
"""

from __future__ import annotations

from pathlib import Path

from local.scenario_ghg_forcing_for_cmip7.bundle import BUNDLE_DIR
from local.scenario_ghg_forcing_for_cmip7.rerun_projections import (
    HISTORICAL_INTERIM_ROOT,
    MAGICC_OUTPUT_DB_DIR,
    RERUN_DIR,
    SEASONALITY_NOTEBOOK,
    get_regression_file,
    get_regression_info_file,
    get_save_pcs_cell,
    get_save_regression_cell,
    path_from_bundle_root,
    rerun,
)

VL_CF_MARKER = ("REMIND-MAgPIE 3.5-4.11", "VL-cf", "vl-cf")
"""VL-CF as (model, scenario, CMIP scenario name)

From the injected parameters of the VL-CF run's notebooks.
"""

VL_CF_HARMONISATION_YEAR = 2015
"""Year in which VL-CF is harmonised to history (from the VL-CF run's notebooks)"""

CO2_SEASONALITY_VL_CF_FILE = (
    RERUN_DIR
    / "seasonality"
    / "modelling-based-projection_co2_seasonality-all-time_vl-cf.nc"
)
"""The re-run VL-CF CO2 seasonality, relative to the bundle's root"""

CO2_SEASONALITY_VL_CF_PC_FILE = (
    RERUN_DIR / "seasonality" / "co2_seasonality-change-pc_vl-cf.csv"
)
"""The re-run VL-CF annual CO2 seasonality change PC, relative to the bundle's root"""

CO2_SEASONALITY_VL_CF_REGRESSION_KEY = "co2-seasonality-change-vl-cf"
"""Key of the VL-CF CO2 seasonality change PC's regression in the file names

VL-CF's regression differs from the markers' one:
the notebook regresses against MAGICC's NPP up to the harmonisation year,
which is 2015 for VL-CF rather than 2023.
"""

RERUN_CO2_SEASONALITY_VL_CF_FILES = (
    CO2_SEASONALITY_VL_CF_PC_FILE,
    get_regression_file(CO2_SEASONALITY_VL_CF_REGRESSION_KEY),
    get_regression_info_file(CO2_SEASONALITY_VL_CF_REGRESSION_KEY),
)
"""What the VL-CF re-run writes that we use, relative to the bundle's root"""


def rerun_co2_seasonality_vl_cf(bundle_dir: Path = BUNDLE_DIR) -> Path:
    """
    Re-run the CO2 seasonality projection (notebook `1011`) for VL-CF

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
        rerun_name="1011_scale-seasonality-based-on-magicc-npp_co2_vl-cf",
        injected_parameters=[
            'ghg = "co2"',
            f'scenario_info_markers = "{";".join(VL_CF_MARKER)}"',
            f"harmonisation_year = {VL_CF_HARMONISATION_YEAR}",
            f'magicc_output_db_dir = "{MAGICC_OUTPUT_DB_DIR.as_posix()}"',
            'magicc_db_backend_str = "feather"',
            "historical_data_seasonality_lat_gradient_info_root"
            f' = "{path_from_bundle_root(HISTORICAL_INTERIM_ROOT, bundle_dir)}"',
            f'out_file = "{CO2_SEASONALITY_VL_CF_FILE.as_posix()}"',
        ],
        extra_cells=[
            get_save_regression_cell(CO2_SEASONALITY_VL_CF_REGRESSION_KEY),
            get_save_pcs_cell(CO2_SEASONALITY_VL_CF_PC_FILE, {0: "pc0_extended"}),
        ],
        bundle_dir=bundle_dir,
    )
