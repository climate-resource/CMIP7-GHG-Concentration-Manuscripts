"""
Generation of the CO2 methods figure

The pieces this figure shares with the N2O methods figure
live in [local.historical_ghg_forcing_for_cmip7.plotting][].
What is here is the data this figure loads,
the panels it has and how they are laid out.

Each gas has two methods figures.
The main one follows the method from the observations to the extended components.
The appendix one holds the panels which show the working along the way:
the interpolation at its best and worst, how much variance each EOF explains,
the principal components the observations give
and the regressions their extension leans on.
Both are drawn by the one function, because they are drawn from the same data.
"""
# Differences from ch4
# - seasonality has PCs
#   - regression against composite back to 1850, constant before
#   - seasonality delta has to be plotted too: it is per latitude

from __future__ import annotations

import json
from pathlib import Path

import cartopy.crs as ccrs
import matplotlib.axes
import numpy as np
import pandas as pd
import seaborn as sns
import xarray as xr
import yaml
from loguru import logger

from local.cmip_ghg_generation import (
    DEFAULT_BUNDLE_DIR,
    DEFAULT_ORIGINAL_RUN_NOTEBOOKS_DIR,
    MODIFIED_NOTEBOOKS_DIR,
    ensure_bundle_available,
    ensure_bundle_environment,
    ensure_executed_notebook_available,
    run_notebook_from_bundle_dir,
    write_modified_notebook,
)
from local.historical_ghg_forcing_for_cmip7.layout import (
    Panel,
    Row,
    create_figure,
    manuscript_style,
    merge_axes,
    save_figure,
)
from local.historical_ghg_forcing_for_cmip7.plotting import (
    BROKEN_SPLIT,
    LATITUDE_COLOUR_MAP,
    LATITUDE_KEY,
    LATITUDE_NORMALISATION,
    MAP_ASPECT,
    add_colour_bar,
    add_latitude_legend,
    add_network_group,
    compact_existing_legend,
    eof_palette,
    get_beside_legend_kwargs,
    get_decimal_year,
    get_interpolated_input_coverage_info,
    ghg,
    label_name,
    plot_coverage_and_interpolated,
    plot_global_mean_extension,
    plot_global_mean_from_obs_network,
    plot_lat_gradient_pieces_from_obs_network,
    plot_observation_counts,
    plot_pc_timeseries_regression,
    plot_pcs_extended,
    plot_seasonality_from_obs_network,
    plot_station_locations,
    plot_station_timeseries,
    plot_variance_explained,
)
from local.historical_ghg_forcing_for_cmip7.variance_explained import (
    DecompositionToSave,
    get_variance_explained,
)

ALL_DATA_WITH_BINS_FILE = Path("manuscript-outputs") / "co2_all-data-with-bins.csv"
"""Where the re-run notebook saves the data we want

Relative to the bundle's root directory,
because that is the notebook's working directory.
"""

PRIMAP_REGRESSION_DATA_FILE = (
    Path("manuscript-outputs") / "co2_primap-regression-data.nc"
)
"""Where the re-run notebook saves the PRIMAP regression data we want

Relative to the bundle's root directory,
because that is the notebook's working directory.
"""

PRIMAP_REGRESSION_YEARS_FILE = (
    Path("manuscript-outputs") / "co2_primap-regression-years.json"
)
"""Where the re-run notebook saves the PRIMAP regression years information we want

Relative to the bundle's root directory,
because that is the notebook's working directory.
"""

CO2_SEASONALITY_CHANGE_COMPOSITE_FILE = (
    Path("manuscript-outputs") / "co2_seasonality-change-composite.nc"
)
"""Where the re-run notebook saves the seasonality composite regression timeseries

Relative to the bundle's root directory,
because that is the notebook's working directory.
"""

CO2_SEASONALITY_CHANGE_COMPOSITE_REGRESSION_YEARS_FILE = (
    Path("manuscript-outputs")
    / "co2_seasonality-change-composite-regression-years.json"
)
"""Where the re-run notebook saves the years the composite regression fills

Relative to the bundle's root directory,
because that is the notebook's working directory.
"""
# PC0_OPTIMISED_YEARS_FILE = Path("manuscript-outputs") / "co2_pc0-optimised-years.json"
# """Where the re-run notebook saves the PC0 optimised years information we want
#
# Relative to the bundle's root directory,
# because that is the notebook's working directory.
# """

LAT_GRADIENT_DECOMPOSITION = DecompositionToSave(
    eofs_pcs_variable="lat_gradient_full_eofs_pcs",
    variance_explained_file=(
        Path("manuscript-outputs") / "co2_lat-gradient-variance-explained.csv"
    ),
    full_eofs_pcs_file=(
        Path("manuscript-outputs") / "co2_lat-gradient-full-eofs-pcs.nc"
    ),
)
"""The latitudinal gradient decomposition, as notebook 1202 leaves it"""

SEASONALITY_CHANGE_DECOMPOSITION = DecompositionToSave(
    eofs_pcs_variable="seasonality_change_full_eofs_pcs",
    variance_explained_file=(
        Path("manuscript-outputs") / "co2_seasonality-change-variance-explained.csv"
    ),
    full_eofs_pcs_file=(
        Path("manuscript-outputs") / "co2_seasonality-change-full-eofs-pcs.nc"
    ),
)
"""The seasonality change decomposition, as notebook 1202 leaves it

CO2 is the only gas whose seasonality changes over time,
so it is the only gas with this decomposition to show.
"""

DECOMPOSITIONS_NOTEBOOK = (
    Path("calculate_co2_monthly_fifteen_degree_pieces")
    / "only"
    / (
        "1202_co2_observational-network"
        "-global-mean-latitudinal-gradient-seasonality.ipynb"
    )
)
"""Notebook which calculates both of CO2's decompositions

Both come out of the one notebook, so both are fetched from one re-run.
"""

TITLES = {
    "timeseries": "Observation network values",
    "counts": "Obs. counts",
    # The empty last line is where the map's legend goes,
    # see plot_station_locations
    "locations": "Obs. locations\n",
    "gm": "Obs. global-mean",
    "seasonality": "Obs. seasonality",
    "gm-ext": "Extended global-mean",
    "seasonality-eof": "Obs. seasonality change EOF",
    "seasonality-pc-ext": "Extended seasonality change PC",
    "lat-grad-eof": "Obs. lat. gradient EOFs",
    "lat-grad-pc-ext": "Extended lat. gradient PCs",
}
"""Title of each panel of the main figure"""

ROWS = (
    Row(
        panels=(Panel("timeseries", colour_bar=True),),
        height=0.95,
    ),
    Row(
        panels=(
            Panel("counts", colour_bar=True),
            Panel("locations", aspect=MAP_ASPECT, projection=ccrs.PlateCarree()),
        ),
        height=0.75,
    ),
    Row(
        panels=(
            Panel("gm"),
            Panel("seasonality"),
            Panel("lat-grad-eof"),
        ),
        height=0.75,
    ),
    Row(
        panels=(
            Panel("seasonality-eof"),
            Panel(
                "seasonality-pc-ext",
                width=1.6,
                broken=True,
                broken_split=BROKEN_SPLIT,
            ),
        ),
        height=0.75,
    ),
    Row(
        panels=(
            Panel("gm-ext", broken=True, broken_split=BROKEN_SPLIT),
            Panel("lat-grad-pc-ext", broken=True, broken_split=BROKEN_SPLIT),
        ),
        height=0.8,
    ),
)
"""Layout of the main figure's panels, top to bottom

The rows follow the steps of the method,
so the panels are labelled in reading order.

- The observation network: what was measured and when,
  across the whole width of the figure, because it is the figure's starting point
  and has the most in it.
- How many observations go into each bin, and where they were taken.
- What the observations give for each component:
  the global-mean, the seasonality and the latitudinal gradient's EOFs.
- The seasonality change: its EOF and its principal component,
  extended back in time.
  CO2 is the only gas whose seasonality changes over time,
  so this row is the only thing this figure has which the others do not.
- The global-mean and the latitudinal gradient's principal components,
  extended back in time, which is where the method ends up.
  The latitudinal gradient goes on the right,
  so its legend (which sits beside it) has the edge of the figure to itself.

Each row is laid out independently of the others,
see [local.historical_ghg_forcing_for_cmip7.layout][],
so rows need not have the same number of panels.
"""

APPENDIX_TITLES = {
    "interpolated-most": "Interpolation: most inputs",
    "interpolated-least": "Interpolation: fewest inputs",
    "seasonality-variance": "Seasonality change EOF\nvariance explained",
    "seasonality-pc": "Obs. seasonality change PC",
    "seasonality-pc-composite": "Seasonality change PC0\nvs. composite",
    "lat-grad-variance": "Lat. gradient EOF\nvariance explained",
    "lat-grad-pc": "Obs. lat. gradient PCs",
    # Which emissions is on the panel's horizontal axis
    "lat-grad-pc-emms": "Lat. gradient PC0\nvs. emissions",
}
"""Title of each panel of the appendix figure"""

APPENDIX_ROWS = (
    Row(
        panels=tuple(
            Panel(name, aspect=MAP_ASPECT, projection=ccrs.PlateCarree())
            for name in ("interpolated-most", "interpolated-least")
        ),
    ),
    Row(
        panels=(
            Panel("seasonality-variance"),
            Panel("seasonality-pc"),
            Panel("seasonality-pc-composite"),
        ),
        height=0.85,
    ),
    Row(
        panels=(
            Panel("lat-grad-variance"),
            Panel("lat-grad-pc"),
            Panel("lat-grad-pc-emms"),
        ),
        height=0.85,
    ),
)
"""Layout of the appendix figure's panels, top to bottom

- The interpolation, with the most and the fewest inputs.
- The seasonality change: how much of it each EOF explains,
  the principal component the observations give,
  and the regression its extension leans on.
- The latitudinal gradient, in the same order.
"""


def get_co2_all_data_with_bins(
    bundle_dir: Path = DEFAULT_BUNDLE_DIR,
    original_run_notebooks_dir: Path = DEFAULT_ORIGINAL_RUN_NOTEBOOKS_DIR,
    force_rerun: bool = False,
) -> pd.DataFrame:
    """
    Get the co2 observation network data, as it went into the binning

    This re-runs the original run's binning notebook if it needs to.
    That is slow the first time (the original run's environment has to be
    downloaded and installed), so the result is re-used if it is already there.

    Parameters
    ----------
    bundle_dir
        Directory in which to keep the original run's bundle

    original_run_notebooks_dir
        The original run's `notebooks-executed` directory

        Only used if we don't already have a copy of the notebook we need.

    force_rerun
        Re-run the notebook even if its output is already there

    Returns
    -------
        The observation network data, with the latitudinal and longitudinal
        bin of each observation added
    """
    out_file = bundle_dir / ALL_DATA_WITH_BINS_FILE
    if out_file.exists() and not force_rerun:
        logger.info(f"Using existing {out_file}")
        return pd.read_csv(out_file)

    base_notebook = (
        Path("calculate_co2_monthly_fifteen_degree_pieces")
        / "only"
        / "1200_co2_bin-observational-network.ipynb"
    )

    start_from = ensure_executed_notebook_available(
        base_notebook,
        original_run_notebooks_dir=original_run_notebooks_dir,
    )
    ensure_bundle_available(
        files_to_get=(
            "pyproject.toml",
            "pixi.lock",
            "v1.0.0-config-raw.yaml",
        ),
        files_to_get_tarred=(
            "src.tar.gz",
            "data--interim.tar.gz",
        ),
        bundle_dir=bundle_dir,
    )
    ensure_bundle_environment(bundle_dir)

    save_cell = f"""
# Added for the CMIP7 GHG manuscript.
# The original run never saved `all_data_with_bins` out,
# but we need it to show the observation network in the manuscript.
from pathlib import Path

manuscript_out_file = Path("{ALL_DATA_WITH_BINS_FILE.as_posix()}")
manuscript_out_file.parent.mkdir(exist_ok=True, parents=True)
all_data_with_bins.to_csv(manuscript_out_file, index=False)
manuscript_out_file
"""

    notebook_name = base_notebook.stem
    ipynb_to_run = bundle_dir / "notebooks-rerun" / f"{notebook_name}.ipynb"
    to_run = write_modified_notebook(
        start_from=start_from,
        out_py=MODIFIED_NOTEBOOKS_DIR / f"{notebook_name}.py",
        out_ipynb=ipynb_to_run,
        extra_cells=[save_cell],
        step_config_id="only",
    )
    run_notebook_from_bundle_dir(
        to_run,
        ipynb_to_run,
        bundle_dir=bundle_dir,
    )

    return pd.read_csv(out_file)


def get_co2_primap_regression_data(
    bundle_dir: Path = DEFAULT_BUNDLE_DIR,
    original_run_notebooks_dir: Path = DEFAULT_ORIGINAL_RUN_NOTEBOOKS_DIR,
    force_rerun: bool = False,
) -> pd.DataFrame:
    """
    Get the PRIMAP data used for the co2 PC regression

    Parameters
    ----------
    bundle_dir
        Directory in which to keep the original run's bundle

    original_run_notebooks_dir
        The original run's `notebooks-executed` directory

        Only used if we don't already have a copy of the notebook we need.

    force_rerun
        Re-run the notebook even if its output is already there

    Returns
    -------
    :
        PRIMAP data used for the regression
    """
    out_file = bundle_dir / PRIMAP_REGRESSION_DATA_FILE
    if out_file.exists() and not force_rerun:
        logger.info(f"Using existing {out_file}")
        return xr.load_dataset(out_file)

    re_run_lat_grad_pc_extension_notebook(bundle_dir, original_run_notebooks_dir)
    return xr.load_dataset(out_file)


def get_co2_primap_regression_years(
    bundle_dir: Path = DEFAULT_BUNDLE_DIR,
    original_run_notebooks_dir: Path = DEFAULT_ORIGINAL_RUN_NOTEBOOKS_DIR,
    force_rerun: bool = False,
) -> list[int]:
    """
    Get the years in which PC0 is extended using a regression against emissions

    I.e. the years PRIMAP covers,
    other than those the observational network covers.
    In the years before PRIMAP starts,
    PC0 is held constant (by holding the emissions constant).

    Parameters
    ----------
    bundle_dir
        Directory in which to keep the original run's bundle

    original_run_notebooks_dir
        The original run's `notebooks-executed` directory

        Only used if we don't already have a copy of the notebook we need.

    force_rerun
        Re-run the notebook even if its output is already there

    Returns
    -------
    :
        Years in which the PRIMAP regression was used
    """
    out_file = bundle_dir / PRIMAP_REGRESSION_YEARS_FILE
    if out_file.exists() and not force_rerun:
        logger.info(f"Using existing {out_file}")
        with open(out_file) as fh:
            res = json.load(fh)

        return res

    re_run_lat_grad_pc_extension_notebook(bundle_dir, original_run_notebooks_dir)
    with open(out_file) as fh:
        res = json.load(fh)

    return res


def re_run_lat_grad_pc_extension_notebook(
    bundle_dir: Path = DEFAULT_BUNDLE_DIR,
    original_run_notebooks_dir: Path = DEFAULT_ORIGINAL_RUN_NOTEBOOKS_DIR,
) -> None:
    """
    Re-run the pc extension notebook
    """
    base_notebook = (
        Path("calculate_co2_monthly_fifteen_degree_pieces")
        / "only"
        / "1203_co2_extend-lat-gradient-pcs.ipynb"
    )

    start_from = ensure_executed_notebook_available(
        base_notebook,
        original_run_notebooks_dir=original_run_notebooks_dir,
    )
    ensure_bundle_available(
        files_to_get=(
            "pyproject.toml",
            "pixi.lock",
            "v1.0.0-config-raw.yaml",
        ),
        files_to_get_tarred=(
            "src.tar.gz",
            "data--interim.tar.gz",
        ),
        bundle_dir=bundle_dir,
    )
    ensure_bundle_environment(bundle_dir)

    # Make sure primap is there too
    primap_notebook = Path("retrieve_misc_data") / "only" / "0002_primap.ipynb"

    start_from_primap = ensure_executed_notebook_available(
        primap_notebook,
        original_run_notebooks_dir=original_run_notebooks_dir,
    )
    ipynb_to_run_primap = (
        bundle_dir / "notebooks-rerun" / f"{primap_notebook.name}.ipynb"
    )
    to_run_primap = write_modified_notebook(
        start_from=start_from_primap,
        out_py=MODIFIED_NOTEBOOKS_DIR / f"{primap_notebook.name}.py",
        out_ipynb=ipynb_to_run_primap,
        extra_cells=[],
        step_config_id="only",
    )
    run_notebook_from_bundle_dir(
        to_run_primap,
        ipynb_to_run_primap,
        bundle_dir=bundle_dir,
    )

    save_cell_primap_data = f"""
# Added for the CMIP7 GHG manuscript.
# The original run never saved these pieces out,
# but we need them for our plotting
from pathlib import Path

primap_regression_data_file = Path("{PRIMAP_REGRESSION_DATA_FILE.as_posix()}")
primap_regression_data_file.parent.mkdir(exist_ok=True, parents=True)
primap_regression_data.pint.dequantify().to_netcdf(primap_regression_data_file)
primap_regression_data_file
"""

    # Not `years_to_fill_with_regression`:
    # that also has the years before PRIMAP starts,
    # in which the emissions (so PC0) are just held constant.
    save_cell_primap_years = f"""
import json

primap_regression_years = np.setdiff1d(
    primap_fossil_co2_emissions["year"],
    pc0_obs_network_regression["year"],
)
primap_regression_years_file = Path("{PRIMAP_REGRESSION_YEARS_FILE.as_posix()}")
primap_regression_years_file.parent.mkdir(exist_ok=True, parents=True)
with open(primap_regression_years_file, "w") as fh:
    json.dump([int(v) for v in primap_regression_years], fh)

primap_regression_years_file
"""

    notebook_name = base_notebook.stem
    ipynb_to_run = bundle_dir / "notebooks-rerun" / f"{notebook_name}.ipynb"
    to_run = write_modified_notebook(
        start_from=start_from,
        out_py=MODIFIED_NOTEBOOKS_DIR / f"{notebook_name}.py",
        out_ipynb=ipynb_to_run,
        extra_cells=[
            save_cell_primap_data,
            save_cell_primap_years,
        ],
        step_config_id="only",
    )
    run_notebook_from_bundle_dir(
        to_run,
        ipynb_to_run,
        bundle_dir=bundle_dir,
    )


def get_co2_seasonality_change_composite_timeseries(
    bundle_dir: Path = DEFAULT_BUNDLE_DIR,
    original_run_notebooks_dir: Path = DEFAULT_ORIGINAL_RUN_NOTEBOOKS_DIR,
    force_rerun: bool = False,
) -> xr.Dataset:
    """
    Get the composite the seasonality change PC is regressed against

    Only in the years the observational network covers,
    i.e. the years used for the regression.

    Parameters
    ----------
    bundle_dir
        Directory in which to keep the original run's bundle

    original_run_notebooks_dir
        The original run's `notebooks-executed` directory

        Only used if we don't already have a copy of the notebook we need.

    force_rerun
        Re-run the notebook even if its output is already there

    Returns
    -------
    :
        Composite timeseries
    """
    out_file = bundle_dir / CO2_SEASONALITY_CHANGE_COMPOSITE_FILE
    if out_file.exists() and not force_rerun:
        logger.info(f"Using existing {out_file}")
        return xr.load_dataset(out_file)

    re_run_seasonality_change_pc_extension_notebook(
        bundle_dir, original_run_notebooks_dir
    )
    return xr.load_dataset(out_file)


def get_co2_seasonality_change_composite_regression_years(
    bundle_dir: Path = DEFAULT_BUNDLE_DIR,
    original_run_notebooks_dir: Path = DEFAULT_ORIGINAL_RUN_NOTEBOOKS_DIR,
    force_rerun: bool = False,
) -> list[int]:
    """
    Get the years in which the seasonality change PC comes from the composite regression

    I.e. the years the composite covers,
    other than those the observational network covers.
    In the years before the composite starts,
    the PC is held constant (by holding the composite constant).

    Parameters
    ----------
    bundle_dir
        Directory in which to keep the original run's bundle

    original_run_notebooks_dir
        The original run's `notebooks-executed` directory

        Only used if we don't already have a copy of the notebook we need.

    force_rerun
        Re-run the notebook even if its output is already there

    Returns
    -------
    :
        Years in which the composite regression was used
    """
    out_file = bundle_dir / CO2_SEASONALITY_CHANGE_COMPOSITE_REGRESSION_YEARS_FILE
    if out_file.exists() and not force_rerun:
        logger.info(f"Using existing {out_file}")
        with open(out_file) as fh:
            res = json.load(fh)

        return res

    re_run_seasonality_change_pc_extension_notebook(
        bundle_dir, original_run_notebooks_dir
    )
    with open(out_file) as fh:
        res = json.load(fh)

    return res


def re_run_seasonality_change_pc_extension_notebook(
    bundle_dir: Path = DEFAULT_BUNDLE_DIR,
    original_run_notebooks_dir: Path = DEFAULT_ORIGINAL_RUN_NOTEBOOKS_DIR,
) -> None:
    """
    Re-run the seasonality change pc extension notebook
    """
    base_notebook = (
        Path("calculate_co2_monthly_fifteen_degree_pieces")
        / "only"
        / "1205_co2_extend-seasonality-change-pcs.ipynb"
    )

    start_from = ensure_executed_notebook_available(
        base_notebook,
        original_run_notebooks_dir=original_run_notebooks_dir,
    )
    ensure_bundle_available(
        files_to_get=(
            "pyproject.toml",
            "pixi.lock",
            "v1.0.0-config-raw.yaml",
        ),
        files_to_get_tarred=(
            "src.tar.gz",
            "data--interim.tar.gz",
            "data--raw--hadcrut5.tar.gz",
        ),
        bundle_dir=bundle_dir,
    )
    ensure_bundle_environment(bundle_dir)

    save_cell_timeseries_data = f"""
# Added for the CMIP7 GHG manuscript.
# The original run never saved these pieces out,
# but we need them for our plotting
from pathlib import Path

co2_seasonality_change_composite_file = Path("{CO2_SEASONALITY_CHANGE_COMPOSITE_FILE.as_posix()}")
co2_seasonality_change_composite_file.parent.mkdir(exist_ok=True, parents=True)
regression_timeseries_same_years.pint.dequantify().to_netcdf(co2_seasonality_change_composite_file)
co2_seasonality_change_composite_file
"""  # noqa: E501

    save_cell_composite_regression_years = f"""
import json

composite_regression_years = np.setdiff1d(
    regression_timeseries["year"],
    pc0_obs_network["year"],
)
composite_regression_years_file = Path("{CO2_SEASONALITY_CHANGE_COMPOSITE_REGRESSION_YEARS_FILE.as_posix()}")
composite_regression_years_file.parent.mkdir(exist_ok=True, parents=True)
with open(composite_regression_years_file, "w") as fh:
    json.dump([int(v) for v in composite_regression_years], fh)

composite_regression_years_file
"""  # noqa: E501

    notebook_name = base_notebook.stem
    ipynb_to_run = bundle_dir / "notebooks-rerun" / f"{notebook_name}.ipynb"
    to_run = write_modified_notebook(
        start_from=start_from,
        out_py=MODIFIED_NOTEBOOKS_DIR / f"{notebook_name}.py",
        out_ipynb=ipynb_to_run,
        extra_cells=[
            save_cell_timeseries_data,
            save_cell_composite_regression_years,
        ],
        step_config_id="only",
    )
    run_notebook_from_bundle_dir(
        to_run,
        ipynb_to_run,
        bundle_dir=bundle_dir,
    )


def plot_seasonality_change_from_obs_network(
    seasonality_change: xr.Dataset,
    axes: dict[str, matplotlib.axes.Axes],
    principal_components_key: str = "principal-components",
    eofs_key: str = "eofs",
) -> matplotlib.axes.Axes:
    """
    Plot seasonality change derived from the observation network

    There is one series here for each of the twelve latitudinal bins.
    Latitude is shown with the same colour map, over the same range,
    as the timeseries panel uses, so a colour means the same latitude
    everywhere.
    """
    da_pcs = seasonality_change[principal_components_key]
    pdf_pcs = da_pcs.to_pandas().stack().rename("value").to_frame().reset_index()
    sns.scatterplot(
        pdf_pcs,
        x="year",
        y="value",
        hue="eof",
        palette=eof_palette(pdf_pcs["eof"].unique()),
        ax=axes["pc"],
    )
    axes["pc"].set_ylabel(f"[{da_pcs.attrs['units']}]", fontsize="small")
    axes["pc"].set_xlabel("year", fontsize="small")
    axes["pc"].tick_params(labelsize="small")
    compact_existing_legend(axes["pc"], loc="best")

    sc_eof = seasonality_change[eofs_key]
    pdf_l = []
    for (eof, lat), sc_eof_lat in sc_eof.groupby(["eof", "lat"]):
        tmp = sc_eof_lat.squeeze().to_pandas().rename("value").to_frame().reset_index()
        tmp["lat"] = lat
        tmp["eof"] = eof
        pdf_l.append(tmp)

    pdf = pd.concat(pdf_l)
    axes["eof"].scatter(
        pdf["month"],
        pdf["value"],
        c=pdf["lat"],
        cmap=LATITUDE_COLOUR_MAP,
        norm=LATITUDE_NORMALISATION,
        s=5.0,
        linewidths=0.0,
    )
    axes["eof"].set_ylabel(f"[{sc_eof.attrs['units']}]", fontsize="small")
    axes["eof"].set_xlabel("month", fontsize="small")
    axes["eof"].set_xticks(np.arange(1, 12 + 1, 3))
    axes["eof"].tick_params(labelsize="small")
    # Beside the panel rather than on it: the points fill most of the panel,
    # so there is nowhere on it the legend would not cover some of them
    add_latitude_legend(axes["eof"], **get_beside_legend_kwargs(axes["eof"]))


@manuscript_style
def generate_co2_methods_figure(  # noqa: PLR0915
    outfile: Path,
    appendix_outfile: Path,
    bundle_dir: Path,
    original_run_notebooks_dir: Path = DEFAULT_ORIGINAL_RUN_NOTEBOOKS_DIR,
    force_rerun: bool = False,
) -> tuple[Path, Path]:
    """
    Generate the co2 methods figure

    Parameters
    ----------
    outfile
        File in which to write the main figure

    appendix_outfile
        File in which to write the appendix figure

    bundle_dir
        Directory in which to keep the original run's bundle

    original_run_notebooks_dir
        The original run's `notebooks-executed` directory

    force_rerun
        Re-generate the figures, even if the output files already exist

    Returns
    -------
    :
        `outfile` and `appendix_outfile`
    """
    if outfile.exists() and appendix_outfile.exists() and not force_rerun:
        logger.info(f"Using existing {outfile} and {appendix_outfile}")
        return outfile, appendix_outfile

    all_data_with_bins = add_network_group(
        get_co2_all_data_with_bins(
            bundle_dir=bundle_dir,
            original_run_notebooks_dir=original_run_notebooks_dir,
            force_rerun=force_rerun,
        )
    )

    fig, main_axes = create_figure(ROWS)
    appendix_fig, appendix_axes = create_figure(APPENDIX_ROWS)
    # Everything below draws each panel into whichever figure it is in
    axes = merge_axes(main_axes, appendix_axes)

    timeseries_scatter = plot_station_timeseries(all_data_with_bins, axes["timeseries"])
    plot_station_locations(all_data_with_bins, axes["locations"])
    counts_mesh = plot_observation_counts(all_data_with_bins, axes["counts"])

    latitude_colour_bar = add_colour_bar(
        fig,
        timeseries_scatter,
        cax=axes["timeseries-colour-bar"],
        label=r"latitude [$^{\circ}$N]",
        ticks=LATITUDE_KEY,
    )
    # The points are drawn see-through so they don't hide each other,
    # but the colour bar should show the colours at full strength
    latitude_colour_bar.solids.set_alpha(1.0)

    add_colour_bar(
        fig,
        counts_mesh,
        cax=axes["counts-colour-bar"],
        label="Number of input\ndata points",
        ticks=np.arange(1, int(counts_mesh.norm.vmax) + 1),
    )

    # Both time axes cover the same period, even though the panels differ in width
    x_limits = (
        get_decimal_year(all_data_with_bins).min() - 1.0,
        get_decimal_year(all_data_with_bins).max() + 1.0,
    )
    for panel in ("timeseries", "counts"):
        axes[panel].set_xlim(x_limits)

    axes["counts"].set_xlabel("year", fontsize="small")

    interpolated_obs_file = (
        bundle_dir / "data/interim/co2/co2_observational-network_interpolated.nc"
    )
    interpolated_obs = xr.load_dataset(interpolated_obs_file)
    most_least_coverage = get_interpolated_input_coverage_info(
        all_data_with_bins, interpolated_obs
    )
    for key in ["most", "least"]:
        # No colour bar on these two: they are here to show where the input
        # points are and how far the interpolation has to reach between them,
        # which is a question about the shape of the field rather than about
        # its values. A colour bar on a map is also the most expensive thing
        # we can add, because a map's aspect ratio is fixed, so width taken
        # from it costs it height too, and its whole row with it.
        plot_coverage_and_interpolated(
            input_data=all_data_with_bins,
            interpolated=interpolated_obs,
            year_month=most_least_coverage[key],
            ax=axes[f"interpolated-{key}"],
        )

    global_mean_from_obs_network = xr.load_dataset(
        bundle_dir / "data/interim/co2/co2_observational-network_global-annual-mean.nc"
    )
    plot_global_mean_from_obs_network(global_mean_from_obs_network, axes["gm"])

    seasonality_from_obs_network = xr.load_dataset(
        bundle_dir / "data/interim/co2/co2_observational-network_seasonality.nc",
    )
    plot_seasonality_from_obs_network(
        seasonality_from_obs_network,
        axes["seasonality"],
        assumed_units="dimensionless",
    )

    seasonality_change_from_obs_network = xr.load_dataset(
        bundle_dir
        / "data/interim/co2/co2_observational-network_seasonality-change-eofs.nc"
    )
    plot_seasonality_change_from_obs_network(
        seasonality_change_from_obs_network,
        {
            "pc": axes["seasonality-pc"],
            "eof": axes["seasonality-eof"],
        },
    )

    lat_gradient_from_obs_network = xr.load_dataset(
        bundle_dir
        / "data/interim/co2/co2_observational-network_latitudinal-gradient-eofs.nc",
    )
    pcs_extended = xr.load_dataset(
        bundle_dir / "data/interim/co2/co2_allyears-lat-gradient-eofs-pcs.nc"
    )
    primap_regression_data = get_co2_primap_regression_data(
        bundle_dir=bundle_dir,
        original_run_notebooks_dir=original_run_notebooks_dir,
        force_rerun=force_rerun,
    )
    with open(
        bundle_dir / "data/interim/co2/co2_pc0-co2-fossil-emissions-regression.yaml"
    ) as fh:
        regression_info = yaml.safe_load(fh)
    # Flip sign of PC0 so it is more intuitive
    with xr.set_options(keep_attrs=True):
        lat_gradient_from_obs_network = lat_gradient_from_obs_network * xr.where(
            lat_gradient_from_obs_network["eof"] == 0, -1, 1
        )
        pcs_extended = pcs_extended * xr.where(pcs_extended["eof"] == 0, -1, 1)

    regression_info["m"][0] *= -1
    regression_info["c"][0] *= -1

    plot_lat_gradient_pieces_from_obs_network(
        lat_gradient_from_obs_network,
        {
            "pcs": axes["lat-grad-pc"],
            "eofs": axes["lat-grad-eof"],
        },
    )

    lat_gradient_variance_explained, seasonality_change_variance_explained = (
        get_variance_explained(
            DECOMPOSITIONS_NOTEBOOK,
            (LAT_GRADIENT_DECOMPOSITION, SEASONALITY_CHANGE_DECOMPOSITION),
            bundle_dir=bundle_dir,
            original_run_notebooks_dir=original_run_notebooks_dir,
            force_rerun=force_rerun,
        )
    )
    plot_variance_explained(
        lat_gradient_variance_explained,
        axes["lat-grad-variance"],
        # Taken from the EOFs the original run kept,
        # rather than hard-coded, so the panel can't disagree with its neighbours
        n_eofs_used=lat_gradient_from_obs_network["eof"].size,
    )
    plot_variance_explained(
        seasonality_change_variance_explained,
        axes["seasonality-variance"],
        n_eofs_used=seasonality_change_from_obs_network["eof"].size,
    )

    global_mean_extended = xr.load_dataset(
        bundle_dir / "data/interim/co2/co2_global-annual-mean_allyears.nc"
    )
    menking_et_al = pd.read_csv(
        bundle_dir / "data/interim/menking-et-al-2025/menking_et_al_2025.csv"
    )
    menking_et_al = menking_et_al[menking_et_al["gas"] == ghg(all_data_with_bins)]
    menking_et_al_lat_l = menking_et_al["latitude"].unique()
    if len(menking_et_al_lat_l) != 1:
        raise AssertionError(menking_et_al_lat_l)
    menking_et_al_lat = menking_et_al_lat_l[0]

    mauna_loa_merged = pd.read_csv(
        bundle_dir / "data/interim/mauna_loa/merged_ice_core.csv"
    ).rename({"time": "year"}, axis="columns")
    mauna_loa_start = 1959
    mauna_loa_merged = mauna_loa_merged[mauna_loa_merged["year"] >= mauna_loa_start]
    plot_global_mean_extension(
        global_mean_extended,
        axes["gm-ext-l"],
        axes["gm-ext-r"],
        input_sources={
            "Mauna Loa + Law Dome (merged)": mauna_loa_merged,
            f"Menking et al. ({menking_et_al_lat:.2f}" + r"$^{\circ}$N)": menking_et_al,
        },
    )

    # Re-run notebook to save composite timeseries
    co2_seasonality_change_composite = get_co2_seasonality_change_composite_timeseries(
        bundle_dir=bundle_dir,
        original_run_notebooks_dir=original_run_notebooks_dir,
        force_rerun=force_rerun,
    )
    with open(
        bundle_dir / "data/interim/co2/co2_seasonality-change_temp-conc-regression.yaml"
    ) as fh:
        regression_info_seasonality_change_raw = yaml.safe_load(fh)

    regression_info_seasonality_change = regression_info_seasonality_change_raw[
        "regression_result"
    ]

    plot_pc_timeseries_regression(
        seasonality_change_from_obs_network,
        co2_seasonality_change_composite,
        # Over two lines, so it is no wider than the panel above it
        timeseries_name=label_name("Temperature - co2 concentration\ncomposite"),
        regression_info=regression_info_seasonality_change,
        ax=axes["seasonality-pc-composite"],
        x_unit="dimensionless",
    )

    seasonality_change_extended = xr.load_dataset(
        bundle_dir / "data/interim/co2/co2_allyears-seasonality-change-eofs-pcs.nc"
    )

    obs_based_years = seasonality_change_from_obs_network["year"].values
    composite_regression_years = np.array(
        get_co2_seasonality_change_composite_regression_years(
            bundle_dir=bundle_dir,
            original_run_notebooks_dir=original_run_notebooks_dir,
            force_rerun=force_rerun,
        )
    )
    seasonality_change_pc_constant_years = seasonality_change_extended["year"][
        ~np.isin(seasonality_change_extended["year"], obs_based_years)
        & ~np.isin(seasonality_change_extended["year"], composite_regression_years)
    ].values

    plot_pcs_extended(
        seasonality_change_extended,
        axes["seasonality-pc-ext-l"],
        axes["seasonality-pc-ext-r"],
        split_year=1800,
        pieces={
            0: {
                "Obs.": obs_based_years,
                "Composite regression": composite_regression_years,
                "Simple extrapolation": seasonality_change_pc_constant_years,
            },
        },
        # The panel's halves are narrow,
        # which is not enough room for matplotlib's choice of year labels.
        max_ticks_per_half=3,
    )

    plot_pc_timeseries_regression(
        lat_gradient_from_obs_network,
        primap_regression_data,
        # Over two lines, so it is no wider than the panel above it
        timeseries_name=label_name("co2 emissions of\ngeological origin"),
        regression_info=regression_info,
        ax=axes["lat-grad-pc-emms"],
        x_unit="GtC / yr",
    )

    obs_based_years = lat_gradient_from_obs_network["year"].values
    primap_regression_years = np.array(
        get_co2_primap_regression_years(
            bundle_dir=bundle_dir,
            original_run_notebooks_dir=original_run_notebooks_dir,
            force_rerun=force_rerun,
        )
    )

    pc0_constant_years = pcs_extended["year"].values[
        ~np.isin(pcs_extended["year"], obs_based_years)
        & ~np.isin(pcs_extended["year"], primap_regression_years)
    ]
    pc1_constant_years = pcs_extended["year"].values[
        ~np.isin(pcs_extended["year"], obs_based_years)
    ]

    plot_pcs_extended(
        pcs_extended,
        axes["lat-grad-pc-ext-l"],
        axes["lat-grad-pc-ext-r"],
        split_year=1800,
        pieces={
            0: {
                "Obs.": obs_based_years,
                "Emissions regression": primap_regression_years,
                "Simple extrapolation": pc0_constant_years,
            },
            1: {
                "Obs.": obs_based_years,
                "Simple extrapolation": pc1_constant_years,
            },
        },
        # The panel's halves are narrow,
        # which is not enough room for matplotlib's choice of year labels.
        max_ticks_per_half=3,
    )

    return (
        save_figure(fig, main_axes, ROWS, TITLES, outfile),
        save_figure(
            appendix_fig,
            appendix_axes,
            APPENDIX_ROWS,
            APPENDIX_TITLES,
            appendix_outfile,
        ),
    )
