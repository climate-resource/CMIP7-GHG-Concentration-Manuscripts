"""
Generation of the N2O methods figure

The pieces this figure shares with the CH4 methods figure
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

from __future__ import annotations

from pathlib import Path

import cartopy.crs as ccrs
import numpy as np
import pandas as pd
import xarray as xr
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
    LATITUDE_KEY,
    MAP_ASPECT,
    add_colour_bar,
    add_network_group,
    get_decimal_year,
    get_interpolated_input_coverage_info,
    ghg,
    plot_coverage_and_interpolated,
    plot_global_mean_extension,
    plot_global_mean_from_obs_network,
    plot_lat_gradient_pieces_from_obs_network,
    plot_observation_counts,
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

ALL_DATA_WITH_BINS_FILE = Path("manuscript-outputs") / "n2o_all-data-with-bins.csv"
"""Where the re-run notebook saves the data we want

Relative to the bundle's root directory,
because that is the notebook's working directory.
"""

LAT_GRADIENT_DECOMPOSITION = DecompositionToSave(
    eofs_pcs_variable="full_eofs_pcs",
    variance_explained_file=(
        Path("manuscript-outputs") / "n2o_lat-gradient-variance-explained.csv"
    ),
    full_eofs_pcs_file=(
        Path("manuscript-outputs") / "n2o_lat-gradient-full-eofs-pcs.nc"
    ),
)
"""The latitudinal gradient decomposition, as notebook 1002 leaves it"""

LAT_GRADIENT_NOTEBOOK = (
    Path("calculate_n2o_monthly_fifteen_degree_pieces")
    / "only"
    / "1002_n2o_global-mean-latitudinal-gradient-seasonality.ipynb"
)
"""Notebook which calculates the latitudinal gradient decomposition"""

TITLES = {
    "timeseries": "Observation network values",
    "counts": "Obs. counts",
    # The empty last line is where the map's legend goes,
    # see plot_station_locations
    "locations": "Obs. locations\n",
    "gm": "Obs. global-mean",
    "seasonality": "Obs. seasonality",
    "gm-ext": "Extended global-mean",
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
- The global-mean and the latitudinal gradient's principal components,
  extended back in time, which is where the method ends up.
  The latitudinal gradient goes on the right,
  so its legend (which sits beside it) has the edge of the figure to itself.

The same layout as CO2's main methods figure,
less CO2's row for the seasonality change,
so the panels are the same size in both and can be read against each other.

Each row is laid out independently of the others,
see [local.historical_ghg_forcing_for_cmip7.layout][],
so rows need not have the same number of panels.
"""

APPENDIX_TITLES = {
    "interpolated-most": "Interpolation: most inputs",
    "interpolated-least": "Interpolation: fewest inputs",
    "lat-grad-variance": "Lat. gradient EOF\nvariance explained",
    "lat-grad-pc": "Obs. lat. gradient PCs",
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
            Panel("lat-grad-variance"),
            Panel("lat-grad-pc"),
        ),
        height=0.85,
    ),
)
"""Layout of the appendix figure's panels, top to bottom

- The interpolation, with the most and the fewest inputs.
- The latitudinal gradient: how much of it each EOF explains
  and the principal components the observations give.
"""


def get_n2o_all_data_with_bins(
    bundle_dir: Path = DEFAULT_BUNDLE_DIR,
    original_run_notebooks_dir: Path = DEFAULT_ORIGINAL_RUN_NOTEBOOKS_DIR,
    force_rerun: bool = False,
) -> pd.DataFrame:
    """
    Get the N2O observation network data, as it went into the binning

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
        Path("calculate_n2o_monthly_fifteen_degree_pieces")
        / "only"
        / "1000_n2o_bin-observational-network.ipynb"
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


@manuscript_style
def generate_n2o_methods_figure(
    outfile: Path,
    appendix_outfile: Path,
    bundle_dir: Path,
    original_run_notebooks_dir: Path = DEFAULT_ORIGINAL_RUN_NOTEBOOKS_DIR,
    force_rerun: bool = False,
) -> tuple[Path, Path]:
    """
    Generate the N2O methods figure

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
        `outfile` and `appendix_outfile`
    """
    if outfile.exists() and appendix_outfile.exists() and not force_rerun:
        logger.info(f"Using existing {outfile} and {appendix_outfile}")
        return outfile, appendix_outfile

    all_data_with_bins = add_network_group(
        get_n2o_all_data_with_bins(
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
        bundle_dir / "data/interim/n2o/n2o_observational-network_interpolated.nc"
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
        bundle_dir / "data/interim/n2o/n2o_observational-network_global-annual-mean.nc"
    )
    plot_global_mean_from_obs_network(global_mean_from_obs_network, axes["gm"])

    seasonality_from_obs_network = xr.load_dataset(
        bundle_dir / "data/interim/n2o/n2o_observational-network_seasonality.nc",
    )
    plot_seasonality_from_obs_network(
        seasonality_from_obs_network,
        axes["seasonality"],
        assumed_units="dimensionless",
    )

    lat_gradient_from_obs_network = xr.load_dataset(
        bundle_dir
        / "data/interim/n2o/n2o_observational-network_latitudinal-gradient-eofs.nc",
    )
    pcs_extended = xr.load_dataset(
        bundle_dir / "data/interim/n2o/n2o_allyears-lat-gradient-eofs-pcs.nc"
    )
    # Flip both PC signs
    with xr.set_options(keep_attrs=True):
        lat_gradient_from_obs_network = -lat_gradient_from_obs_network
        pcs_extended = -pcs_extended

    plot_lat_gradient_pieces_from_obs_network(
        lat_gradient_from_obs_network,
        {
            "pcs": axes["lat-grad-pc"],
            "eofs": axes["lat-grad-eof"],
        },
    )

    (lat_gradient_variance_explained,) = get_variance_explained(
        LAT_GRADIENT_NOTEBOOK,
        (LAT_GRADIENT_DECOMPOSITION,),
        bundle_dir=bundle_dir,
        original_run_notebooks_dir=original_run_notebooks_dir,
        force_rerun=force_rerun,
    )
    plot_variance_explained(
        lat_gradient_variance_explained,
        axes["lat-grad-variance"],
        # Taken from the EOFs the original run kept,
        # rather than hard-coded, so the panel can't disagree with its neighbours
        n_eofs_used=lat_gradient_from_obs_network["eof"].size,
    )

    global_mean_extended = xr.load_dataset(
        bundle_dir / "data/interim/n2o/n2o_global-annual-mean_allyears.nc"
    )
    menking_et_al = pd.read_csv(
        bundle_dir / "data/interim/menking-et-al-2025/menking_et_al_2025.csv"
    )
    menking_et_al = menking_et_al[menking_et_al["gas"] == ghg(all_data_with_bins)]
    plot_global_mean_extension(
        global_mean_extended,
        axes["gm-ext-l"],
        axes["gm-ext-r"],
        input_sources={
            "Menking et al.": menking_et_al,
        },
    )

    obs_based_years = lat_gradient_from_obs_network["year"].values
    pc_constant_years = pcs_extended["year"].values[
        ~np.isin(pcs_extended["year"], obs_based_years)
    ]
    plot_pcs_extended(
        pcs_extended,
        axes["lat-grad-pc-ext-l"],
        axes["lat-grad-pc-ext-r"],
        split_year=1950,
        pieces={
            0: {
                "Obs.": obs_based_years,
                "Constant": pc_constant_years,
            },
            1: {
                "Obs.": obs_based_years,
                "Constant": pc_constant_years,
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
