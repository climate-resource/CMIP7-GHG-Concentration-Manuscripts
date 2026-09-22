"""
Generation of the N2O methods figure

The pieces this figure shares with the CH4 methods figure
live in [local.historical_ghg_forcing_for_cmip7.plotting][].
What is here is the data this figure loads,
the panels it has and how they are laid out.
"""

from __future__ import annotations

from pathlib import Path

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
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
    label_panels,
    lay_out_figure,
)
from local.historical_ghg_forcing_for_cmip7.plotting import (
    BROKEN_SPLIT,
    LAT_BIN_BOUNDS,
    MAP_ASPECT,
    add_colour_bar,
    add_network_group,
    get_decimal_year,
    get_interpolated_input_coverage_info,
    get_only_data_variable,
    ghg,
    label_name,
    plot_coverage_and_interpolated,
    plot_flying_carpet,
    plot_global_mean_extension,
    plot_global_mean_from_obs_network,
    plot_lat_gradient_pieces_from_obs_network,
    plot_monthly_means,
    plot_observation_counts,
    plot_pcs_extended,
    plot_seasonality_from_obs_network,
    plot_station_locations,
    plot_station_timeseries,
    plot_variance_explained,
    plot_yearly_means,
)

ALL_DATA_WITH_BINS_FILE = Path("manuscript-outputs") / "n2o_all-data-with-bins.csv"
"""Where the re-run notebook saves the data we want

Relative to the bundle's root directory,
because that is the notebook's working directory.
"""

LAT_GRADIENT_FULL_EOFS_PCS_FILE = (
    Path("manuscript-outputs") / "n2o_lat-gradient-full-eofs-pcs.nc"
)
"""Where the re-run notebook saves every EOF and PC of the latitudinal gradient

The original run only saves the EOFs it keeps, so this is the only place
the full decomposition is written down. We keep the whole thing, not just
the variance explained we derive from it, so the decomposition can be looked
at again without paying for another re-run.

Relative to the bundle's root directory,
because that is the notebook's working directory.
"""

LAT_GRADIENT_VARIANCE_EXPLAINED_FILE = (
    Path("manuscript-outputs") / "n2o_lat-gradient-variance-explained.csv"
)
"""Where the re-run notebook saves the variance each latitudinal gradient EOF explains

Relative to the bundle's root directory,
because that is the notebook's working directory.
"""

TITLES = {
    "timeseries": "Observation network values",
    "counts": "Obs. counts",
    "locations": "Obs. locations",
    "interpolated-most": "Interpolation: most inputs",
    "interpolated-least": "Interpolation: fewest inputs",
    "gm": "Obs. global-mean",
    "seasonality": "Obs. seasonality",
    "lat-grad-variance": "Lat. gradient EOFs variance explained",
    "lat-grad-eof": "Obs. lat. gradient EOFs",
    "lat-grad-pc": "Obs. lat. gradient PCs",
    "gm-ext": "Extended global-mean",
    "lat-grad-pc-ext": "Extended lat. gradient PCs",
    "monthly": "Monthly spatial-means",
    "yearly": "Yearly spatial-means",
    "flying-carpet": "Native resolution",
}
"""Title of each panel"""

ROWS = (
    Row(
        panels=(
            Panel("timeseries", width=2.0, colour_bar=True),
            Panel("counts", colour_bar=True),
        ),
        height=3.6,
    ),
    Row(
        panels=tuple(
            Panel(name, aspect=MAP_ASPECT, projection=ccrs.PlateCarree())
            for name in ("locations", "interpolated-most", "interpolated-least")
        ),
    ),
    Row(
        panels=(
            Panel("gm"),
            Panel("seasonality"),
            Panel("gm-ext", width=2.0, broken=True, broken_split=BROKEN_SPLIT),
        ),
        height=2.3,
    ),
    Row(
        panels=(
            Panel("lat-grad-variance"),
            Panel("lat-grad-eof"),
            Panel("lat-grad-pc"),
            Panel("lat-grad-pc-ext", broken=True, broken_split=BROKEN_SPLIT),
        ),
        height=2.3,
    ),
    Row(
        panels=(
            Panel("monthly"),
            Panel("yearly", width=1.5, broken=True, broken_split=BROKEN_SPLIT),
            Panel(
                "flying-carpet",
                aspect=1.0,
                projection="3d",
                colour_bar=True,
                colour_bar_height=0.6,
            ),
        ),
        height=2.8,
    ),
)
"""Layout of the figure's panels, top to bottom

The rows follow the steps of the method,
so the panels are labelled in reading order.

- The observational network: what was measured and when.
- The maps: where it was measured and how we interpolate between measurements.
  Every map is in the same row: a map's shape is fixed,
  so its row's height follows from how many maps share the row's width,
  and with all of them together we only pay for that once.
- The global-mean and the seasonality, and the global-mean extended back in time.
- The latitudinal gradient: how much of it each EOF explains,
  the EOFs themselves, their principal components,
  and those principal components extended back in time.
- The outputs, including the flying carpet,
  which is square and so sets its row's height.

Each row is laid out independently of the others,
see [local.historical_ghg_forcing_for_cmip7.layout][],
so rows need not have the same number of panels.
"""


def get_n2o_all_data_with_bins(
    bundle_dir: Path = DEFAULT_BUNDLE_DIR,
    original_run_notebooks_dir: Path = DEFAULT_ORIGINAL_RUN_NOTEBOOKS_DIR,
    force_rerun: bool = False,
) -> pd.DataFrame:
    """
    Get the N2O observational network data, as it went into the binning

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
        The observational network data, with the latitudinal and longitudinal
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
# but we need it to show the observational network in the manuscript.
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


def get_n2o_lat_gradient_variance_explained(
    bundle_dir: Path = DEFAULT_BUNDLE_DIR,
    original_run_notebooks_dir: Path = DEFAULT_ORIGINAL_RUN_NOTEBOOKS_DIR,
    force_rerun: bool = False,
) -> pd.DataFrame:
    """
    Get the variance explained by each of the latitudinal gradient's EOFs

    The original run only saves the EOFs it keeps, so the ones it drops --
    which are exactly the ones this tells us we can afford to drop --
    only ever exist inside the notebook which calculates them.
    This re-runs that notebook to get them back.

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
        Fraction of the variance explained by each EOF
    """
    out_file = bundle_dir / LAT_GRADIENT_VARIANCE_EXPLAINED_FILE
    full_eofs_pcs_file = bundle_dir / LAT_GRADIENT_FULL_EOFS_PCS_FILE
    if out_file.exists() and full_eofs_pcs_file.exists() and not force_rerun:
        logger.info(f"Using existing {out_file}")
        return pd.read_csv(out_file)

    base_notebook = (
        Path("calculate_n2o_monthly_fifteen_degree_pieces")
        / "only"
        / "1002_n2o_global-mean-latitudinal-gradient-seasonality.ipynb"
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
# The original run only ever saved the EOFs it keeps,
# but we want to show how much of the variance every EOF explains,
# so that keeping only the first few can be justified.
from pathlib import Path

import numpy as np
import pandas as pd

# The EOFs are the right singular vectors of the residuals, so they are
# orthonormal, which makes the principal components uncorrelated:
# the principal components' cross-product is `diag(D) ** 2`,
# i.e. the singular values squared with nothing off the diagonal.
pcs = full_eofs_pcs["principal-components"].transpose("year", "eof").data.m
singular_values_squared = pcs.T @ pcs

off_diagonal = singular_values_squared - np.diag(np.diag(singular_values_squared))
if not np.allclose(off_diagonal, 0.0, atol=1e-10 * np.trace(singular_values_squared)):
    msg = "The principal components are not uncorrelated, so this is not an SVD"
    raise AssertionError(msg)

variance_explained = np.diag(singular_values_squared) / np.trace(
    singular_values_squared
)

full_eofs_pcs_out_file = Path("{LAT_GRADIENT_FULL_EOFS_PCS_FILE.as_posix()}")
full_eofs_pcs_out_file.parent.mkdir(exist_ok=True, parents=True)
full_eofs_pcs.pint.dequantify().to_netcdf(full_eofs_pcs_out_file)

manuscript_out_file = Path("{LAT_GRADIENT_VARIANCE_EXPLAINED_FILE.as_posix()}")
manuscript_out_file.parent.mkdir(exist_ok=True, parents=True)
pd.DataFrame(
    {{
        "eof": full_eofs_pcs["eof"].values,
        "variance_explained_fraction": variance_explained,
    }}
).to_csv(manuscript_out_file, index=False)
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


def generate_n2o_methods_figure(  # noqa: PLR0915
    outfile: Path,
    bundle_dir: Path,
    original_run_notebooks_dir: Path = DEFAULT_ORIGINAL_RUN_NOTEBOOKS_DIR,
    force_rerun: bool = False,
) -> Path:
    """
    Generate the N2O methods figure

    Parameters
    ----------
    outfile
        File in which to write the figure

    bundle_dir
        Directory in which to keep the original run's bundle

    original_run_notebooks_dir
        The original run's `notebooks-executed` directory

    force_rerun
        Re-generate the figure, even if the output file already exists

    Returns
    -------
        `outfile`
    """
    outfile.unlink()
    if outfile.exists() and not force_rerun:
        logger.info(f"Using existing {outfile}")
        return outfile

    all_data_with_bins = add_network_group(
        get_n2o_all_data_with_bins(
            bundle_dir=bundle_dir,
            original_run_notebooks_dir=original_run_notebooks_dir,
            force_rerun=force_rerun,
        )
    )

    fig, axes = create_figure(ROWS)

    timeseries_scatter = plot_station_timeseries(all_data_with_bins, axes["timeseries"])
    plot_station_locations(all_data_with_bins, axes["locations"])
    counts_mesh = plot_observation_counts(all_data_with_bins, axes["counts"])

    latitude_colour_bar = add_colour_bar(
        fig,
        timeseries_scatter,
        cax=axes["timeseries-colour-bar"],
        label=r"latitude [$^{\circ}$N]",
        ticks=LAT_BIN_BOUNDS[::2],
    )
    # The points are drawn see-through so they don't hide each other,
    # but the colour bar should show the colours at full strength
    latitude_colour_bar.solids.set_alpha(1.0)

    add_colour_bar(
        fig,
        counts_mesh,
        cax=axes["counts-colour-bar"],
        label="Number of input data points",
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

    plot_variance_explained(
        get_n2o_lat_gradient_variance_explained(
            bundle_dir=bundle_dir,
            original_run_notebooks_dir=original_run_notebooks_dir,
            force_rerun=force_rerun,
        ),
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
        # Both PCs are flat bands which sit near the top and the bottom
        # of this panel, so the room for the legend is in the middle.
        legend_loc="center left",
        # This panel is a quarter of the figure wide,
        # which is not enough room for matplotlib's choice of year labels.
        max_ticks_per_half=3,
    )

    native_resolution = xr.load_dataset(
        bundle_dir / "data/interim/n2o/n2o_fifteen-degree_monthly.nc"
    )
    max_year = int(native_resolution["year"].max())
    flying_carpet_mesh = plot_flying_carpet(
        native_resolution.sel(year=range(max_year - 9, max_year + 1)),
        axes["flying-carpet"],
    )

    gm_monthly = xr.load_dataset(
        bundle_dir / "data/interim/n2o/n2o_global-mean_monthly.nc"
    )
    gm_monthly = gm_monthly.assign_coords(lat=["Global"])
    hm_monthly = xr.load_dataset(
        bundle_dir / "data/interim/n2o/n2o_hemispheric-mean_monthly.nc"
    )
    sh_lat = -45.0
    hm_monthly = hm_monthly.assign_coords(
        lat=[
            "Southern hemisphere" if v == sh_lat else "Northern hemisphere"
            for v in hm_monthly["lat"]
        ]
    )
    pda = xr.concat([gm_monthly, hm_monthly], "lat")
    pda = pda.rename({"lat": "Region"})
    max_year = int(pda["year"].max())
    plot_monthly_means(
        pda.sel(year=range(max_year - 4, max_year + 1)), ax=axes["monthly"]
    )

    gm_yearly = xr.load_dataset(
        bundle_dir / "data/interim/n2o/n2o_global-mean_annual-mean.nc"
    )
    gm_yearly = gm_yearly.assign_coords(lat=["Global"])
    hm_yearly = xr.load_dataset(
        bundle_dir / "data/interim/n2o/n2o_hemispheric-mean_annual-mean.nc"
    )
    hm_yearly = hm_yearly.assign_coords(
        lat=[
            "Southern hemisphere" if v == sh_lat else "Northern hemisphere"
            for v in hm_yearly["lat"]
        ]
    )
    pda = xr.concat([gm_yearly, hm_yearly], "lat")
    pda = pda.rename({"lat": "Region"})
    plot_yearly_means(pda, ax_left=axes["yearly-l"], ax_right=axes["yearly-r"])

    add_colour_bar(
        fig,
        flying_carpet_mesh,
        cax=axes["flying-carpet-colour-bar"],
        # The panel's own vertical axis carries no label, so this says
        # both what is plotted and what its units are.
        label=label_name(
            f"{ghg(all_data_with_bins)} "
            f"[{get_only_data_variable(native_resolution).attrs['units']}]"
        ),
        label_on_top=True,
    )

    label_panels(ROWS, axes, TITLES)
    # Last, because it needs to know how much room everything takes up
    lay_out_figure(fig, axes, ROWS)

    outfile.parent.mkdir(exist_ok=True, parents=True)
    logger.info(f"Writing {outfile}")
    fig.savefig(outfile)
    plt.close(fig)

    return outfile
