"""
Generation of the N2O methods figure

The pieces this figure shares with the CH4 methods figure
live in [local.historical_ghg_forcing_for_cmip7.plotting][].
What is here is the data this figure loads,
the panels it has and how they are laid out.
"""

from __future__ import annotations

import string
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
from local.historical_ghg_forcing_for_cmip7.plotting import (
    FIGURE_WIDTH,
    LAT_BIN_BOUNDS,
    MAP_PANELS,
    ROW_HEIGHT,
    add_colour_bar,
    add_colour_bar_beside,
    add_network_group,
    close_broken_axis_pairs,
    fit_rows_to_fixed_aspect_panels,
    get_decimal_year,
    get_interpolated_input_coverage_info,
    get_only_data_variable,
    ghg,
    label_name,
    plot_coverage_and_interpolated,
    plot_flying_carpet,
    plot_global_mean_extension,
    plot_global_mean_from_obs_network,
    plot_lat_gradient_pcs_extended,
    plot_lat_gradient_pieces_from_obs_network,
    plot_monthly_means,
    plot_observation_counts,
    plot_seasonality_from_obs_network,
    plot_station_locations,
    plot_station_timeseries,
    plot_yearly_means,
    tuck_colour_bars_against_their_panels,
    unit,
)

ALL_DATA_WITH_BINS_FILE = Path("manuscript-outputs") / "n2o_all-data-with-bins.csv"
"""Where the re-run notebook saves the data we want

Relative to the bundle's root directory,
because that is the notebook's working directory.
"""

PANELS = (
    ("timeseries", "Observational network values"),
    ("locations", "Obs. locations"),
    ("counts", "Obs. counts"),
    ("interpolated-most", "Interpolation: most inputs"),
    ("interpolated-least", "Interpolation: fewest inputs"),
    ("gm", "Obs. global-mean"),
    ("seasonality", "Obs. seasonality"),
    ("lat-grad-eof", "Obs. lat. gradient EOFs"),
    ("lat-grad-pc", "Obs. lat. gradient PCs"),
    ("gm-ext-l", "Extended global-mean"),
    ("lat-grad-pc-ext-l", "Extended lat. gradient. PCs"),
    ("flying-carpet", "Native resolution"),
    ("monthly", "Monthly spatial-means"),
    ("yearly-l", "Yearly spatial-means"),
)
"""The figure's panels, in the order in which they are labelled

The labels follow the order of the steps in the method,
which is not the order in which the panels are laid out.
"""

MOSAIC = [
    ["timeseries", "timeseries", "interpolated-most", "gm"],
    ["timeseries", "timeseries", "interpolated-least", "seasonality"],
    ["locations", "counts", "lat-grad-eof", "lat-grad-pc"],
    ["gm-ext-l", "gm-ext-r", "lat-grad-pc-ext-l", "lat-grad-pc-ext-r"],
    ["flying-carpet", "monthly", "yearly-l", "yearly-r"],
]
"""Layout of the figure's panels

Four columns, because the columns have to be wide enough
to carry a panel's labels and its colour bar.
With more (hence narrower) columns, the layout engine runs out of room
and gives up, which is what leaves the panels sitting
in matplotlib's default positions rather than in the positions we asked for.
If a panel needs a width that four columns cannot express,
it is cheaper to move panels between rows than to split the columns further.

The two halves of each broken axis (the panels whose names end in `-l` and
`-r`) each take a column of their own here,
then have the gap between them closed up once the figure has been laid out,
see [close_broken_axis_pairs][].

There is no point setting height ratios here:
the rows which hold a map are shrunk onto the height their map needs
while the figure is being laid out,
see [fit_rows_to_fixed_aspect_panels][].
"""

FIGURE_SIZE = (FIGURE_WIDTH, len(MOSAIC) * ROW_HEIGHT)
"""Size of the figure in inches

The width is set by the page,
the height by how many rows of panels this gas needs.
"""

BROKEN_AXIS_PAIRS = (
    ("gm-ext-l", "gm-ext-r"),
    ("lat-grad-pc-ext-l", "lat-grad-pc-ext-r"),
    ("yearly-l", "yearly-r"),
)
"""The figure's broken axes, each as its left half and its right half

These are the panels which show a long record and a short one at once,
so the two halves are one panel as far as a reader is concerned.
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
        Re-run the original run's notebook even if its output is already there

    Returns
    -------
        `outfile`
    """
    all_data_with_bins = add_network_group(
        get_n2o_all_data_with_bins(
            bundle_dir=bundle_dir,
            original_run_notebooks_dir=original_run_notebooks_dir,
            force_rerun=force_rerun,
        )
    )

    fig, axes = plt.subplot_mosaic(
        MOSAIC,
        figsize=FIGURE_SIZE,
        per_subplot_kw={
            **{panel: {"projection": ccrs.PlateCarree()} for panel in MAP_PANELS},
            "flying-carpet": {"projection": "3d"},
        },
        layout="constrained",
    )

    timeseries_scatter = plot_station_timeseries(all_data_with_bins, axes["timeseries"])
    plot_station_locations(all_data_with_bins, axes["locations"])
    counts_mesh = plot_observation_counts(all_data_with_bins, axes["counts"])

    latitude_colour_bar = add_colour_bar(
        fig,
        timeseries_scatter,
        ax=axes["timeseries"],
        label=r"latitude [$^{\circ}$N]",
        ticks=LAT_BIN_BOUNDS[::2],
    )
    # The points are drawn see-through so they don't hide each other,
    # but the colour bar should show the colours at full strength
    latitude_colour_bar.solids.set_alpha(1.0)

    counts_colour_bar = add_colour_bar(
        fig,
        counts_mesh,
        ax=axes["counts"],
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
    coverage_colour_bars_axes = []
    for key in ["most", "least"]:
        ax = axes[f"interpolated-{key}"]
        coverage_mesh = plot_coverage_and_interpolated(
            input_data=all_data_with_bins,
            interpolated=interpolated_obs,
            year_month=most_least_coverage[key],
            ax=ax,
        )
        colour_bar = add_colour_bar(
            fig,
            coverage_mesh,
            ax=ax,
            label=f"[{unit(all_data_with_bins)}]",
        )

        coverage_colour_bars_axes.append([colour_bar, ax])

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
    plot_lat_gradient_pieces_from_obs_network(
        lat_gradient_from_obs_network,
        {
            "pcs": axes["lat-grad-pc"],
            "eofs": axes["lat-grad-eof"],
        },
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

    pcs_extended = xr.load_dataset(
        bundle_dir / "data/interim/n2o/n2o_allyears-lat-gradient-eofs-pcs.nc"
    )
    plot_lat_gradient_pcs_extended(
        pcs_extended,
        axes["lat-grad-pc-ext-l"],
        axes["lat-grad-pc-ext-r"],
        split_year=1950,
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

    for label, (panel, title) in zip(string.ascii_lowercase, PANELS):
        axes[panel].set_title(
            f"$\\bf{{({label})}}$ {title}",
            loc="left",
            fontsize="medium",
        )

    # The figure's colour bars, each with the panel it belongs to.
    # Any colour bar we add has to be listed here too,
    # otherwise it is left stranded next to its neighbour's panel.
    colour_bars = [
        (latitude_colour_bar, axes["timeseries"]),
        (counts_colour_bar, axes["counts"]),
        *((cb, ax) for cb, ax in coverage_colour_bars_axes),
    ]

    # These are last, and in this order,
    # because they each need to know how much space everything before them
    # has taken up, and they freeze the layout.
    fit_rows_to_fixed_aspect_panels(fig, axes)
    close_broken_axis_pairs(
        fig, [(axes[left], axes[right]) for left, right in BROKEN_AXIS_PAIRS]
    )
    add_colour_bar_beside(
        fig,
        flying_carpet_mesh,
        axes["flying-carpet"],
        # The panel's own vertical axis carries no label, so this says
        # both what is plotted and what its units are.
        label=label_name(
            f"{ghg(all_data_with_bins)} "
            f"[{get_only_data_variable(native_resolution).attrs['units']}]"
        ),
    )
    tuck_colour_bars_against_their_panels(fig, colour_bars)

    outfile.parent.mkdir(exist_ok=True, parents=True)
    logger.info(f"Writing {outfile}")
    fig.savefig(outfile)
    plt.close(fig)

    return outfile
