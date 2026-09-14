"""
Generation of the CH4 methods figure

The pieces this figure shares with the N2O methods figure
live in [local.historical_ghg_forcing_for_cmip7.plotting][].
What is here is the data this figure loads,
the panels it has and how they are laid out.
"""
# Differences from N2O
# - regression against PRIMAP to get back to ice core overlap 1948.
#   i.e. put the regression panel back in, combine the extended PCs panel ?
#   (Or keep the split, add another column)
# - first PC optimisation with ice cores
# - first PC constant before ice core overlap
# - note that Law Dome is different lat on the extended global-mean panel

from __future__ import annotations

import shutil
import string
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
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
from local.historical_ghg_forcing_for_cmip7.plotting import (
    FIGURE_WIDTH,
    LAT_BIN_BOUNDS,
    ROW_HEIGHT,
    add_colour_bar,
    add_colour_bar_beside,
    add_network_group,
    close_broken_axis_pairs,
    create_panels,
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
    plot_lat_gradient_pcs_emissions_regression,
    plot_lat_gradient_pcs_extended,
    plot_lat_gradient_pieces_from_obs_network,
    plot_monthly_means,
    plot_observation_counts,
    plot_seasonality_from_obs_network,
    plot_station_locations,
    plot_station_timeseries,
    plot_yearly_means,
    tuck_colour_bars_against_their_panels,
)
from local.paths import DATA_RAW_DIR

ALL_DATA_WITH_BINS_FILE = Path("manuscript-outputs") / "ch4_all-data-with-bins.csv"
"""Where the re-run notebook saves the data we want

Relative to the bundle's root directory,
because that is the notebook's working directory.
"""

PRIMAP_REGRESSION_DATA_FILE = (
    Path("manuscript-outputs") / "ch4_primap-regression-data.nc"
)
"""Where the re-run notebook saves the PRIMAP regression data we want

Relative to the bundle's root directory,
because that is the notebook's working directory.
"""

PRIMAP_REGRESSION_YEARS_FILE = (
    Path("manuscript-outputs") / "ch4_primap-regression-years.json"
)
"""Where the re-run notebook saves the PRIMAP regression years information we want

Relative to the bundle's root directory,
because that is the notebook's working directory.
"""

LAW_DOME_SMOOTHED_DATA_FILE = Path("manuscript-outputs") / "ch4_law-dome-smoothed.csv"
"""Where the re-run notebook saves the Law Dome data we want"""

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
    ("lat-grad-pc-emms", "Lat. gradient PCs against geological emissions"),
    ("flying-carpet", "Native resolution"),
    ("monthly", "Monthly spatial-means"),
    ("yearly-l", "Yearly spatial-means"),
)
"""The figure's panels, in the order in which they are labelled

The labels follow the order of the steps in the method,
which is not the order in which the panels are laid out.
"""

MOSAIC = [
    ["timeseries"] * 6 + ["gm"] * 3 + ["seasonality"] * 3,
    ["timeseries"] * 6 + ["lat-grad-eof"] * 3 + ["lat-grad-pc"] * 3,
    ["locations"] * 3
    + ["counts"] * 3
    + ["interpolated-most"] * 3
    + ["interpolated-least"] * 3,
    ["gm-ext-l"] * 3
    + ["gm-ext-r"] * 3
    + ["lat-grad-pc-ext-l"] * 3
    + ["lat-grad-pc-ext-r"] * 3,
    ["lat-grad-pc-emms"] * 4
    + ["flying-carpet"] * 2
    + ["monthly"] * 3
    + ["yearly-l"] * 2
    + ["yearly-r"] * 1,
]
"""Layout of the figure's panels

This figure has a panel the N2O figure does not.
It fits in the same five rows because the broken axes on the last row
are given only the width they need rather than a quarter of the figure each.

Twelve columns, which is four panels wide:
a panel takes three columns, and the ones which do not need three
(the flying carpet, whose 3D axes matplotlib always draws square,
however wide a cell we give it) take fewer,
so the panel beside them can have what is left over.

Four panels to a row is the most the page's width will take.
It is a limit on the number of panels, not on their width:
each one needs room for its labels, its ticks and its colour bar
whatever size it is drawn at, and a fifth panel in a row
leaves the layout engine no way to satisfy them all,
at which point it gives up and every panel lands
wherever matplotlib would have put it without us.
So a panel can only be made wider by taking width from its neighbours,
never by squeezing another panel into the row.

Every map is in the same row.
A map's aspect ratio is fixed, so a map is only half as tall as it is wide,
and its row is shrunk onto it (see [fit_rows_to_fixed_aspect_panels][]),
which takes every panel sharing that row down with it.
With all three maps in one row, that happens to one row rather than three,
and the row it happens to is the one where it helps:
all four panels there carry a latitude axis, and shrinking the row
is what lines those axes up with each other.

No panel may span the middle of the figure.
Column spans have to nest: a panel may span the left half or the right half,
but a panel which straddles the middle ties every column
to every other one, which the maps' fixed aspect ratio then over-constrains.

The two halves of each broken axis (the panels whose names end in `-l` and
`-r`) each take columns of their own here,
then have the gap between them closed up once the figure has been laid out,
see [close_broken_axis_pairs][].

There is no point setting height ratios here:
the row which holds the maps is shrunk onto the height they need
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


def get_ch4_all_data_with_bins(
    bundle_dir: Path = DEFAULT_BUNDLE_DIR,
    original_run_notebooks_dir: Path = DEFAULT_ORIGINAL_RUN_NOTEBOOKS_DIR,
    force_rerun: bool = False,
) -> pd.DataFrame:
    """
    Get the CH4 observational network data, as it went into the binning

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
        Path("calculate_ch4_monthly_fifteen_degree_pieces")
        / "only"
        / "1100_ch4_bin-observational-network.ipynb"
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


def get_ch4_primap_regression_data(
    bundle_dir: Path = DEFAULT_BUNDLE_DIR,
    original_run_notebooks_dir: Path = DEFAULT_ORIGINAL_RUN_NOTEBOOKS_DIR,
    force_rerun: bool = False,
) -> pd.DataFrame:
    """
    Get the PRIMAP data used for the CH4 PC regression

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
    out_file = bundle_dir / PRIMAP_REGRESSION_DATA_FILE
    if out_file.exists() and not force_rerun:
        logger.info(f"Using existing {out_file}")
        return xr.load_dataset(out_file)

    re_run_pc_extension_notebook(bundle_dir, original_run_notebooks_dir)
    return xr.load_dataset(out_file)


def re_run_pc_extension_notebook(
    bundle_dir: Path = DEFAULT_BUNDLE_DIR,
    original_run_notebooks_dir: Path = DEFAULT_ORIGINAL_RUN_NOTEBOOKS_DIR,
) -> None:
    """
    Re-run the pc extension notebook
    """
    base_notebook = (
        Path("calculate_ch4_monthly_fifteen_degree_pieces")
        / "only"
        / "1103_ch4_extend-pcs.ipynb"
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

    # Put needed file in right directory
    # TODO: try and reduce hard-coding here
    for source_file, target_file in (
        (
            (
                DATA_RAW_DIR
                / "historical-ghg-forcing-for-cmip7/zenodo-missing/law-dome_ch4_smoothed_median.csv"  # noqa: E501
            ),
            bundle_dir / "data/interim/law_dome/law-dome_ch4_smoothed_median.csv",
        ),
        (
            (
                DATA_RAW_DIR
                / "historical-ghg-forcing-for-cmip7/zenodo-missing/neem_with_location.csv"  # noqa: E501
            ),
            bundle_dir / "data/interim/neem/neem_with_location.csv",
        ),
    ):
        target_file.parent.mkdir(exist_ok=True, parents=True)
        shutil.copy2(source_file, target_file)

    save_cell_primap_data = f"""
# Added for the CMIP7 GHG manuscript.
# The original run never saved these pieces out,
# but we need them for our plotting
from pathlib import Path

years_to_fill_with_regression
primap_regression_data_file = Path("{PRIMAP_REGRESSION_DATA_FILE.as_posix()}")
primap_regression_data_file.parent.mkdir(exist_ok=True, parents=True)
primap_regression_data.pint.dequantify().to_netcdf(primap_regression_data_file)
primap_regression_data_file
"""

    save_cell_primap_years = f"""
years_to_fill_with_regression
years_to_fill_with_regression_file = Path("{PRIMAP_REGRESSION_YEARS_FILE.as_posix()}")
years_to_fill_with_regression_file.parent.mkdir(exist_ok=True, parents=True)
with open(years_to_fill_with_regression_file, "w") as fh:
    json.dump([int(v) for v in years_to_fill_with_regression], fh)

years_to_fill_with_regression_file
"""

    notebook_name = base_notebook.stem
    ipynb_to_run = bundle_dir / "notebooks-rerun" / f"{notebook_name}.ipynb"
    to_run = write_modified_notebook(
        start_from=start_from,
        out_py=MODIFIED_NOTEBOOKS_DIR / f"{notebook_name}.py",
        out_ipynb=ipynb_to_run,
        extra_cells=[save_cell_primap_data, save_cell_primap_years],
        step_config_id="only",
    )
    run_notebook_from_bundle_dir(
        to_run,
        ipynb_to_run,
        bundle_dir=bundle_dir,
    )


def generate_ch4_methods_figure(  # noqa: PLR0915
    outfile: Path,
    bundle_dir: Path,
    original_run_notebooks_dir: Path = DEFAULT_ORIGINAL_RUN_NOTEBOOKS_DIR,
    force_rerun: bool = False,
) -> Path:
    """
    Generate the CH4 methods figure

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
        get_ch4_all_data_with_bins(
            bundle_dir=bundle_dir,
            original_run_notebooks_dir=original_run_notebooks_dir,
            force_rerun=force_rerun,
        )
    )

    fig, axes = create_panels(MOSAIC, FIGURE_SIZE)

    axes["timeseries"].set_ylim([1400, 2200])
    timeseries_scatter = plot_station_timeseries(
        all_data_with_bins, axes["timeseries"], inset_y0=0.1
    )
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
        bundle_dir / "data/interim/ch4/ch4_observational-network_interpolated.nc"
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
        bundle_dir / "data/interim/ch4/ch4_observational-network_global-annual-mean.nc"
    )
    plot_global_mean_from_obs_network(global_mean_from_obs_network, axes["gm"])

    seasonality_from_obs_network = xr.load_dataset(
        bundle_dir / "data/interim/ch4/ch4_observational-network_seasonality.nc",
    )
    plot_seasonality_from_obs_network(
        seasonality_from_obs_network,
        axes["seasonality"],
        assumed_units="dimensionless",
    )

    lat_gradient_from_obs_network = xr.load_dataset(
        bundle_dir
        / "data/interim/ch4/ch4_observational-network_latitudinal-gradient-eofs.nc",
    )
    plot_lat_gradient_pieces_from_obs_network(
        lat_gradient_from_obs_network,
        {
            "pcs": axes["lat-grad-pc"],
            "eofs": axes["lat-grad-eof"],
        },
    )

    global_mean_extended = xr.load_dataset(
        bundle_dir / "data/interim/ch4/ch4_global-annual-mean_allyears.nc"
    )
    law_dome_smoothed = pd.read_csv(
        DATA_RAW_DIR
        / "historical-ghg-forcing-for-cmip7/zenodo-missing/law-dome_ch4_smoothed_median.csv"  # noqa: E501
    )
    law_dome_lat_l = law_dome_smoothed["latitude"].unique()
    if len(law_dome_lat_l) > 1:
        raise AssertionError
    law_dome_lat = law_dome_lat_l[0]

    epica = pd.read_csv(
        DATA_RAW_DIR
        / "historical-ghg-forcing-for-cmip7/zenodo-missing/epica_with_location.csv"
    )
    epica = epica[epica["year"] < law_dome_smoothed["year"].min()]
    epica_lat_l = epica["latitude"].unique()
    if len(epica_lat_l) > 1:
        raise AssertionError
    epica_lat = epica_lat_l[0]

    plot_global_mean_extension(
        global_mean_extended,
        axes["gm-ext-l"],
        axes["gm-ext-r"],
        input_sources={
            (
                f"Law Dome (smoothed, {law_dome_lat:.2f}" + r"$^{\circ}$N)"
            ): law_dome_smoothed,
            (f"EPICA ({epica_lat:.2f}" + r"$^{\circ}$N)"): epica,
        },
    )

    pcs_extended = xr.load_dataset(
        bundle_dir / "data/interim/ch4/ch4_allyears-lat-gradient-eofs-pcs.nc"
    )
    # TODO - Up to here:
    # - need to add plot of regression against PRIMAP
    #   (this is what the lat-grad-pc-emms panel is waiting for)
    # - need to colour the different bits of the PC extension
    #   - direct from obs. network
    #   - based on regression against PRIMAP
    #   - optimised to match Law Dome and NEEM
    primap_regression_data = get_ch4_primap_regression_data(
        bundle_dir=bundle_dir,
        original_run_notebooks_dir=original_run_notebooks_dir,
        force_rerun=force_rerun,
    )
    with open(
        bundle_dir / "data/interim/ch4/ch4_pc0-ch4-fossil-emissions-regression.yaml"
    ) as fh:
        regression_info = yaml.safe_load(fh)

    plot_lat_gradient_pcs_emissions_regression(
        lat_gradient_from_obs_network,
        primap_regression_data,
        emissions_name=label_name("ch4 emissions of geological origin"),
        regression_info=regression_info,
        ax=axes["lat-grad-pc-emms"],
    )
    plot_lat_gradient_pcs_extended(
        pcs_extended,
        axes["lat-grad-pc-ext-l"],
        axes["lat-grad-pc-ext-r"],
        split_year=1930,
    )

    native_resolution = xr.load_dataset(
        bundle_dir / "data/interim/ch4/ch4_fifteen-degree_monthly.nc"
    )
    max_year = int(native_resolution["year"].max())
    flying_carpet_mesh = plot_flying_carpet(
        native_resolution.sel(year=range(max_year - 9, max_year + 1)),
        axes["flying-carpet"],
    )

    gm_monthly = xr.load_dataset(
        bundle_dir / "data/interim/ch4/ch4_global-mean_monthly.nc"
    )
    gm_monthly = gm_monthly.assign_coords(lat=["Global"])
    hm_monthly = xr.load_dataset(
        bundle_dir / "data/interim/ch4/ch4_hemispheric-mean_monthly.nc"
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
        bundle_dir / "data/interim/ch4/ch4_global-mean_annual-mean.nc"
    )
    gm_yearly = gm_yearly.assign_coords(lat=["Global"])
    hm_yearly = xr.load_dataset(
        bundle_dir / "data/interim/ch4/ch4_hemispheric-mean_annual-mean.nc"
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
