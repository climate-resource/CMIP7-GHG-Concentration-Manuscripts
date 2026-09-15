"""
Generation of the CO2 methods figure

The pieces this figure shares with the N2O methods figure
live in [local.historical_ghg_forcing_for_cmip7.plotting][].
What is here is the data this figure loads,
the panels it has and how they are laid out.
"""
# Differences from co2
# - global-mean extension
#   - harmonised Mauna Loa merged back to 1959
#   - harmonised Menkin et al before then (note latitude)
# - seasonality has PCs
#   - regression against composite back to 1850, constant before
#   - seasonality delta has to be plotted too: it is per latitude

from __future__ import annotations

import json
from pathlib import Path

import cartopy.crs as ccrs
import matplotlib.axes
import matplotlib.pyplot as plt
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
    label_panels,
    lay_out_figure,
)
from local.historical_ghg_forcing_for_cmip7.plotting import (
    LAT_BIN_BOUNDS,
    LATITUDE_COLOUR_MAP,
    LATITUDE_NORMALISATION,
    MAP_ASPECT,
    add_colour_bar,
    add_latitude_legend,
    add_network_group,
    compact_existing_legend,
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

# PC0_OPTIMISED_YEARS_FILE = Path("manuscript-outputs") / "co2_pc0-optimised-years.json"
# """Where the re-run notebook saves the PC0 optimised years information we want
#
# Relative to the bundle's root directory,
# because that is the notebook's working directory.
# """

TITLES = {
    "timeseries": "Observational network values",
    "counts": "Obs. counts",
    "locations": "Obs. locations",
    "interpolated-most": "Interpolation: most inputs",
    "interpolated-least": "Interpolation: fewest inputs",
    "gm": "Obs. global-mean",
    "seasonality": "Obs. seasonality",
    "seasonality-eof": "Obs. seasonality change EOF",
    "seasonality-pc": "Obs. seasonality change PC",
    "lat-grad-eof": "Obs. lat. gradient EOFs",
    "lat-grad-pc": "Obs. lat. gradient PCs",
    "gm-ext": "Extended global-mean",
    "seasonality-eof-composite": "Seasonality PC0 against composite",  # might need to put a new line in this title
    "seasonality-pc-ext": "Extended seasonality PC",
    "lat-grad-pc-emms": "Lat. gradient PC0 against geological emissions",
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
    # Might need to shuffle from here done to get things looking ok
    Row(
        panels=(
            Panel("gm"),
            Panel("seasonality"),
            Panel("seasonality-eof"),
            Panel("seasonality-pc"),
        ),
        height=2.3,
    ),
    Row(
        panels=(
            Panel("lat-grad-eof"),
            Panel("lat-grad-pc"),
            Panel("gm-ext", width=1.5, broken=True),
            Panel("seasonality-eof-composite"),
            Panel("seasonality-pc-ext", width=1.5, broken=True),
        ),
        height=2.3,
    ),
    Row(
        panels=(
            Panel("lat-grad-pc-emms"),
            Panel("lat-grad-pc-ext", width=1.5, broken=True),
            Panel("monthly"),
        ),
        height=2.3,
    ),
    Row(
        panels=(
            Panel("yearly", width=1.5, broken=True),
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
- Decomposition into a global-mean, seasonality and latitudinal gradient.
- Extending each of those back in time.
- The outputs, including the flying carpet,
  which is square and so sets its row's height.

Each row is laid out independently of the others,
see [local.historical_ghg_forcing_for_cmip7.layout][],
so rows need not have the same number of panels.
"""


def get_co2_all_data_with_bins(
    bundle_dir: Path = DEFAULT_BUNDLE_DIR,
    original_run_notebooks_dir: Path = DEFAULT_ORIGINAL_RUN_NOTEBOOKS_DIR,
    force_rerun: bool = False,
) -> pd.DataFrame:
    """
    Get the co2 observational network data, as it went into the binning

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

    re_run_pc_extension_notebook(bundle_dir, original_run_notebooks_dir)
    return xr.load_dataset(out_file)


def get_co2_primap_regression_years(
    bundle_dir: Path = DEFAULT_BUNDLE_DIR,
    original_run_notebooks_dir: Path = DEFAULT_ORIGINAL_RUN_NOTEBOOKS_DIR,
    force_rerun: bool = False,
) -> pd.DataFrame:
    """
    Get the years in which PC0 is extended using a regression against emissions

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

    re_run_pc_extension_notebook(bundle_dir, original_run_notebooks_dir)
    with open(out_file) as fh:
        res = json.load(fh)

    return res


def re_run_pc_extension_notebook(
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

    save_cell_primap_years = f"""
import json

years_to_fill_with_regression_file = Path("{PRIMAP_REGRESSION_YEARS_FILE.as_posix()}")
years_to_fill_with_regression_file.parent.mkdir(exist_ok=True, parents=True)
with open(years_to_fill_with_regression_file, "w") as fh:
    json.dump([int(v) for v in regression_years], fh)

years_to_fill_with_regression_file
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


def plot_seasonality_change_from_obs_network(
    seasonality_change: xr.Dataset,
    axes: dict[str, matplotlib.axes.Axes],
    principal_components_key: str = "principal-components",
    eofs_key: str = "eofs",
) -> matplotlib.axes.Axes:
    """
    Plot seasonality change derived from the observational network

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
        s=30.0,
        linewidths=0.0,
    )
    axes["eof"].set_ylabel(f"[{sc_eof.attrs['units']}]", fontsize="small")
    axes["eof"].set_xlabel("month", fontsize="small")
    axes["eof"].set_xticks(np.arange(1, 12 + 1, 3))
    axes["eof"].tick_params(labelsize="small")
    add_latitude_legend(axes["eof"], pdf["lat"].unique(), loc="best")


def generate_co2_methods_figure(
    outfile: Path,
    bundle_dir: Path,
    original_run_notebooks_dir: Path = DEFAULT_ORIGINAL_RUN_NOTEBOOKS_DIR,
    force_rerun: bool = False,
) -> Path:
    """
    Generate the co2 methods figure

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
        get_co2_all_data_with_bins(
            bundle_dir=bundle_dir,
            original_run_notebooks_dir=original_run_notebooks_dir,
            force_rerun=force_rerun,
        )
    )

    fig, axes = create_figure(ROWS)

    # axes["timeseries"].set_ylim([1400, 2200])
    timeseries_scatter = plot_station_timeseries(
        all_data_with_bins, axes["timeseries"], inset_y0=0.1
    )
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

    # TODO: seasonality change bits
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

    global_mean_extended = xr.load_dataset(
        bundle_dir / "data/interim/co2/co2_global-annual-mean_allyears.nc"
    )
    # law_dome_smoothed = pd.read_csv(
    #     DATA_RAW_DIR
    #     / "historical-ghg-forcing-for-cmip7/zenodo-missing/law-dome_co2_smoothed_median.csv"  # noqa: E501
    # )
    # law_dome_lat_l = law_dome_smoothed["latitude"].unique()
    # if len(law_dome_lat_l) > 1:
    #     raise AssertionError
    # law_dome_lat = law_dome_lat_l[0]
    #
    # epica = pd.read_csv(
    #     DATA_RAW_DIR
    #     / "historical-ghg-forcing-for-cmip7/zenodo-missing/epica_with_location.csv"
    # )
    # epica = epica[epica["year"] < law_dome_smoothed["year"].min()]
    # epica_lat_l = epica["latitude"].unique()
    # if len(epica_lat_l) > 1:
    #     raise AssertionError
    # epica_lat = epica_lat_l[0]
    # TODO: global-mean extension components

    plot_global_mean_extension(
        global_mean_extended,
        axes["gm-ext-l"],
        axes["gm-ext-r"],
        input_sources={
            # (
            #     f"Law Dome (smoothed, {law_dome_lat:.2f}" + r"$^{\circ}$N)"
            # ): law_dome_smoothed,
            # (f"EPICA ({epica_lat:.2f}" + r"$^{\circ}$N)"): epica,
        },
    )

    plot_lat_gradient_pcs_emissions_regression(
        lat_gradient_from_obs_network,
        primap_regression_data,
        emissions_name=label_name("co2 emissions of geological origin"),
        regression_info=regression_info,
        ax=axes["lat-grad-pc-emms"],
        x_unit="GtC / yr",
    )

    primap_years_all = np.arange(1750, 2024)

    obs_based_years = lat_gradient_from_obs_network["year"].values
    # TODO: remove hard-coding
    primap_regression_years = primap_years_all[
        ~np.isin(primap_years_all, obs_based_years)
    ]

    pc0_constant_years = pcs_extended["year"].values[
        ~np.isin(pcs_extended["year"], obs_based_years)
        & ~np.isin(pcs_extended["year"], primap_regression_years)
    ]
    pc1_constant_years = pcs_extended["year"].values[
        ~np.isin(pcs_extended["year"], obs_based_years)
    ]

    plot_lat_gradient_pcs_extended(
        pcs_extended,
        axes["lat-grad-pc-ext-l"],
        axes["lat-grad-pc-ext-r"],
        split_year=1930,
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
    )

    native_resolution = xr.load_dataset(
        bundle_dir / "data/interim/co2/co2_fifteen-degree_monthly.nc"
    )
    max_year = int(native_resolution["year"].max())
    flying_carpet_mesh = plot_flying_carpet(
        native_resolution.sel(year=range(max_year - 9, max_year + 1)),
        axes["flying-carpet"],
    )

    gm_monthly = xr.load_dataset(
        bundle_dir / "data/interim/co2/co2_global-mean_monthly.nc"
    )
    gm_monthly = gm_monthly.assign_coords(lat=["Global"])
    hm_monthly = xr.load_dataset(
        bundle_dir / "data/interim/co2/co2_hemispheric-mean_monthly.nc"
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
        bundle_dir / "data/interim/co2/co2_global-mean_annual-mean.nc"
    )
    gm_yearly = gm_yearly.assign_coords(lat=["Global"])
    hm_yearly = xr.load_dataset(
        bundle_dir / "data/interim/co2/co2_hemispheric-mean_annual-mean.nc"
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
