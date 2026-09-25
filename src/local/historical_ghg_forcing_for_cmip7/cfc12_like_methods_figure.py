"""
Generation of the methods figure for the gases processed like CFC-12

The pieces this figure shares with the CH4 methods figure
live in [local.historical_ghg_forcing_for_cmip7.plotting][].
What is here is the data this figure loads,
the panels it has and how they are laid out.

One function serves every gas in this group,
because the method is the same for all of them
(see the CFC-12-like section of the manuscript's methods).
"""
# Differences from ch4
# - the global-mean is, for many gases, overridden by a reference source,
#   so that source has to appear on the observational network panel
#   and the network's own global-mean becomes an input to the extension
# - the lat. gradient PC is regressed against total emissions,
#   not emissions of geological origin
# - only ever one lat. gradient EOF, hence one PC

from __future__ import annotations

from pathlib import Path
from typing import Any

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
import yaml
from loguru import logger

from local.cmip_ghg_generation import (
    BUNDLE_CONFIG_FILE,
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
    plot_pc_timeseries_regression,
    plot_pcs_extended,
    plot_seasonality_from_obs_network,
    plot_station_locations,
    plot_station_timeseries,
    plot_variance_explained,
    plot_yearly_means,
)
from local.historical_ghg_forcing_for_cmip7.variance_explained import (
    DecompositionToSave,
    get_variance_explained,
)

CFC12_LIKE_GASES = (
    "c2f6",
    "c3f8",
    "ccl4",
    "cf4",
    "cfc11",
    "cfc113",
    "cfc114",
    "cfc115",
    "cfc12",
    "ch2cl2",
    "ch3br",
    "ch3ccl3",
    "ch3cl",
    "chcl3",
    "halon1211",
    "halon1301",
    "halon2402",
    "hcfc141b",
    "hcfc142b",
    "hcfc22",
    "hfc125",
    "hfc134a",
    "hfc143a",
    "hfc152a",
    "hfc227ea",
    "hfc23",
    "hfc236fa",
    "hfc245fa",
    "hfc32",
    "hfc365mfc",
    "hfc4310mee",
    "nf3",
    "sf6",
    "so2f2",
)
"""Gases which are processed the way CFC-12 is

These are the step config IDs of the original run's
`calculate_sf6_like_monthly_fifteen_degree_pieces` step,
i.e. exactly the gases this figure can be drawn for.
"""

GLOBAL_MEAN_SPLIT_MARGIN = 10
"""Years to leave before the pre-industrial year when breaking the global-mean

The extended global-mean is flat at the pre-industrial value
for every year before the pre-industrial year,
and everything the panel is about
(the point it is anchored to, the years filled by the fit, the sources)
happens after it.
So the axis is broken just before that year
rather than at some fixed year which suits no gas in particular:
the flat run-up goes in the left half, and the rest gets the right half.
This is how much of the flat run-up is kept on the right,
so the pre-industrial marker is not drawn on top of the break.
"""

BUNDLE_CONFIG_STEP = "calculate_sf6_like_monthly_fifteen_degree_pieces"
"""Step of the original run which produced these gases' pieces

Its config is where the pre-industrial value and year of each gas is recorded,
which is the one thing the figure needs
that is not written into a data file somewhere.
"""

HISTORICAL_EMISSIONS_FILE = Path("manuscript-outputs") / "historical-emissions.csv"
"""Where the re-run notebook saves the historical emissions we want

Relative to the bundle's root directory,
because that is the notebook's working directory.

This holds every gas, not one, so it is re-run once for the whole group.
"""

GLOBAL_MEAN_SUPPLEMENT_SOURCES = {
    "wmo-2022-ozone-assessment-ch7/wmo_2022_ozone_assessment_ch7.csv": (
        "WMO (2022)",
        (
            "cfc11",
            "cfc12",
            "cfc113",
            "cfc114",
            "cfc115",
            "ccl4",
            "ch3ccl3",
            "halon1211",
            "halon1301",
            "halon2402",
            "ch3br",
            "ch3cl",
        ),
    ),
    "western-et-al-2024/western_et_al_2024.csv": (
        "Western et al. (2024)",
        ("hcfc141b", "hcfc142b", "hcfc22"),
    ),
    "velders-et-al-2022/velders_et_al_2022.csv": (
        "Velders et al. (2022)",
        (
            "hfc32",
            "hfc125",
            "hfc134a",
            "hfc143a",
            "hfc152a",
            "hfc227ea",
            "hfc236fa",
            "hfc245fa",
            "hfc365mfc",
            "hfc4310mee",
        ),
    ),
    "adam-et-al-2024/adam_et_al_2024.csv": ("Adam et al. (2024)", ("hfc23",)),
    "trudinger-et-al-2016/trudinger_et_al_2016.csv": (
        "Trudinger et al. (2016)",
        ("cf4", "c2f6", "c3f8"),
    ),
}
"""Global-mean source used in place of the observational network, by file

Each entry is the source's label and the gases it is used for,
keyed by the source's file relative to the bundle's `data/interim` directory.

This mirrors `get_global_mean_supplement_config`
in the original run's `local/global_mean_extension.py`.
The gases which appear in none of these entries
take their global-mean from the observational network alone.
"""

SOURCE_BIBKEYS = {
    "WMO (2022)": "wmo_2022_ozone_ch7",
    "Western et al. (2024)": "western_2024",
    "Velders et al. (2022)": "velders_2022",
    "Adam et al. (2024)": "adam_2024",
    "Trudinger et al. (2016)": "trudinger_2016",
}
"""Bibtex key of each global-mean source, by its label

Kept next to [`GLOBAL_MEAN_SUPPLEMENT_SOURCES`][]
so a source can't be added there without being given a key here.
"""

LAT_GRADIENT_NOTEBOOK = (
    Path("calculate_sf6_like_monthly_fifteen_degree_pieces")
    / "{gas}"
    / (
        "1302_sf6-like_observational-network"
        "-global-mean-latitudinal-gradient-seasonality.ipynb"
    )
)
"""Notebook which calculates the latitudinal gradient decomposition

This step runs once per gas, so the gas goes in the path.
"""


def get_lat_gradient_decomposition(gas: str) -> DecompositionToSave:
    """
    Get the latitudinal gradient decomposition to pull out for a gas

    Parameters
    ----------
    gas
        Gas of interest

    Returns
    -------
    :
        The decomposition, as notebook 1302 leaves it
    """
    return DecompositionToSave(
        eofs_pcs_variable="lat_gradient_full_eofs_pcs",
        variance_explained_file=(
            Path("manuscript-outputs") / f"{gas}_lat-gradient-variance-explained.csv"
        ),
        full_eofs_pcs_file=(
            Path("manuscript-outputs") / f"{gas}_lat-gradient-full-eofs-pcs.nc"
        ),
    )


TITLES = {
    "timeseries": "Observation network values",
    "counts": "Obs. counts",
    "locations": "Obs. locations",
    "interpolated-most": "Interpolation: most inputs",
    "interpolated-least": "Interpolation: fewest inputs",
    "gm": "Obs. global-mean",
    "seasonality": "Obs. seasonality",
    "lat-grad-variance": "Lat. gradient EOFs variance explained",
    "lat-grad-eof": "Obs. lat. gradient EOF",
    "lat-grad-pc": "Obs. lat. gradient PC",
    "gm-ext": "Extended global-mean",
    "lat-grad-pc-emms": "Lat. gradient PC0 against total emissions",
    "lat-grad-pc-ext": "Extended lat. gradient PC",
    "monthly": "Monthly spatial-means",
    "yearly": "Yearly spatial-means",
    "flying-carpet": "Native resolution",
}
"""Title of each panel

Singular where CH4 is plural: these gases only ever use one
latitudinal gradient EOF, so there is only ever one PC to go with it.
"""

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
            Panel("lat-grad-variance", width=0.8),
            Panel("lat-grad-eof"),
            Panel("lat-grad-pc"),
            Panel("lat-grad-pc-emms"),
            Panel(
                "lat-grad-pc-ext",
                width=1.2,
                broken=True,
                broken_split=BROKEN_SPLIT,
            ),
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

The same layout as CH4 uses, see
[local.historical_ghg_forcing_for_cmip7.ch4_methods_figure][],
because these gases follow the same steps.
"""


def interim_dir(gas: str, bundle_dir: Path) -> Path:
    """
    Get the directory which holds a gas' interim data

    Parameters
    ----------
    gas
        Gas of interest

    bundle_dir
        Directory which holds the original run's bundle

    Returns
    -------
        Directory which holds `gas`' interim data
    """
    return bundle_dir / "data" / "interim" / gas


def get_cfc12_like_all_data_with_bins(gas: str, bundle_dir: Path) -> pd.DataFrame:
    """
    Get a gas' observational network data, as it went into the binning

    Unlike CO2, CH4 and N2O, the original run saved this out for these gases,
    so there is no notebook to re-run here.

    Parameters
    ----------
    gas
        Gas of interest

    bundle_dir
        Directory which holds the original run's bundle

    Returns
    -------
        The observational network data, with the latitudinal and longitudinal
        bin of each observation added
    """
    return pd.read_csv(
        interim_dir(gas, bundle_dir)
        / f"{gas}_observational-network_all-data-with-bin-information.csv"
    )


def get_historical_emissions(
    gas: str,
    bundle_dir: Path = DEFAULT_BUNDLE_DIR,
    original_run_notebooks_dir: Path = DEFAULT_ORIGINAL_RUN_NOTEBOOKS_DIR,
    force_rerun: bool = False,
) -> pd.DataFrame:
    """
    Get the historical emissions a gas' latitudinal gradient PC is regressed against

    These come from the original run's historical emissions compilation,
    which the original run kept in `data/processed` rather than `data/interim`,
    so it is not in the bundle and the notebook has to be re-run to get it.
    The notebook compiles every gas at once, so the result is shared
    by every gas in this group.

    Parameters
    ----------
    gas
        Gas of interest

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
        Historical emissions of `gas`, with a `year` column

    Raises
    ------
    AssertionError
        The compilation has no emissions for `gas`
    """
    out_file = bundle_dir / HISTORICAL_EMISSIONS_FILE
    if not out_file.exists() or force_rerun:
        re_run_historical_emissions_notebook(bundle_dir, original_run_notebooks_dir)
    else:
        logger.info(f"Using existing {out_file}")

    all_emissions = pd.read_csv(out_file)
    res = all_emissions[all_emissions["variable"] == f"Emissions|{gas}"]
    if res.empty:
        msg = f"No historical emissions for {gas=}, check {out_file}"
        raise AssertionError(msg)

    return res.rename({"time": "year"}, axis="columns")


def re_run_historical_emissions_notebook(
    bundle_dir: Path = DEFAULT_BUNDLE_DIR,
    original_run_notebooks_dir: Path = DEFAULT_ORIGINAL_RUN_NOTEBOOKS_DIR,
) -> None:
    """
    Re-run the historical emissions compilation notebook

    Parameters
    ----------
    bundle_dir
        Directory in which to keep the original run's bundle

    original_run_notebooks_dir
        The original run's `notebooks-executed` directory
    """
    base_notebook = (
        Path("compile_historical_emissions")
        / "only"
        / "0109_compile-complete-dataset.ipynb"
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
# The original run wrote this into `data/processed`,
# which is not in the bundle, but we need it
# for the latitudinal gradient PC regression panel.
from pathlib import Path

manuscript_out_file = Path("{HISTORICAL_EMISSIONS_FILE.as_posix()}")
manuscript_out_file.parent.mkdir(exist_ok=True, parents=True)
out.to_csv(manuscript_out_file, index=False)
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


def get_step_config(gas: str, bundle_dir: Path) -> dict[str, Any]:
    """
    Get the original run's config for a gas, for the step which processed it

    Parameters
    ----------
    gas
        Gas of interest

    bundle_dir
        Directory which holds the original run's bundle

    Returns
    -------
    :
        `gas`' config for [`BUNDLE_CONFIG_STEP`][], as loaded from the YAML

    Raises
    ------
    AssertionError
        The bundle's config has no entry for `gas`
    """
    with open(bundle_dir / BUNDLE_CONFIG_FILE) as fh:
        config = yaml.safe_load(fh)

    for step_config in config[BUNDLE_CONFIG_STEP]:
        if step_config["step_config_id"] == gas:
            return step_config  # type: ignore[no-any-return]

    msg = f"No {BUNDLE_CONFIG_STEP} config for {gas=}"
    raise AssertionError(msg)


def get_pre_industrial(gas: str, bundle_dir: Path) -> tuple[int, float]:
    """
    Get the pre-industrial point a gas' global-mean is extended back to

    Parameters
    ----------
    gas
        Gas of interest

    bundle_dir
        Directory which holds the original run's bundle

    Returns
    -------
    :
        The pre-industrial year and the value reached in it
    """
    pre_industrial = get_step_config(gas, bundle_dir)["pre_industrial"]

    return int(pre_industrial["year"]), float(pre_industrial["value"][0])


def supplement_replaces_obs_network(
    supplement: pd.DataFrame, max_year_extended: int
) -> bool:
    """
    Get whether a global-mean source replaces the observational network's outright

    This is the rule the original run used
    (`1304_sf6-like_create-global-annual-mean`):
    a source which reaches the end of the dataset is used as-is,
    while one which stops short of it is harmonised
    to the observational network's global-mean
    and the network's global-mean is used from there on.

    Parameters
    ----------
    supplement
        The source's data, with a `year` column

    max_year_extended
        Last year of the extended global-mean

    Returns
    -------
    :
        `True` if the source replaces the observational network's global-mean,
        `False` if it is harmonised to it
    """
    return bool(supplement["year"].max() >= max_year_extended)


def get_global_mean_supplement(
    gas: str, bundle_dir: Path
) -> tuple[str, pd.DataFrame] | None:
    """
    Get the reference global-mean used in place of the observational network's

    Parameters
    ----------
    gas
        Gas of interest

    bundle_dir
        Directory which holds the original run's bundle

    Returns
    -------
    :
        The source's label and its data, or `None` if this gas' global-mean
        comes from the observational network alone
    """
    for source_file, (label, gases) in GLOBAL_MEAN_SUPPLEMENT_SOURCES.items():
        if gas not in gases:
            continue

        supplement = pd.read_csv(bundle_dir / "data" / "interim" / source_file)

        return label, supplement[supplement["gas"] == gas]

    return None


def get_flat_regression_info(
    regression_info: dict[str, list[object]],
) -> dict[str, list[object]]:
    """
    Get regression information with a scalar gradient and intercept

    The CFC-12-like part of the original run wrote the gradient out
    one list deeper than the CO2, CH4 and N2O parts did,
    so the value has to be unwrapped before it can be used as a number.

    Parameters
    ----------
    regression_info
        Regression information, as loaded from the original run's YAML

    Returns
    -------
        `regression_info`, with a scalar value for each of `m` and `c`
    """
    res = {}
    for key, (value, units) in regression_info.items():
        while isinstance(value, list):
            if len(value) != 1:
                msg = f"Expected a single value for {key}, got {value}"
                raise AssertionError(msg)

            value = value[0]  # noqa: PLW2901

        res[key] = [value, units]

    return res


def clip_to_years(
    pdf: pd.DataFrame, max_year: float, year_column: str = "year"
) -> pd.DataFrame:
    """
    Clip a source to the years the figure covers

    Some reference sources run well past the end of the dataset
    (WMO, for one, carries projections out to 2099),
    and a panel which drew those would suggest the dataset has them too.

    Parameters
    ----------
    pdf
        Data to clip

    max_year
        Last year to keep

    year_column
        Column which holds the year

    Returns
    -------
        `pdf`, with the years after `max_year` dropped
    """
    return pdf[pdf[year_column] <= max_year]


def generate_cfc12_like_methods_figure(  # noqa: PLR0915
    gas: str,
    outfile: Path,
    bundle_dir: Path,
    original_run_notebooks_dir: Path = DEFAULT_ORIGINAL_RUN_NOTEBOOKS_DIR,
    force_rerun: bool = False,
) -> Path:
    """
    Generate the methods figure for a gas which is processed like CFC-12

    Parameters
    ----------
    gas
        Gas to draw the figure for

        Must be one of [`CFC12_LIKE_GASES`][].

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
    :
        `outfile`

    Raises
    ------
    AssertionError
        `gas` is not processed like CFC-12,
        or its latitudinal gradient has more than one EOF
    """
    if outfile.exists() and not force_rerun:
        logger.info(f"Using existing {outfile}")
        return outfile

    if gas not in CFC12_LIKE_GASES:
        msg = f"{gas=} is not processed like CFC-12, expected one of {CFC12_LIKE_GASES}"
        raise AssertionError(msg)

    gas_dir = interim_dir(gas, bundle_dir)

    all_data_with_bins = add_network_group(
        get_cfc12_like_all_data_with_bins(gas, bundle_dir)
    )

    global_mean_supplement = get_global_mean_supplement(gas, bundle_dir)

    fig, axes = create_figure(ROWS)

    # The reference global-mean, where there is one, is what the rest of the
    # figure is built on, so it belongs on the panel which says what went in.
    global_mean_sources = (
        {global_mean_supplement[0]: global_mean_supplement[1]}
        if global_mean_supplement is not None
        else None
    )
    timeseries_scatter = plot_station_timeseries(
        all_data_with_bins,
        axes["timeseries"],
        global_mean_sources=global_mean_sources,
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

    interpolated_obs = xr.load_dataset(
        gas_dir / f"{gas}_observational-network_interpolated.nc"
    )
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
        gas_dir / f"{gas}_observational-network_global-annual-mean.nc"
    )
    plot_global_mean_from_obs_network(global_mean_from_obs_network, axes["gm"])

    seasonality_from_obs_network = xr.load_dataset(
        gas_dir / f"{gas}_observational-network_seasonality.nc",
    )
    plot_seasonality_from_obs_network(
        seasonality_from_obs_network,
        axes["seasonality"],
        assumed_units="dimensionless",
    )

    lat_gradient_from_obs_network = xr.load_dataset(
        gas_dir / f"{gas}_observational-network_latitudinal-gradient-eofs.nc",
    )
    expected_n_eofs = 1
    if len(lat_gradient_from_obs_network["eof"]) != expected_n_eofs:
        msg = (
            f"{gas} has {len(lat_gradient_from_obs_network['eof'])} EOFs, "
            f"these gases should only ever have {expected_n_eofs}"
        )
        raise AssertionError(msg)

    pcs_extended = xr.load_dataset(gas_dir / f"{gas}_allyears-lat-gradient-eofs-pcs.nc")
    with open(gas_dir / f"{gas}_pc0-total-emissions-regression.yaml") as fh:
        regression_info = get_flat_regression_info(yaml.safe_load(fh))

    # Flip the PC sign so it is more intuitive:
    # the EOF then runs from negative in the south to positive in the north,
    # i.e. the way round these gases' latitudinal gradient actually points,
    # and the PC grows as emissions grow.
    with xr.set_options(keep_attrs=True):
        lat_gradient_from_obs_network = -lat_gradient_from_obs_network
        pcs_extended = -pcs_extended

    regression_info["m"][0] *= -1
    regression_info["c"][0] *= -1

    plot_lat_gradient_pieces_from_obs_network(
        lat_gradient_from_obs_network,
        {
            "pcs": axes["lat-grad-pc"],
            "eofs": axes["lat-grad-eof"],
        },
    )

    (lat_gradient_variance_explained,) = get_variance_explained(
        Path(str(LAT_GRADIENT_NOTEBOOK).format(gas=gas)),
        (get_lat_gradient_decomposition(gas),),
        step_config_id=gas,
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
        gas_dir / f"{gas}_global-annual-mean_allyears.nc"
    )
    max_year_extended = int(global_mean_extended["year"].max())

    input_sources = {}
    if global_mean_supplement is not None:
        supplement_label, supplement = global_mean_supplement
        input_sources[supplement_label] = clip_to_years(supplement, max_year_extended)

        if supplement_replaces_obs_network(supplement, max_year_extended):
            # This source replaces the observational network's global-mean
            # outright, so the network's own global-mean is an input to the
            # extension rather than the thing being extended, and the panel
            # has to show it for the reader to see what was set aside.
            input_sources["Obs. global-mean"] = (
                get_only_data_variable(global_mean_from_obs_network)
                .to_pandas()
                .rename("value")
                .to_frame()
                .reset_index()
            )

    # The extension is flat at the pre-industrial value up to the
    # pre-industrial year, then the sources take over from the first year
    # they cover. Everything in between came from a fit, and nothing in the
    # panel says so unless we say it.
    pre_industrial_year, pre_industrial_value = get_pre_industrial(gas, bundle_dir)
    if global_mean_supplement is None:
        composite_start_year = int(global_mean_from_obs_network["year"].min())
    else:
        composite_start_year = int(
            clip_to_years(global_mean_supplement[1], max_year_extended)["year"].min()
        )

    fit_period = (
        (pre_industrial_year + 1, composite_start_year - 1)
        if composite_start_year > pre_industrial_year + 1
        else None
    )

    plot_global_mean_extension(
        global_mean_extended,
        axes["gm-ext-l"],
        axes["gm-ext-r"],
        input_sources=input_sources,
        split_year=pre_industrial_year - GLOBAL_MEAN_SPLIT_MARGIN,
        pre_industrial=(pre_industrial_year, pre_industrial_value),
        fit_period=fit_period,
    )

    historical_emissions = get_historical_emissions(
        gas,
        bundle_dir=bundle_dir,
        original_run_notebooks_dir=original_run_notebooks_dir,
        force_rerun=force_rerun,
    )
    emissions_unit_l = historical_emissions["unit"].unique()
    if len(emissions_unit_l) != 1:
        raise AssertionError(emissions_unit_l)
    emissions_unit = emissions_unit_l[0]

    # The regression only ever saw the years the observational network covers
    obs_based_years = lat_gradient_from_obs_network["year"].values
    regression_years = np.intersect1d(
        obs_based_years, historical_emissions["year"].values
    )
    regression_emissions = historical_emissions[
        historical_emissions["year"].isin(regression_years)
    ]
    regression_emissions_xr = xr.Dataset(
        {
            "value": xr.DataArray(
                regression_emissions["value"].values,
                dims=("year",),
                coords={"year": regression_emissions["year"].values},
                attrs={"units": emissions_unit},
            )
        }
    )

    plot_pc_timeseries_regression(
        lat_gradient_from_obs_network,
        regression_emissions_xr,
        timeseries_name=label_name(f"{gas} total emissions"),
        regression_info=regression_info,
        ax=axes["lat-grad-pc-emms"],
        x_unit=emissions_unit,
    )

    # The PC is filled with the regression everywhere the network does not
    # reach. Before the emissions dataset starts the emissions themselves are
    # held constant, so the PC is too, which is worth saying in its own right.
    emissions_start_year = int(historical_emissions["year"].min())
    extension_years = pcs_extended["year"].values[
        ~np.isin(pcs_extended["year"], obs_based_years)
    ]
    plot_pcs_extended(
        pcs_extended,
        axes["lat-grad-pc-ext-l"],
        axes["lat-grad-pc-ext-r"],
        split_year=1930,
        pieces={
            0: {
                "Obs.": obs_based_years,
                "Emissions regression": extension_years[
                    extension_years >= emissions_start_year
                ],
                "Constant": extension_years[extension_years < emissions_start_year],
            },
        },
        # This panel now shares its row with four others,
        # which is not enough room for matplotlib's choice of year labels.
        max_ticks_per_half=3,
    )

    native_resolution = xr.load_dataset(gas_dir / f"{gas}_fifteen-degree_monthly.nc")
    max_year = int(native_resolution["year"].max())
    flying_carpet_mesh = plot_flying_carpet(
        native_resolution.sel(year=range(max_year - 9, max_year + 1)),
        axes["flying-carpet"],
    )

    gm_monthly = xr.load_dataset(gas_dir / f"{gas}_global-mean_monthly.nc")
    gm_monthly = gm_monthly.assign_coords(lat=["Global"])
    hm_monthly = xr.load_dataset(gas_dir / f"{gas}_hemispheric-mean_monthly.nc")
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

    gm_yearly = xr.load_dataset(gas_dir / f"{gas}_global-mean_annual-mean.nc")
    gm_yearly = gm_yearly.assign_coords(lat=["Global"])
    hm_yearly = xr.load_dataset(gas_dir / f"{gas}_hemispheric-mean_annual-mean.nc")
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
