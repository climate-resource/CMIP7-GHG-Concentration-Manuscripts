"""
Shared pieces of the historical GHG manuscript's methods figures

The CH4 and N2O methods figures tell the same story about different gases,
so they are the same figure with different data in it.
Everything which does not depend on which gas we are plotting lives here;
what is left in each gas' module is the data it loads,
the panels it has and how they are laid out.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from typing import TYPE_CHECKING

import cartopy.crs as ccrs
import matplotlib.axes
import matplotlib.colorbar
import matplotlib.colors
import matplotlib.figure
import matplotlib.lines
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np
import openscm_units
import pandas as pd
import seaborn as sns
import xarray as xr

from local.xarray_time import convert_time_to_year_month, convert_year_month_to_time

if TYPE_CHECKING:
    import matplotlib.cm
    import matplotlib.collections

LAT_BIN_BOUNDS = np.arange(-90, 91, 15)
"""Bounds of the latitudinal bins used by the original run

This mirrors `local.binning.LAT_BIN_BOUNDS` in the original run.
"""

LAT_BIN_CENTRES = (LAT_BIN_BOUNDS[:-1] + LAT_BIN_BOUNDS[1:]) / 2.0
"""Centres of the latitudinal bins used by the original run
"""

LON_BIN_BOUNDS = np.arange(-180, 181, 60)
"""Bounds of the longitudinal bins used by the original run

This mirrors `local.binning.LON_BIN_BOUNDS` in the original run.
"""

LON_BIN_CENTRES = (LON_BIN_BOUNDS[:-1] + LON_BIN_BOUNDS[1:]) / 2.0
"""Centres of the longitudinal bins used by the original run
"""

LAT_AXIS_LIMITS = (-91, 91)
"""Limits of the latitude axis

Used by both the map and the counts panel, so the two can be read against
each other. A little room past the poles, so the polar stations have
somewhere to sit.
"""

MAP_ASPECT = (LAT_AXIS_LIMITS[1] - LAT_AXIS_LIMITS[0]) / 360.0
"""Height of a map's data box divided by its width

The maps are plate carrée, so one degree is the same length in either direction
and the map is as tall, relative to its width,
as its latitude range is relative to the 360 degrees of longitude.
"""

BROKEN_SPLIT = 0.35
"""Share of a broken panel's width to give to its left half

Every broken panel in these figures is a timeseries
which has been extended back in time,
so its left half is a long, flat run-up
and its right half is where the values actually move.
Splitting a broken panel evenly spends half of it on the run-up,
so we give the left half rather less than half.
"""

OKABE_ITO = {
    "orange": "#e69f00",
    "sky blue": "#56b4e9",
    "bluish green": "#009e73",
    "yellow": "#f0e442",
    "blue": "#0072b2",
    "vermillion": "#d55e00",
    "reddish purple": "#cc79a7",
    "black": "#000000",
}
"""The Okabe-Ito qualitative palette

Every colour these figures use to tell one series from another comes from
here, so that the figures are colour-blind safe. Naming the colours rather
than writing their hex codes at each use is what makes that checkable:
a colour which is not in here is a colour someone has to justify.

The yellow is the one entry these figures avoid: it is the palette's lightest
colour by a long way, and these panels draw small scatter markers and thin
lines on white.
"""

EOF_COLOURS = (
    OKABE_ITO["bluish green"],
    OKABE_ITO["reddish purple"],
    OKABE_ITO["orange"],
    OKABE_ITO["sky blue"],
)
"""Colour of each EOF, in order

These are from the Okabe-Ito palette, i.e. they are colour-blind safe,
as the observational network's colours are.

They are the entries that palette has left once the networks have taken
theirs: an EOF and a network can share a panel's figure, so a colour has to
mean one or the other, not both. The two we usually keep are first, and are
the two furthest apart of what remains, in lightness as well as in hue,
so they stay apart in greyscale too.

Naming them here, rather than letting each panel pick up whatever its
plotting library offers, is what keeps an EOF the same colour across the
panels which show it, and stops the figures moving if a default cycle changes.
"""

EXTENSION_SOURCE_COLOURS = {
    "Obs.": OKABE_ITO["blue"],
    "Constant": OKABE_ITO["orange"],
    "Simple extrapolation": OKABE_ITO["orange"],
    "Emissions regression": OKABE_ITO["bluish green"],
    "Composite regression": OKABE_ITO["reddish purple"],
    "Ice core optimised": OKABE_ITO["reddish purple"],
}
"""Colour to use for each source of an extended principal component

These are from the Okabe-Ito palette, i.e. they are colour-blind safe.

A source is how a year's value was arrived at, so the observations get the
one colour that is not an assumption, and the two names for holding a value
constant -- which is what "Constant" and "Simple extrapolation" both are,
under the names their gases give them -- share a colour, because they are
the same assumption.

Two sources may share a colour when no panel shows both of them, which is how
the reddish purple serves CH4's ice core optimisation and CO2's composite
regression. These colours also stand for other things in other panels of the
same figure. Both are deliberate: a panel's legend says what its colours mean,
and spending a distinct colour per meaning across a figure this size would run
the palette out long before the panels do.
"""

REGION_COLOURS = {
    "Global": OKABE_ITO["blue"],
    "Northern hemisphere": OKABE_ITO["vermillion"],
    "Southern hemisphere": OKABE_ITO["bluish green"],
}
"""Colour to use for each region of a spatial mean

The monthly and the yearly spatial-mean panels sit in the same row and show
the same three regions, one zoomed in on the end of the other. They are read
against each other, so they are given the same colours from here rather than
each picking its own up from whatever is drawing it.
"""

SOURCE_SEQUENCE_COLOURS = (
    OKABE_ITO["blue"],
    OKABE_ITO["vermillion"],
    OKABE_ITO["bluish green"],
    OKABE_ITO["orange"],
    OKABE_ITO["sky blue"],
)
"""Colours for a panel whose series are named by the caller, in order

Some panels are told what their series are at the call site: the sources
which went into a global-mean, or the stretches an extended component was
built from. Those names vary by gas, and some of them are built at run time,
so there is nothing stable to key a colour off. They take these colours in
the order the caller lists them instead.

The reddish purple is left out: it marks the fit period on the very panel
these serve, so a source drawn in it would look like one.
"""

UNUSED_EOF_COLOUR = "#b0b0b0"
"""Colour for an EOF which the workflow calculates but does not keep

The EOFs we keep are shown in their own colour, so they can be picked out
in the panels which show them; the ones we drop are all the same grey,
because the point of showing them is how little they add.
"""


NETWORK_GROUPS = {
    "NOAA": "NOAA",
    "AGAGE": "AGAGE",
    "GAGE": "AGAGE",
    "ALE": "AGAGE",
}
"""How the networks are grouped for plotting

ALE, GAGE and AGAGE are successive instruments at what are largely
the same physical sites, so we show them as one network.
The group goes by the name of the current instrument,
because the full name of all three is wider
than the panel whose legend has to carry it.
"""

MOVING_SUFFIX = " (moving)"
"""Suffix which marks a network group whose stations move

Some networks report from ships as well as from fixed sites.
A ship is a different kind of observation from a fixed site
(it has no one location to put on the map),
so it is shown as its own group rather than folded in with the fixed sites.
"""

NETWORK_GROUP_COLOURS = {
    "NOAA": OKABE_ITO["blue"],
    "NOAA (moving)": "#3a3b3a",
    "AGAGE": OKABE_ITO["vermillion"],
}
"""Colour to use for each group of observational networks

These are from the Okabe-Ito palette, i.e. they are colour-blind safe.
The exception is the moving group, which is drawn in grey:
it is a qualifier on NOAA rather than a network in its own right,
and there is no fourth Okabe-Ito colour left which reads as such.

This covers every group either gas has.
A gas whose data has no shipboard observations
simply never asks for the moving group.
"""

NETWORK_GROUP_MARKERS = {
    "NOAA": "o",
    "NOAA (moving)": ".",
    "AGAGE": "^",
}
"""Marker to use for each group of observational networks

Groups are distinguished by marker as well as colour
so the figure still works in greyscale.
"""

NETWORK_GROUP_MARKER_SIZES = {
    "NOAA": 50.0,
    "NOAA (moving)": 10.0,
    "AGAGE": 25.0,
}
"""Marker size to use for each group of observational networks on the map

Several sites host both groups, so the markers are drawn
as concentric outlines, largest first, to keep both visible.
"""

NETWORK_GROUP_COLUMN = "network_group"
"""Column we add to hold the network group"""

LATITUDE_NORMALISATION = matplotlib.colors.Normalize(vmin=-90.0, vmax=90.0)
"""Range of latitudes the latitude colour map is stretched over

Fixed to the whole globe, rather than to the latitudes a gas happens
to have been measured at, so that a colour means the same latitude
in every panel and in both gases' figures.
"""

LATITUDE_COLOUR_MAP = "coolwarm"
"""Colour map used to show latitude

Diverging about the equator, and distinguishable with the common forms
of colour vision deficiency. Deliberately different from the colour map
used for the observation counts, so the two scales don't get confused.
"""

INTERPOLATION_COLOUR_MAP = "autumn_r"

LEGEND_MARKER_COLOUR = "0.35"
"""Colour to draw legend markers in where colour means something else

In the timeseries panel colour shows latitude, so the legend
can only speak about the marker shape.
"""

GLOBAL_MEAN_SOURCE_LINESTYLES = ("-", "--", ":")
"""Line styles to tell apart global-mean sources drawn over a network

They all share a colour, so the style is the only thing left to tell them
apart. A gas rarely has more than one, so this is short on purpose.
"""

FIT_PERIOD_COLOUR = OKABE_ITO["reddish purple"]
"""Colour to draw the years of a global-mean which came from a fit

The extension is one line, but not every part of it is the same kind of
thing: some of it is a source, some of it is a fit, some of it is an
assumption. Colouring the line by which is which says so without
needing a second panel.

Not one of the first colours in matplotlib's cycle,
because those go to the extension itself and to the sources.
"""

ASSUMED_CONSTANT_COLOUR = "0.55"
"""Colour to draw the years of a global-mean which are simply assumed

Grey, because these years are not a measurement or a fit to one:
they are the pre-industrial value, held for as long as we assume it held.
"""

PRE_INDUSTRIAL_MARKER_MARGIN = 0.04
"""Room to leave under a pre-industrial marker, as a fraction of the panel"""

PRE_INDUSTRIAL_COLOUR = "k"
"""Colour to mark the pre-industrial value a global-mean is anchored to

Not on the colour cycle the sources are drawn from,
because it is a single point rather than a source.
"""

INSET_CANDIDATE_CORNERS = (
    (0.65, 0.10),
    (0.65, 0.60),
)
"""Corners the inset may be put in, as (x0, y0) in axes co-ordinates

Both on the right, because the inset magnifies the right-hand end of the
record: put it on the left and the lines joining it to the years it magnifies
are dragged across the whole panel, over the data they are meant to help read.
So the choice is between low and high on the right,
which is the choice between a gas whose values are rising
and a gas whose values are falling.

The panel's legend sits in the top left, which is why that corner is not here.

In preference order, so a gas with nothing in either corner
keeps the one the figures have always used.
"""

INSET_CLEARANCE = 0.02
"""Room to leave between the top of the inset and the data above it

As a fraction of the panel's height.
"""

ROBUST_Y_LIMIT_PERCENTILES = (0.1, 99.9)
"""Percentiles of a network's values to scale its panel by, when a few outliers
would otherwise set the scale"""

OUTLIER_HEADROOM_FRACTION = 0.15
"""How far past the robust range an outlier has to be to be left off scale

As a fraction of the robust range itself.
A network which reaches a little past its own 99.9th percentile
is a network with a long tail, and the panel shows all of it.
A network with a single reading half as far again as everything else
is a network with an outlier, and showing it costs every other point
the room it needs to be read at all.
"""

GLOBAL_MEAN_SOURCE_COLOUR = "k"
"""Colour to draw a global-mean source over the observational network in

Colour means latitude in that panel, so a source which has no latitude
cannot take a colour from that scale without lying about where it is from.
Black is not on the scale, and it reads over the scatter in greyscale too.
"""


def ghg(pdf: pd.DataFrame, ghg_col: str = "gas") -> str:
    """
    Get GHG of a pandas.DataFrame
    """
    res_l = pdf[ghg_col].unique()
    if len(res_l) != 1:
        raise AssertionError(res_l)

    return res_l[0]


def unit(pdf: pd.DataFrame, unit_col: str = "unit") -> str:
    """
    Get unit of a pandas.DataFrame
    """
    res_l = pdf[unit_col].unique()
    if len(res_l) != 1:
        raise AssertionError(res_l)

    return res_l[0]


GHG_LABELS = {
    "co2": "CO$_2$",
    "ch4": "CH$_4$",
    "n2o": "N$_2$O",
    "c2f6": "C$_2$F$_6$",
    "c3f8": "C$_3$F$_8$",
    "ccl4": "CCl$_4$",
    "cf4": "CF$_4$",
    "cfc11": "CFC-11",
    "cfc113": "CFC-113",
    "cfc114": "CFC-114",
    "cfc115": "CFC-115",
    "cfc12": "CFC-12",
    "ch2cl2": "CH$_2$Cl$_2$",
    "ch3br": "CH$_3$Br",
    "ch3ccl3": "CH$_3$CCl$_3$",
    "ch3cl": "CH$_3$Cl",
    "chcl3": "CHCl$_3$",
    "halon1211": "Halon-1211",
    "halon1301": "Halon-1301",
    "halon2402": "Halon-2402",
    "hcfc141b": "HCFC-141b",
    "hcfc142b": "HCFC-142b",
    "hcfc22": "HCFC-22",
    "hfc125": "HFC-125",
    "hfc134a": "HFC-134a",
    "hfc143a": "HFC-143a",
    "hfc152a": "HFC-152a",
    "hfc227ea": "HFC-227ea",
    "hfc23": "HFC-23",
    "hfc236fa": "HFC-236fa",
    "hfc245fa": "HFC-245fa",
    "hfc32": "HFC-32",
    "hfc365mfc": "HFC-365mfc",
    "hfc4310mee": "HFC-43-10mee",
    "nf3": "NF$_3$",
    "sf6": "SF$_6$",
    "so2f2": "SO$_2$F$_2$",
    "cc4f8": "c-C$_4$F$_8$",
    "c4f10": "C$_4$F$_{10}$",
    "c5f12": "C$_5$F$_{12}$",
    "c6f14": "C$_6$F$_{14}$",
    "c7f16": "C$_7$F$_{16}$",
    "c8f18": "C$_8$F$_{18}$",
}
"""How each gas' name is written when it is shown to a reader

Keyed by the name the gas goes by in the data,
which is what the figures have to hand.
"""


def label_name(ghg: str) -> str:
    """
    Get the label name for a given string

    Parameters
    ----------
    ghg
        String in which to write out any gas names

        This is a string rather than a bare gas name
        because the gas name is usually already part of a label
        by the time we get here.

    Returns
    -------
        `ghg`, with any gas names written out for a reader
    """
    res = ghg
    # Longest name first, so that a gas whose name starts with another gas'
    # name (cfc11 and cfc113, say) is not half-replaced by the shorter one.
    for old, new in sorted(GHG_LABELS.items(), key=lambda kv: -len(kv[0])):
        res = res.replace(old, new)

    return res


def get_only_data_variable(ds: xr.Dataset) -> xr.DataArray:
    """
    Get only data variable in a dataset
    """
    ds_vs = ds.data_vars
    if len(ds_vs) != 1:
        raise AssertionError(ds_vs)

    res = ds_vs[next(iter(ds.data_vars))]

    return res


def get_interpolated_input_coverage_info(
    all_data_with_bins: pd.DataFrame,
    interpolated_obs: xr.Dataset,
) -> dict[str, tuple[int, int]]:
    """
    Get interpolated input coverage information

    Returns the year and month for the time point
    with the most input points and the least input points.
    """
    # Reproduce what was done in the workflow.
    # Break this out if we need it elsewhere.
    interpolated_obs_nan_free = convert_year_month_to_time(
        convert_time_to_year_month(interpolated_obs).dropna("year")
    )
    # The years we can draw, not the span they cover: a gas whose network
    # thinned out for a few years in the middle has years inside that span
    # which were interpolated to nothing, and picking one of those
    # would ask the map panels for a year the interpolated field has not got.
    nan_free_years = set(
        int(year) for year in interpolated_obs_nan_free["time"].dt.year.values
    )
    all_data_with_bins_relevant = all_data_with_bins[
        all_data_with_bins["year"].isin(nan_free_years)
    ].copy()
    all_data_with_bins_relevant["lat-lon-bin"] = (
        all_data_with_bins["lat_bin"].astype(str)
        + "__"
        + all_data_with_bins["lon_bin"].astype(str)
    )
    input_point_counts = all_data_with_bins_relevant.groupby(
        [
            "year",
            "month",
        ]
    )["lat-lon-bin"].count()

    res = {
        "most": input_point_counts.idxmax(),
        "least": input_point_counts.idxmin(),
    }

    return res


def get_decimal_year(indf: pd.DataFrame) -> pd.Series[float]:
    """
    Get the decimal year of each observation

    The data has separate year and month columns rather than a time axis.

    Parameters
    ----------
    indf
        Data for which to calculate the decimal year

    Returns
    -------
        Decimal year of each row of `indf`
    """
    return indf["year"] + (indf["month"] - 0.5) / 12.0


def add_network_group(
    indf: pd.DataFrame, surf_or_ship_col: str = "surf_or_ship"
) -> pd.DataFrame:
    """
    Add the network group of each observation

    Parameters
    ----------
    indf
        Data to which to add the network group

    surf_or_ship_col
        Column which says whether each observation was taken
        at the surface or from a ship

        Not every gas' data has this column.
        Where it is missing, we take every station to be a fixed site.

    Returns
    -------
        `indf`, with the network group added
    """
    out = indf.copy()
    out[NETWORK_GROUP_COLUMN] = out["network"].map(NETWORK_GROUPS)

    unmapped = out[NETWORK_GROUP_COLUMN].isna()
    if unmapped.any():
        msg = f"No network group for {sorted(set(out.loc[unmapped, 'network']))}"
        raise AssertionError(msg)

    if surf_or_ship_col not in out:
        # No shipboard observations for this gas, so every station is fixed
        # and there is nothing to split out.
        return out

    out[NETWORK_GROUP_COLUMN] = out[NETWORK_GROUP_COLUMN] + out[surf_or_ship_col].apply(
        get_moving_suffix
    )

    return pd.concat(
        [
            unify_station_metadata(station_df)
            for _, station_df in out.groupby(["network", "station"])
        ]
    )


def get_moving_suffix(surf_or_ship: str | float) -> str:
    """
    Get the suffix which marks whether a station moves

    Parameters
    ----------
    surf_or_ship
        Whether the observation was taken at the surface or from a ship

        Null means we were not told, in which case we assume a fixed site:
        a moving platform is the thing the data goes out of its way to flag.

    Returns
    -------
        Suffix to add to the network group's name
    """
    if pd.isnull(surf_or_ship):
        return ""

    if surf_or_ship == "surface":
        return ""

    if surf_or_ship == "shipboard":
        return MOVING_SUFFIX

    raise NotImplementedError(surf_or_ship)


def unify_station_metadata(
    station_df: pd.DataFrame, location_tolerance: float = 0.65
) -> pd.DataFrame:
    """
    Give one station one network group and, where we can, one location

    The data reports a network group and a location per observation
    rather than per station, and those can disagree within a station
    (a station which reports from both a fixed site and a ship,
    or a fixed site whose reported position wobbles by a fraction of a degree
    between reports).
    The map draws one marker per station, so it needs one answer for each.

    Parameters
    ----------
    station_df
        One station's observations

    location_tolerance
        How far a station's reported positions may spread, in degrees,
        and still be treated as one position

        Bigger than the wobble in a fixed site's reported position,
        much smaller than the distance a ship covers.

    Returns
    -------
        `station_df`, with its network group, and where possible its location,
        made consistent
    """
    out = station_df.copy()

    group_counts = out[NETWORK_GROUP_COLUMN].value_counts()
    if group_counts.shape[0] > 1:
        # Majority rules
        out[NETWORK_GROUP_COLUMN] = group_counts.index.values[0]

    location_cols = ["latitude", "longitude"]
    locations = out[location_cols]
    if (locations.std() < location_tolerance).all(axis=None):
        for col in location_cols:
            out.loc[:, col] = out[col].median()

    elif (
        locations.drop_duplicates().shape[0] > 1
        and not out[NETWORK_GROUP_COLUMN].str.endswith(MOVING_SUFFIX).all()
    ):
        msg = (
            "A station which does not move reports from more than one place: "
            f"{out['network'].iloc[0]} {out['station'].iloc[0]}"
        )
        raise NotImplementedError(msg)

    return out


def get_network_groups_largest_first(indf: pd.DataFrame) -> list[str]:
    """
    Get the network groups, ordered so the biggest markers are drawn first

    Parameters
    ----------
    indf
        Data from which to get the network groups

    Returns
    -------
        Network groups, largest marker first
    """
    return sorted(
        indf[NETWORK_GROUP_COLUMN].unique(),
        key=lambda group: NETWORK_GROUP_MARKER_SIZES[group],
        reverse=True,
    )


def add_compact_legend(
    ax: matplotlib.axes.Axes, fontsize: str = "x-small", **kwargs: object
) -> None:
    """
    Add a legend which takes as little of its panel as it can

    The panels are small, so a legend at matplotlib's defaults
    is a serious fraction of the panel it sits in,
    and a legend which is taller than its panel spills out over
    the panel's title and its neighbours.

    Parameters
    ----------
    ax
        Axes to add the legend to

    fontsize
        Size to draw the legend's text at

    **kwargs
        Passed on to `ax.legend`
    """
    ax.legend(
        **{
            "fontsize": fontsize,
            "framealpha": 0.9,
            "handletextpad": 0.2,
            "columnspacing": 1.0,
            "labelspacing": 0.3,
            "borderpad": 0.3,
            "borderaxespad": 0.3,
            # Anything the caller asks for wins over what we ask for here
            **kwargs,
        }
    )


def compact_existing_legend(ax: matplotlib.axes.Axes, **kwargs: object) -> None:
    """
    Re-draw a legend which seaborn has already put on a panel, compactly

    Seaborn draws its own legend, at matplotlib's default size,
    so the only way to get a compact one is to ask for it again
    with the handles seaborn worked out.

    Parameters
    ----------
    ax
        Axes whose legend to re-draw

        If it has no legend, nothing happens.

    **kwargs
        Passed on to [add_compact_legend][]
    """
    legend = ax.get_legend()
    if legend is None:
        return

    handles = legend.legend_handles
    labels = [text.get_text() for text in legend.get_texts()]
    title = legend.get_title().get_text()

    add_compact_legend(ax, handles=handles, labels=labels, title=title, **kwargs)
    ax.get_legend().get_title().set_fontsize("x-small")


def eof_palette(eofs: Iterable[int], n_eofs_used: int | None = None) -> dict[int, str]:
    """
    Get the colour to draw each EOF in

    Parameters
    ----------
    eofs
        EOFs to get colours for

    n_eofs_used
        Number of EOFs the workflow keeps

        EOFs at or beyond this are drawn in
        [`UNUSED_EOF_COLOUR`][local.historical_ghg_forcing_for_cmip7.plotting.UNUSED_EOF_COLOUR].
        If `None`, every EOF is treated as kept.

    Returns
    -------
    :
        Colour for each EOF in `eofs`
    """
    res = {}
    for eof in eofs:
        eof = int(eof)  # noqa: PLW2901
        used = n_eofs_used is None or eof < n_eofs_used
        if used and eof < len(EOF_COLOURS):
            res[eof] = EOF_COLOURS[eof]
        else:
            res[eof] = UNUSED_EOF_COLOUR

    return res


def extension_source_palette(sources: Iterable[str]) -> dict[str, str]:
    """
    Get the colour to draw each source of an extended component in

    Parameters
    ----------
    sources
        Sources to get colours for

    Returns
    -------
    :
        Colour for each source in `sources`

    Raises
    ------
    KeyError
        A source has no colour of its own

        Falling back to whatever the plotting library offers is what this
        exists to stop, so a new source has to be given a colour in
        [`EXTENSION_SOURCE_COLOURS`][local.historical_ghg_forcing_for_cmip7.plotting.EXTENSION_SOURCE_COLOURS]
        rather than picking one up by accident.
    """
    sources = list(sources)
    unknown = [source for source in sources if source not in EXTENSION_SOURCE_COLOURS]
    if unknown:
        msg = (
            f"No colour for {unknown}. "
            f"Add one to EXTENSION_SOURCE_COLOURS, which covers "
            f"{sorted(EXTENSION_SOURCE_COLOURS)}."
        )
        raise KeyError(msg)

    return {source: EXTENSION_SOURCE_COLOURS[source] for source in sources}


def region_palette(regions: Iterable[str]) -> dict[str, str]:
    """
    Get the colour to draw each region of a spatial mean in

    Parameters
    ----------
    regions
        Regions to get colours for

    Returns
    -------
    :
        Colour for each region in `regions`

    Raises
    ------
    KeyError
        A region has no colour of its own

        The monthly and the yearly panels are read against each other, so a
        region which picked its colour up from whatever was drawing it could
        be one colour in one and another in the other.
    """
    regions = list(regions)
    unknown = [region for region in regions if region not in REGION_COLOURS]
    if unknown:
        msg = (
            f"No colour for {unknown}. "
            f"Add one to REGION_COLOURS, which covers {sorted(REGION_COLOURS)}."
        )
        raise KeyError(msg)

    return {region: REGION_COLOURS[region] for region in regions}


def source_sequence_palette(labels: Iterable[str]) -> dict[str, str]:
    """
    Get a colour for each of a panel's series, in the order they are given

    Parameters
    ----------
    labels
        Series of the panel, in the order the reader meets them

    Returns
    -------
    :
        Colour for each label

    Raises
    ------
    ValueError
        There are more labels than
        [`SOURCE_SEQUENCE_COLOURS`][local.historical_ghg_forcing_for_cmip7.plotting.SOURCE_SEQUENCE_COLOURS]
        has colours

        Running off the end would mean either repeating a colour inside one
        panel or falling back to a colour from outside the palette, and
        neither is something to do quietly.
    """
    labels = list(labels)
    if len(labels) > len(SOURCE_SEQUENCE_COLOURS):
        msg = (
            f"Only {len(SOURCE_SEQUENCE_COLOURS)} colours, "
            f"but {len(labels)} series to draw: {labels}"
        )
        raise ValueError(msg)

    return dict(zip(labels, SOURCE_SEQUENCE_COLOURS))


def latitude_colour(latitude: float) -> tuple[float, float, float, float]:
    """
    Get the colour which stands for a latitude

    Parameters
    ----------
    latitude
        Latitude to get the colour of

    Returns
    -------
        Colour to draw `latitude` in
    """
    return plt.get_cmap(LATITUDE_COLOUR_MAP)(LATITUDE_NORMALISATION(latitude))


def add_latitude_legend(
    ax: matplotlib.axes.Axes,
    latitudes: np.typing.ArrayLike,
    ncols: int = 2,
    every: int = 2,
    **kwargs: object,
) -> None:
    """
    Add a legend which says which latitude each colour stands for

    Parameters
    ----------
    ax
        Axes to add the legend to

    latitudes
        Latitudes which appear on `ax`

        Listed north first, so the legend runs the same way up as a map does.

    ncols
        Number of columns to lay the legend out in

    every
        Show every nth latitude rather than all of them

        Colour stands for latitude here rather than for a category,
        so the legend is a scale, and a scale only has to be sampled
        finely enough to be read off:
        an entry for every latitudinal bin is more entries
        than a panel this size has room for,
        and says no more than every second one does.

    **kwargs
        Passed on to [add_compact_legend][]
    """
    add_compact_legend(
        ax,
        fontsize="xx-small",
        ncols=ncols,
        handles=[
            matplotlib.lines.Line2D(
                [],
                [],
                linestyle="none",
                marker="o",
                markersize=3,
                color=latitude_colour(latitude),
                label=f"{latitude:.0f}",
            )
            for latitude in sorted(np.asarray(latitudes), reverse=True)[::every]
        ],
        title=r"lat [$^{\circ}$N]",
        **kwargs,
    )
    ax.get_legend().get_title().set_fontsize("xx-small")


def plot_input_timeseries(
    indf: pd.DataFrame,
    ax: matplotlib.axes.Axes,
    lat_column: str = "lat",
    year_column: str = "year",
    value_column: str = "value",
) -> matplotlib.axes.Axes:
    """
    Plot the input timeseries of a gas which has no observational network

    There are only a handful of latitudes here, one line each,
    rather than the thousands of scattered points a network gives,
    so the panel is a line plot and the latitudes go in a legend
    rather than on a colour bar.
    The colours are the same ones the network panels use,
    so a colour means the same latitude everywhere in the figure.

    Parameters
    ----------
    indf
        Input data to plot

    ax
        Axes to plot on

    lat_column
        Column which holds the latitude of each observation

    year_column
        Column which holds the year of each observation

    value_column
        Column which holds the value of each observation

    Returns
    -------
        `ax`
    """
    latitudes = np.sort(indf[lat_column].unique())[::-1]
    for latitude in latitudes:
        at_latitude = indf[indf[lat_column] == latitude].sort_values(year_column)
        ax.plot(
            at_latitude[year_column],
            at_latitude[value_column],
            color=latitude_colour(latitude),
            linewidth=2,
        )

    ax.set_ylabel(f"[{unit(indf)}]", fontsize="small")
    ax.set_xlabel("year", fontsize="small")
    ax.tick_params(labelsize="small")
    # Every latitude, because there are few enough of them to name them all
    add_latitude_legend(ax, latitudes, ncols=1, every=1, loc="upper left")

    return ax


def plot_station_timeseries(  # noqa: PLR0913
    indf: pd.DataFrame,
    ax: matplotlib.axes.Axes,
    # In axis co-ords
    inset_corner: tuple[float, float] | None = None,
    inset_width: float = 0.3,
    inset_height: float = 0.3,
    global_mean_sources: Mapping[str, pd.DataFrame] | None = None,
    robust_y_limits: bool = True,
    clear_inset: bool = True,
) -> matplotlib.collections.PathCollection:
    """
    Plot the monthly mean measured at each station

    This is a scatter plot rather than a line plot,
    so that months in which more than one station reported are visible.
    Colour shows the station's latitude and the marker shows its network,
    so the latitudinal spread of the network can be read directly.

    Parameters
    ----------
    indf
        Observational network data

    ax
        Axes on which to plot

    inset_corner
        Corner to put the inset in, as (x0, y0) in axes co-ordinates

        If `None`, the emptiest of [`INSET_CANDIDATE_CORNERS`][] is used,
        see [choose_inset_corner][].

    inset_width
        Width of the inset, in axes co-ordinates

    inset_height
        Height of the inset, in axes co-ordinates

    global_mean_sources
        Global-mean timeseries which are used in place of,
        or alongside, this network's own global-mean, labelled by source

        For some gases the global-mean does not come from this network at all,
        and a panel which showed only the network would not say
        where the numbers in the rest of the figure actually came from.
        These are global-means, so they have no latitude
        and are drawn as lines rather than on the latitude colour scale.

    robust_y_limits
        Whether to scale the panel by the bulk of the network's values,
        leaving a few far-out points off scale, see [set_robust_y_limits][]

        Turn this off for a gas whose panel should keep every point
        it has, whatever that costs the rest of them.
        Ignored if the caller has already set the panel's vertical limits.

    clear_inset
        Whether to open the panel out until the inset has the corner
        it sits in to itself, see [y_limits_to_clear_inset][]

        On by default: an inset drawn over the data hides measurements,
        and a taller panel only costs empty space.
        Nothing happens where the corner the inset was given is already
        empty, which is most gases.

    Returns
    -------
        The last scatter drawn, so a colour bar can be added for it
    """
    # The panel's vertical limits are settled first, because where there is
    # room for the inset depends on them, and the inset has to exist before
    # anything can be drawn into it.
    panel_values = [indf["value"].to_numpy(dtype=float)]
    if global_mean_sources is not None:
        # The sources are drawn on this panel too, so they get a say in
        # how tall it is, even though the outlier question is about the
        # network: these are global-means, they do not have outliers.
        panel_values.extend(
            in_panel_years(pdf, indf)["value"].to_numpy(dtype=float)
            for pdf in global_mean_sources.values()
        )

    n_off_scale = 0
    if ax.get_autoscaley_on():
        n_off_scale = set_robust_y_limits(
            ax,
            np.concatenate(panel_values),
            headroom_fraction=(
                OUTLIER_HEADROOM_FRACTION if robust_y_limits else np.inf
            ),
        )

    if inset_corner is None:
        inset_corner = choose_inset_corner(
            indf, ax.get_ylim(), inset_width, inset_height
        )

    inset_x0, inset_y0 = inset_corner

    inset_xlim = (indf["year"].max() - 3, indf["year"].max())
    inset_values = indf[
        (indf["year"] >= inset_xlim[0]) & (indf["year"] <= inset_xlim[1])
    ]
    inset_values_range = inset_values["value"].max() - inset_values["value"].min()
    inset_ylim = (
        np.floor(max(inset_values["value"].min() - 0.05 * inset_values_range, 0.0)),
        np.ceil(inset_values["value"].max() + 0.05 * inset_values_range),
    )
    ax_inset = ax.inset_axes(
        [inset_x0, inset_y0, inset_width, inset_height],
        xlim=inset_xlim,
        ylim=inset_ylim,
        # xticklabels=[],
        # yticklabels=[],
    )
    ax_inset.tick_params(labelsize="small")
    ax_inset.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))

    scatter = None
    for group, group_df in indf.groupby(NETWORK_GROUP_COLUMN):
        for axh in [ax, ax_inset]:
            scatter = axh.scatter(
                get_decimal_year(group_df),
                group_df["value"],
                c=group_df["latitude"],
                cmap=LATITUDE_COLOUR_MAP,
                norm=LATITUDE_NORMALISATION,
                marker=NETWORK_GROUP_MARKERS[group],
                s=10.0,
                alpha=0.4,
                linewidths=0.0,
            )

    if scatter is None:
        msg = "No data to plot"
        raise AssertionError(msg)

    # Colour is doing latitude here, so the legend can only speak about shape
    handles = [
        matplotlib.lines.Line2D(
            [],
            [],
            linestyle="none",
            marker=NETWORK_GROUP_MARKERS[group],
            color=LEGEND_MARKER_COLOUR,
            markersize=5,
            label=group,
        )
        for group in sorted(indf[NETWORK_GROUP_COLUMN].unique())
    ]

    if global_mean_sources is not None:
        for i, (label, pdf) in enumerate(global_mean_sources.items()):
            in_panel = in_panel_years(pdf, indf)
            # Drawn on both axes so the inset tells the same story as the panel
            for axh in [ax, ax_inset]:
                axh.plot(
                    in_panel["year"],
                    in_panel["value"],
                    color=GLOBAL_MEAN_SOURCE_COLOUR,
                    linestyle=GLOBAL_MEAN_SOURCE_LINESTYLES[
                        i % len(GLOBAL_MEAN_SOURCE_LINESTYLES)
                    ],
                    linewidth=1.5,
                    label=label,
                )

            handles.append(
                matplotlib.lines.Line2D(
                    [],
                    [],
                    color=GLOBAL_MEAN_SOURCE_COLOUR,
                    linestyle=GLOBAL_MEAN_SOURCE_LINESTYLES[
                        i % len(GLOBAL_MEAN_SOURCE_LINESTYLES)
                    ],
                    linewidth=1.5,
                    label=label,
                )
            )

    if n_off_scale:
        # In the legend rather than on the panel: it is a note about what
        # is drawn, and the legend is the one place on a panel this full
        # which is guaranteed not to have a measurement under it.
        handles.append(
            matplotlib.lines.Line2D(
                [],
                [],
                linestyle="none",
                marker="none",
                label=(
                    f"{n_off_scale} point{'' if n_off_scale == 1 else 's'} off scale"
                ),
            )
        )

    if clear_inset:
        ax.set_ylim(
            *y_limits_to_clear_inset(
                indf, ax.get_ylim(), inset_corner, inset_width, inset_height
            )
        )

    # The inset is a zoom of this panel, so it cannot show what the panel
    # has put off scale: the box drawn around the years it covers would run
    # off the top of the panel it is a zoom of. Drawn only now, because the
    # box has to be drawn against the limits both axes have ended up with.
    panel_lower, panel_upper = ax.get_ylim()
    inset_lower, inset_upper = ax_inset.get_ylim()
    ax_inset.set_ylim(max(inset_lower, panel_lower), min(inset_upper, panel_upper))
    ax.indicate_inset_zoom(ax_inset, edgecolor="black")

    units = unit(indf)
    # gas = ghg(indf)

    ax.set_ylabel(label_name(f"[{units}]"), fontsize="small")
    ax.set_xlabel("year", fontsize="small")
    ax.tick_params(labelsize="small")

    add_compact_legend(ax, loc="upper left", handles=handles)

    return scatter


def in_panel_years(
    pdf: pd.DataFrame, indf: pd.DataFrame, year_column: str = "year"
) -> pd.DataFrame:
    """
    Cut a source down to the years the observational network panel shows

    These sources usually run well past the network in both directions
    (Trudinger starts in 1901, the network in 2003), and the panel is held
    to the network's years, so the rest is never seen. Keeping it would
    leave it in the sums matplotlib scales the vertical axis by, which
    squashes every point in the panel into a corner to make room for a line
    nobody can see.

    Parameters
    ----------
    pdf
        Source to cut down

    indf
        Observational network data, whose years the panel shows

    year_column
        Column which holds the year

    Returns
    -------
        `pdf`, with the years outside the network's own dropped
    """
    return pdf[
        (pdf[year_column] >= indf[year_column].min())
        & (pdf[year_column] <= indf[year_column].max())
    ]


def y_limits_to_clear_inset(  # noqa: PLR0913
    indf: pd.DataFrame,
    y_limits: tuple[float, float],
    corner: tuple[float, float],
    inset_width: float,
    inset_height: float,
    clearance: float = INSET_CLEARANCE,
    value_column: str = "value",
) -> tuple[float, float]:
    """
    Get the limits a panel needs for its inset to sit clear of the data

    The inset is parked in a corner of the panel it magnifies,
    which only works while that corner is empty.
    Where it is not empty, the panel has to open up to make it so:
    the floor drops for an inset low in the panel,
    the ceiling lifts for one high in it.

    Parameters
    ----------
    indf
        Observational network data drawn on the panel

    y_limits
        Vertical limits the panel has

    corner
        Corner the inset sits in, as (x0, y0) in axes co-ordinates

    inset_width
        Width of the inset, in axes co-ordinates

    inset_height
        Height of the inset, in axes co-ordinates

    clearance
        Room to leave between the inset and the data,
        as a fraction of the panel's height

    value_column
        Column of `indf` which holds the measured values

    Returns
    -------
        Limits which leave the corner to the inset

        `y_limits` unchanged, if the corner is already clear.
    """
    inset_x0, inset_y0 = corner
    lower, upper = y_limits

    # The panel shows the network's own span, give or take the year either
    # side the callers add, so a share of the panel's width is a share of that
    decimal_year = get_decimal_year(indf).to_numpy(dtype=float)
    first_year = decimal_year.min() - 1.0
    last_year = decimal_year.max() + 1.0
    x_fraction = (decimal_year - first_year) / (last_year - first_year)
    under_inset = indf[
        (x_fraction >= inset_x0) & (x_fraction <= inset_x0 + inset_width)
    ][value_column].to_numpy(dtype=float)
    under_inset = under_inset[np.isfinite(under_inset)]
    if under_inset.size < 1:
        return y_limits

    # Which half of the panel the inset's own middle falls in
    panel_middle = 0.5
    if inset_y0 + inset_height / 2.0 < panel_middle:
        # Low in the panel: the data has to end up above the inset's top
        inset_top = inset_y0 + inset_height + clearance
        if inset_top >= 1.0:
            return y_limits

        needed_lower = (under_inset.min() - inset_top * upper) / (1.0 - inset_top)

        # Only ever outwards: the inset needs room, it does not need company
        return (min(lower, needed_lower), upper)

    # High in the panel: the data has to end up below the inset's bottom
    inset_bottom = inset_y0 - clearance
    if inset_bottom <= 0.0:
        return y_limits

    needed_upper = (under_inset.max() - lower * (1.0 - inset_bottom)) / inset_bottom

    return (lower, max(upper, needed_upper))


def choose_inset_corner(  # noqa: PLR0913
    indf: pd.DataFrame,
    y_limits: tuple[float, float],
    inset_width: float,
    inset_height: float,
    candidates: tuple[tuple[float, float], ...] = INSET_CANDIDATE_CORNERS,
    clearance: float = INSET_CLEARANCE,
) -> tuple[float, float]:
    """
    Pick the corner of a panel which costs the least to give to the inset

    The corner the inset goes in should be one the panel is not using,
    and which corner that is depends on the gas: the ones which are still
    accumulating leave the bottom of their panel empty, and the ones which
    have been phased out leave the top empty.

    A corner which is already empty costs nothing, so those win outright.
    Where no corner is empty, the panel has to open up to clear one,
    and the cheapest corner is the one which needs it to open up least.
    Scoring the corners this way rather than by how many points are in them
    is the difference between paying a little dead space and paying a lot:
    a corner with a handful of points hard against it can cost more room
    than a corner with many points which barely reach in.

    Parameters
    ----------
    indf
        Observational network data which will be drawn on the panel

    y_limits
        Vertical limits the panel would have if the inset needed nothing

    inset_width
        Width of the inset, in axes co-ordinates

    inset_height
        Height of the inset, in axes co-ordinates

    candidates
        Corners to choose between, as (x0, y0) in axes co-ordinates

        Ties go to the earliest in this order.

    clearance
        Room to leave between the inset and the data,
        as a fraction of the panel's height

    Returns
    -------
        The chosen corner, as (x0, y0) in axes co-ordinates
    """

    def cost(corner: tuple[float, float]) -> float:
        lower, upper = y_limits_to_clear_inset(
            indf, y_limits, corner, inset_width, inset_height, clearance=clearance
        )

        return (upper - lower) / (y_limits[1] - y_limits[0])

    return min(candidates, key=cost)


def set_robust_y_limits(
    ax: matplotlib.axes.Axes,
    values: np.typing.ArrayLike,
    percentiles: tuple[float, float] = ROBUST_Y_LIMIT_PERCENTILES,
    headroom_fraction: float = OUTLIER_HEADROOM_FRACTION,
    pad_fraction: float = 0.05,
) -> int:
    """
    Scale a panel by the bulk of its values rather than by its extremes

    A handful of outlying observations can take most of a panel's height
    for themselves, leaving everything else squashed into a band too thin
    to read. Where that is happening, this scales the panel to the values
    which are not outliers and reports how many are left outside,
    so the panel can say so rather than quietly dropping them.

    Parameters
    ----------
    ax
        Axes to set the limits of

    values
        Values which are drawn on `ax`

    percentiles
        Percentiles of `values` to treat as the robust range

    headroom_fraction
        How far past the robust range a value has to be to be left off scale,
        as a fraction of the robust range

    pad_fraction
        Room to leave above and below, as a fraction of the range shown

    Returns
    -------
        Number of `values` left outside the limits
    """
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size < 1:
        return 0

    lower_percentile, upper_percentile = np.percentile(finite, percentiles)
    robust_span = upper_percentile - lower_percentile
    lower = finite.min()
    upper = finite.max()
    if robust_span <= 0.0:
        return 0

    # Only where the extremes are far enough out to be costing the panel:
    # a gas whose values simply spread out keeps every one of them.
    if upper - upper_percentile > headroom_fraction * robust_span:
        upper = upper_percentile

    if lower_percentile - lower > headroom_fraction * robust_span:
        lower = lower_percentile

    pad = pad_fraction * (upper - lower)
    ax.set_ylim(lower - pad, upper + pad)

    return int(((finite < lower - pad) | (finite > upper + pad)).sum())


def plot_station_locations(
    indf: pd.DataFrame,
    ax: matplotlib.axes.Axes,
) -> None:
    """
    Plot the location of each station in the observational network

    Parameters
    ----------
    indf
        Observational network data

    ax
        Axes on which to plot

        This must have been created with a cartopy projection.
    """
    ax.coastlines(linewidth=0.4, color="0.55")
    ax.set_global()

    # The bins the observations are binned into,
    # so this panel can be read against the counts panel
    for lat_bound in LAT_BIN_BOUNDS:
        ax.axhline(lat_bound, linewidth=0.4, color="0.85", zorder=0)

    for lon_bound in LON_BIN_BOUNDS:
        ax.axvline(lon_bound, linewidth=0.4, color="0.85", zorder=0)

    stations = indf[
        [NETWORK_GROUP_COLUMN, "station", "latitude", "longitude"]
    ].drop_duplicates()

    # Biggest markers first, so co-located networks nest rather than hide each other
    for group in get_network_groups_largest_first(stations):
        group_stations = stations[stations[NETWORK_GROUP_COLUMN] == group]
        ax.scatter(
            group_stations["longitude"],
            group_stations["latitude"],
            transform=ccrs.PlateCarree(),
            facecolors="none",
            edgecolors=NETWORK_GROUP_COLOURS[group],
            marker=NETWORK_GROUP_MARKERS[group],
            s=NETWORK_GROUP_MARKER_SIZES[group],
            linewidths=1.1,
            label=group,
            zorder=3,
            # Cartopy clips to the projection boundary,
            # which would cut the South Pole station in half
            clip_on=False,
        )

    ax.set_yticks(LAT_BIN_BOUNDS[::2], crs=ccrs.PlateCarree())
    ax.set_ylabel(r"latitude [$^{\circ}$N]", fontsize="small")
    ax.set_xticks(LON_BIN_BOUNDS[::2], crs=ccrs.PlateCarree())
    ax.set_xlabel(r"longitude [$^{\circ}$E]", fontsize="small")
    # The map's aspect is fixed, so it can't fill its cell.
    # Anchoring it to the top keeps it up against its panel label.
    ax.set_anchor("N")
    ax.tick_params(labelsize="small")
    ax.set_ylim(LAT_AXIS_LIMITS)
    # On the panel's title line, which has room to spare and is the only
    # place around this panel that does: the map cannot be drawn over
    # (that would cover stations), its row is only as tall as the map itself
    # so there is no band under it, and the gap beside it is a gap between
    # two columns rather than space belonging to this panel.
    add_compact_legend(
        ax,
        fontsize="xx-small",
        loc="lower right",
        bbox_to_anchor=(1.0, 1.0),
        # All on one row, so the legend is no taller than the title line it shares.
        # The maps are wide enough that it still stays clear of the title.
        ncols=len(get_network_groups_largest_first(stations)),
        handlelength=1.0,
        handletextpad=0.3,
    )
    # The gap beside the map is there whether or not we put the legend in it,
    # so the legend is placed by hand and kept out of the layout engine's sums.
    # Left in them, it would ask for its width from the next panel's column,
    # which every panel in that column would then pay for.
    ax.get_legend().set_in_layout(False)


def plot_observation_counts(
    indf: pd.DataFrame,
    ax: matplotlib.axes.Axes,
) -> matplotlib.collections.QuadMesh:
    """
    Plot the number of input data points in each month-latitudinal bin

    The networks report monthly averages, so each row of the data
    is one station's average for one month, and that is what is counted here.
    Many more measurements sit underneath each of those averages,
    but the networks do not report them, so we cannot count them.
    Counts are summed over all longitudes.

    Parameters
    ----------
    indf
        Observational network data

    ax
        Axes on which to plot

    Returns
    -------
        The mesh which was drawn, so a colour bar can be added for it
    """
    # The bins the observations are binned into
    for lat_bound in LAT_BIN_BOUNDS:
        ax.axhline(lat_bound, linewidth=0.4, color="0.85", zorder=3)

    counts = (
        indf.groupby(["year", "month", "lat_bin"]).size().rename("count").reset_index()
    )
    counts["decimal_year"] = get_decimal_year(counts)

    grid = counts.pivot_table(
        index="lat_bin", columns="decimal_year", values="count", fill_value=0
    )
    # Make sure every latitudinal bin has a row, even the empty ones
    lat_bin_centres = (LAT_BIN_BOUNDS[:-1] + LAT_BIN_BOUNDS[1:]) / 2
    grid = grid.reindex(lat_bin_centres, fill_value=0)

    # Every month in the record, so gaps show up as gaps rather than being closed up
    month_starts = np.arange(
        np.floor(counts["decimal_year"].min()),
        np.ceil(counts["decimal_year"].max()) - 1 / 24,
        1 / 12,
    )
    grid = grid.reindex(
        columns=grid.columns.union(month_starts + 1 / 24), fill_value=0
    ).sort_index(axis="columns")

    x_bounds = np.append(grid.columns.to_numpy() - 1 / 24, grid.columns[-1] + 1 / 24)

    # Counts are small integers, so use a discrete scale rather than a continuous one
    max_count = int(grid.to_numpy().max())
    mesh = ax.pcolormesh(
        x_bounds,
        LAT_BIN_BOUNDS,
        np.ma.masked_equal(grid.to_numpy(), 0),
        cmap=plt.get_cmap("YlOrRd", max_count),
        norm=matplotlib.colors.BoundaryNorm(np.arange(0.5, max_count + 1.0), max_count),
        shading="flat",
    )
    ax.set_yticks(LAT_BIN_BOUNDS[::2])
    ax.set_ylim(LAT_AXIS_LIMITS)
    ax.set_ylabel(r"latitude [$^{\circ}$N]", fontsize="small")
    ax.tick_params(labelsize="small")

    return mesh


def plot_coverage_and_interpolated(
    input_data: pd.DataFrame,
    interpolated: xr.Dataset,
    year_month: tuple[int, int],
    ax: matplotlib.axes.Axes,
    # return type hint is wrong
) -> None:
    """
    Plot the interpolated values and the coverage of input data
    """
    ax.coastlines(linewidth=0.6, color="0.3", zorder=2.0)
    ax.set_global()

    for lat_bound in LAT_BIN_BOUNDS:
        ax.axhline(lat_bound, linewidth=0.4, color="0.85", zorder=0)

    for lon_bound in LON_BIN_BOUNDS:
        ax.axvline(lon_bound, linewidth=0.4, color="0.85", zorder=0)

    lon_grid, lat_grid = np.meshgrid(
        LON_BIN_CENTRES,
        LAT_BIN_CENTRES,
    )

    interpolated_ym = convert_time_to_year_month(interpolated).sel(
        year=year_month[0], month=year_month[1]
    )
    interpolated_ym_da = get_only_data_variable(interpolated_ym)
    # interpolated_ym_vs = interpolated_ym.data_vars
    # if len(interpolated_ym_vs) != 1:
    #     raise AssertionError
    #
    # interpolated_ym_da = interpolated_ym_vs[next(iter(interpolated_ym.data_vars))]

    mesh = ax.pcolormesh(
        lon_grid,
        lat_grid,
        interpolated_ym_da.T,
        shading="auto",
        cmap=INTERPOLATION_COLOUR_MAP,
    )

    input_data_ym = input_data[
        (input_data["year"] == year_month[0]) & (input_data["month"] == year_month[1])
    ]
    ax.scatter(
        input_data_ym["lon_bin"],
        input_data_ym["lat_bin"],
        transform=ccrs.PlateCarree(),
        c="k",
        marker="o",
        s=10.0,
        label="Input point",
        zorder=3,
    )

    ax.set_yticks(LAT_BIN_BOUNDS[::2], crs=ccrs.PlateCarree())
    ax.set_ylabel(r"latitude [$^{\circ}$N]", fontsize="small")
    ax.set_xticks(LON_BIN_BOUNDS[::2], crs=ccrs.PlateCarree())
    ax.set_xlabel(r"longitude [$^{\circ}$E]", fontsize="small")
    # The map's aspect is fixed, so it can't fill its cell.
    # Anchoring it to the top keeps it up against its panel label.
    ax.set_anchor("N")
    ax.tick_params(labelsize="small")
    # Same as the locations map, so all the maps are the same shape
    ax.set_ylim(LAT_AXIS_LIMITS)

    return mesh


def plot_global_mean_from_obs_network(gm: xr.Dataset, ax: matplotlib.axes.Axes) -> None:
    """
    Plot global-mean derived from the observational network
    """
    gm_da = get_only_data_variable(gm)
    ax.scatter(
        gm_da["year"].values.squeeze(),
        gm_da.values.squeeze(),
        # The same colour the extended global-mean takes in its own panel,
        # which sits in the same row: one quantity, one colour.
        color=SOURCE_SEQUENCE_COLOURS[0],
        s=15,
    )
    ax.set_ylabel(f"[{gm_da.attrs['units']}]", fontsize="small")
    ax.tick_params(labelsize="small")


def plot_seasonality_from_obs_network(
    seasonality: xr.Dataset,
    ax: matplotlib.axes.Axes,
    assumed_units: str,
    # assumed_ghg: str,
) -> matplotlib.axes.Axes:
    """
    Plot seasonality derived from the observational network

    There is one series here for each of the twelve latitudinal bins.
    Latitude is shown with the same colour map, over the same range,
    as the timeseries panel uses, so a colour means the same latitude
    everywhere in the figure.
    """
    seasonality_da = get_only_data_variable(seasonality)
    pdf = seasonality_da.to_pandas().stack().rename("value").to_frame().reset_index()
    ax.scatter(
        pdf["month"],
        pdf["value"],
        c=pdf["lat"],
        cmap=LATITUDE_COLOUR_MAP,
        norm=LATITUDE_NORMALISATION,
        s=30.0,
        linewidths=0.0,
    )
    ax.set_ylabel(f"[{seasonality_da.attrs['units']}]", fontsize="small")
    ax.set_xlabel("month", fontsize="small")
    ax.set_xticks(np.arange(1, 12 + 1, 3))
    ax.tick_params(labelsize="small")
    add_latitude_legend(ax, pdf["lat"].unique(), loc="best")

    return ax


def linear_latitudinal_gradient_eof(units: str) -> xr.DataArray:
    """
    Get the latitudinal gradient EOF assumed for gases with no observational network

    The original run assumes a single EOF for these gases,
    linear in latitude before the area-weight of each latitudinal box
    is taken into account, normalised so that its area-weighted mean is zero
    and its northern (southern) hemispheric mean is a half (minus a half).
    That makes the PC which goes with it
    the difference between the hemispheric means.

    This mirrors the EOF built in the original run's
    `calculate_c4f10_like_monthly_fifteen_degree_pieces`
    and `calculate_c8f18_like_monthly_fifteen_degree_pieces` steps.
    The C4F10-like gases have theirs written into a data file,
    so this is really only needed for C8F18, whose EOF was never saved out;
    the C4F10-like figure checks its stored EOF against this one,
    which is what keeps the two honest.

    Parameters
    ----------
    units
        Units to give the EOF

        Which way round the units fall between the EOF and its PC
        depends on the gas, so the caller has to say.

    Returns
    -------
        The assumed latitudinal gradient EOF, by latitude

    Raises
    ------
    AssertionError
        The EOF does not have the properties it is normalised to have
    """
    lat_bin_weights = np.diff(np.sin(LAT_BIN_BOUNDS * np.pi / 180.0))
    n_bins = len(LAT_BIN_CENTRES)
    eof = np.linspace(-1.0, 1.0, n_bins) / lat_bin_weights

    northern = slice(n_bins // 2, n_bins)
    southern = slice(0, n_bins // 2)

    def hemispheric_mean(half: slice) -> float:
        return float(
            np.sum(lat_bin_weights[half] * eof[half]) / np.sum(lat_bin_weights[half])
        )

    eof = eof / (2.0 * hemispheric_mean(northern))

    # The same checks the original run made, for the same reason:
    # the normalisation is what lets the PC be read as a hemispheric difference,
    # so it is worth knowing straight away if it has gone wrong.
    if not np.isclose((lat_bin_weights * eof).sum(), 0.0):
        msg = "The spatial mean of the latitudinal gradient EOF is not zero"
        raise AssertionError(msg)

    for half, half_name, expected in (
        (northern, "northern", 0.5),
        (southern, "southern", -0.5),
    ):
        if not np.isclose(hemispheric_mean(half), expected):
            msg = (
                f"The {half_name} hemispheric mean of the latitudinal gradient EOF "
                f"is {hemispheric_mean(half)}, expected {expected}"
            )
            raise AssertionError(msg)

    return xr.DataArray(
        eof,
        dims=("lat",),
        coords={"lat": LAT_BIN_CENTRES},
        attrs={"units": units},
    )


def plot_lat_gradient_pieces_from_obs_network(
    lat_gradient_info: xr.Dataset,
    axes: Mapping[str, matplotlib.axes.Axes],
    pcs_name: str = "principal-components",
    eofs_name: str = "eofs",
) -> dict[str, matplotlib.axes.Axes]:
    """
    Plot seasonality derived from the observational network
    """
    da_pcs = lat_gradient_info[pcs_name]
    pdf_pcs = da_pcs.to_pandas().stack().rename("value").to_frame().reset_index()
    sns.scatterplot(
        pdf_pcs,
        x="year",
        y="value",
        hue="eof",
        palette=eof_palette(pdf_pcs["eof"].unique()),
        ax=axes["pcs"],
    )
    axes["pcs"].set_ylabel(f"[{da_pcs.attrs['units']}]", fontsize="small")
    axes["pcs"].set_xlabel("year", fontsize="small")
    axes["pcs"].tick_params(labelsize="small")
    compact_existing_legend(axes["pcs"], loc="best")

    plot_lat_gradient_eofs(lat_gradient_info[eofs_name], axes["eofs"])

    return axes


def plot_lat_gradient_eofs(
    da_eofs: xr.DataArray, ax: matplotlib.axes.Axes
) -> matplotlib.axes.Axes:
    """
    Plot the latitudinal gradient's EOFs

    Parameters
    ----------
    da_eofs
        EOFs to plot, with an `eof` and a `lat` dimension

    ax
        Axes to plot on

    Returns
    -------
        `ax`
    """
    pdf_eofs = da_eofs.to_pandas().stack().rename("value").to_frame().reset_index()
    sns.scatterplot(
        pdf_eofs,
        x="value",
        y="lat",
        hue="eof",
        palette=eof_palette(pdf_eofs["eof"].unique()),
        ax=ax,
    )
    ax.set_yticks(LAT_BIN_BOUNDS[::2])
    ax.set_ylabel(r"latitude [$^{\circ}$N]", fontsize="small")
    ax.set_xlabel(f"[{da_eofs.attrs['units']}]", fontsize="small")
    ax.tick_params(labelsize="small")
    compact_existing_legend(ax, loc="best")

    return ax


def plot_variance_explained(  # noqa: PLR0913
    variance_explained: pd.DataFrame,
    ax: matplotlib.axes.Axes,
    n_eofs_used: int,
    eof_col: str = "eof",
    value_col: str = "variance_explained_fraction",
    cumulative_colour: str = "#3a3b3a",
) -> matplotlib.axes.Axes:
    """
    Plot how much of the variance each EOF explains

    The bars are what each EOF explains on its own, the line is the running
    total. Together they say how much is lost by keeping only the first
    `n_eofs_used` EOFs, which is the choice this panel is here to justify.

    Both are drawn against the same vertical axis. A running total deserves
    its own axis when it would otherwise be squashed, but here the first EOF
    explains almost everything, so the total starts at the top of the first
    bar and has nowhere to go but the last few percent to 100 --
    which is the panel's whole point, and is lost if the two are rescaled
    against each other.

    This is deliberately agnostic about what was decomposed,
    so it serves the latitudinal gradient and the seasonality change alike.

    Parameters
    ----------
    variance_explained
        Fraction of the variance explained by each EOF

    ax
        Axes to plot on

    n_eofs_used
        Number of EOFs the workflow keeps

        These are drawn in their own colour, the rest in
        [`UNUSED_EOF_COLOUR`][local.historical_ghg_forcing_for_cmip7.plotting.UNUSED_EOF_COLOUR].

    eof_col
        Column of `variance_explained` which holds the EOF

    value_col
        Column of `variance_explained` which holds the fraction explained

    cumulative_colour
        Colour to draw the running total in

    Returns
    -------
    :
        `ax`
    """
    pdf = variance_explained.sort_values(eof_col)
    eofs = pdf[eof_col].to_numpy()
    percentages = 100.0 * pdf[value_col].to_numpy()

    palette = eof_palette(eofs, n_eofs_used=n_eofs_used)
    ax.bar(
        eofs,
        percentages,
        color=[palette[int(eof)] for eof in eofs],
        width=0.7,
        label="Individual",
    )
    ax.plot(
        eofs,
        np.cumsum(percentages),
        color=cumulative_colour,
        marker="o",
        markersize=3.0,
        linewidth=1.0,
        label="Cumulative",
    )

    ax.set_xlabel("EOF", fontsize="small")
    ax.set_ylabel("[%]", fontsize="small")
    ax.set_xticks(thin_ticks(eofs, max_ticks=6))
    ax.set_ylim(0.0, 105.0)
    ax.tick_params(labelsize="small")
    # Named in the order the reader meets them, which is not the order
    # matplotlib collects them in: it groups the artists by type,
    # which puts the line before the bars it is the running total of.
    handles = {
        handle.get_label(): handle for handle in ax.get_legend_handles_labels()[0]
    }
    labels = ["Individual", "Cumulative"]
    add_compact_legend(
        ax,
        handles=[handles[label] for label in labels],
        labels=labels,
        loc="center right",
    )

    return ax


def add_break_lines_and_setup(  # noqa: PLR0913
    ax_left: matplotlib.axes.Axes,
    ax_right: matplotlib.axes.Axes,
    min_year: int,
    split_year: int,
    max_year: int,
    units: str,
    legend_loc: str = "best",
    max_ticks_per_half: int | None = None,
) -> None:
    """
    Add break lines to axes and do other general setup for broken axes

    Parameters
    ----------
    ax_left
        Left half of the broken axis

    ax_right
        Right half of the broken axis

    min_year
        First year the left half covers

    split_year
        Year the axis is broken at

    max_year
        Last year the right half covers

    units
        Units of the values plotted

    legend_loc
        Where to put the pair's legend

        Passed on to `ax.legend`, so `"best"` leaves it to matplotlib.
        Worth naming where matplotlib's answer is not a good one:
        it places the legend where the data is thinnest,
        which on a panel whose data hugs one edge
        is on top of the axis' own tick labels.

        The legend goes on the left half, which is only part of a panel wide,
        so a legend with long labels in it is wider than the axes it sits in.
        `"best"` centres such a legend, which leaves it hanging out of
        both sides of the half and over the neighbouring panel's tick labels.
        Anchoring it to a corner instead keeps the overhang on one side,
        over the panel's own other half.

    max_ticks_per_half
        The most x tick labels to leave on each half

        Matplotlib picks ticks for the range it is given, not for the room
        it has to show them in, so a narrow half ends up with its year labels
        run together. If `None`, matplotlib's answer is kept.
    """
    # One legend for the pair, on the left half:
    # the two halves are one panel as far as a reader is concerned.
    if ax_right.get_legend() is not None:
        ax_right.get_legend().remove()

    if ax_left.get_legend() is None:
        add_compact_legend(ax_left, loc=legend_loc)
    else:
        compact_existing_legend(ax_left, loc=legend_loc)

    # A legend which is wider than the half it sits in
    # has to draw over the other half to be read at all.
    # Axes are drawn in the order they were added,
    # so without this the right half's background
    # paints over whatever hangs into it,
    # which cuts the legend off mid-label.
    ax_left.set_zorder(ax_right.get_zorder() + 1)
    ax_left.patch.set_visible(False)

    ax_left.set_ylabel(f"[{units}]", fontsize="small")
    ax_right.set_ylabel("")
    for ax in (ax_left, ax_right):
        ax.set_xlabel("year", fontsize="small")
        ax.tick_params(labelsize="small")
    # ax.set_xlabel("year")
    ax_left.set_xlim(xmin=min_year, xmax=split_year)
    ax_right.set_xlim(xmin=split_year, xmax=max_year)

    # Hide the spines
    ax_left.spines["right"].set_visible(False)
    ax_right.spines["left"].set_visible(False)

    # The two halves are one panel, so they share one vertical scale,
    # and the right half's copy of it would say nothing the left half
    # has not already said. Hiding it is not only tidier:
    # a panel's tick labels are the most expensive thing about it
    # as far as the layout engine is concerned,
    # and the halves are the only panels which can give that cost up.
    ax_right.set_ylim(ax_left.get_ylim())
    ax_right.yaxis.set_visible(False)

    # Potential option for cutout lines.
    # The important
    # thing to know here is that in axes coordinates, which are always
    # between 0-1, spine endpoints are at these locations (0, 0), (0, 1),
    # (1, 0), and (1, 1).  Thus, we just need to put the diagonals in the
    # appropriate corners of each of our axes, and so long as we use the
    # right transform and disable clipping.

    clear_ticks_near_break(ax_left, at="right")
    clear_ticks_near_break(ax_right, at="left")

    if max_ticks_per_half is not None:
        # After clearing, not before: the ticks near the break go whatever
        # happens, so thinning first would spend part of the budget on ticks
        # which are about to be dropped anyway.
        for ax in (ax_left, ax_right):
            x_limits = ax.get_xlim()
            ax.set_xticks(thin_ticks(ax.get_xticks(), max_ticks=max_ticks_per_half))
            ax.set_xlim(x_limits)

    d = 0.015  # how big to make the diagonal lines in axes coordinates
    # arguments to pass plot, just so we don't keep repeating them
    kwargs = dict(transform=ax_left.transAxes, color="k", clip_on=False)
    ax_left.plot((1 - d, 1 + d), (-d, +d), **kwargs)
    ax_left.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)

    kwargs.update(transform=ax_right.transAxes)  # switch to the bottom axes
    ax_right.plot((-d, +d), (1 - d, 1 + d), **kwargs)
    ax_right.plot((-d, +d), (-d, +d), **kwargs)


def clear_ticks_near_break(
    ax: matplotlib.axes.Axes, at: str, keep_clear: float = 0.06
) -> None:
    """
    Drop the tick labels which sit right against a broken axis' break

    The two halves of a broken axis are drawn hard up against each other,
    so a tick label at the very end of one half
    lands on top of a tick label at the very start of the other.

    Parameters
    ----------
    ax
        Half of a broken axis

    at
        Which end of `ax` the break is at, `"left"` or `"right"`

    keep_clear
        How much of the axis to keep clear of tick labels,
        as a fraction of the axis' range

    Raises
    ------
    ValueError
        `at` is neither `"left"` nor `"right"`
    """
    if at not in ("left", "right"):
        msg = f"at must be 'left' or 'right', got {at!r}"
        raise ValueError(msg)

    x_min, x_max = ax.get_xlim()
    keep_clear_from = keep_clear * (x_max - x_min)

    if at == "right":
        keep = [tick for tick in ax.get_xticks() if tick < x_max - keep_clear_from]
    else:
        keep = [tick for tick in ax.get_xticks() if tick > x_min + keep_clear_from]

    ax.set_xticks(keep)
    # set_xticks resets the limits it was given the ticks for
    ax.set_xlim(x_min, x_max)


def plot_global_mean_extension(  # noqa: PLR0913
    gm: xr.Dataset,
    ax_left: matplotlib.axes.Axes,
    ax_right: matplotlib.axes.Axes,
    input_sources=Mapping[str, pd.DataFrame],
    split_year: int = 1950,
    legend_loc: str = "upper left",
    pre_industrial: tuple[int, float] | None = None,
    fit_period: tuple[int, int] | None = None,
) -> None:
    """
    Plot global-mean derived from the observational network

    Parameters
    ----------
    gm
        Extended global-mean

    ax_left
        Left half of the broken axis

    ax_right
        Right half of the broken axis

    input_sources
        Sources which went into the extension, labelled by source

    split_year
        Year the axis is broken at

    legend_loc
        Where to put the panel's legend

    pre_industrial
        Year and value of the pre-industrial point the extension is anchored to

        The extension is flat at this value for every year up to this one,
        so it is the left-hand end of everything the panel shows.

    fit_period
        First and last year filled by fitting between the pre-industrial point
        and the first year the data sources cover

        These years are neither observed nor taken from a reference source,
        which is worth marking on a panel where they are otherwise
        indistinguishable from the years which are.
        They are drawn as their own stretch of the line,
        as are the years before them, which are simply assumed
        to sit at the pre-industrial value.
    """
    gm_da = get_only_data_variable(gm)
    gm_years = gm_da["year"].values.squeeze()
    gm_values = gm_da.values.squeeze()

    # The extended global-mean is the panel's subject and the sources are what
    # it was built from, so it takes the first colour and they follow it.
    extended_label = "Extended global-mean"
    palette = source_sequence_palette([extended_label, *input_sources])

    for i, ax in enumerate((ax_left, ax_right)):
        ax.plot(
            gm_years,
            gm_values,
            label=extended_label if i < 1 else None,
            color=palette[extended_label],
            linewidth=2,
            # s=30,
        )

        for label, pdf in input_sources.items():
            ax.plot(
                pdf["year"],
                pdf["value"],
                label=label if i < 1 else None,
                color=palette[label],
                linewidth=2,
                # s=30,
            )

    # Over the top of the line, so the stretches which are not a source
    # are marked out on the line itself rather than beside it.
    # Each piece runs up to where the next one starts, so the line
    # changes colour without breaking.
    if pre_industrial is not None:
        pre_industrial_year, pre_industrial_value = pre_industrial
        assumed_constant = gm_years <= pre_industrial_year
        for i, ax in enumerate((ax_left, ax_right)):
            ax.plot(
                gm_years[assumed_constant],
                gm_values[assumed_constant],
                color=ASSUMED_CONSTANT_COLOUR,
                linewidth=2,
                # The pair shares one legend, which is built on the left half
                label="Assumed constant" if i < 1 else None,
            )

    if fit_period is not None:
        first_fit_year, last_fit_year = fit_period
        # Out to the year either side of the fitted ones, which are the
        # points the fit was anchored to: the stretch then runs from the
        # pre-industrial marker to the first year a source covers,
        # with no sliver of the base colour left showing at either end.
        fitted = (gm_years >= first_fit_year - 1) & (gm_years <= last_fit_year + 1)
        for i, ax in enumerate((ax_left, ax_right)):
            ax.plot(
                gm_years[fitted],
                gm_values[fitted],
                color=FIT_PERIOD_COLOUR,
                linewidth=2,
                label="Fit period" if i < 1 else None,
            )

    if pre_industrial is not None:
        for i, ax in enumerate((ax_left, ax_right)):
            # Both halves: which half the pre-industrial year lands in
            # depends on the gas, and each half shows only what is inside
            # its own limits.
            ax.plot(
                pre_industrial_year,
                pre_industrial_value,
                marker="*",
                markersize=9,
                linestyle="none",
                color=PRE_INDUSTRIAL_COLOUR,
                label="Pre-industrial value" if i < 1 else None,
                zorder=5,
            )

    add_break_lines_and_setup(
        ax_left,
        ax_right,
        gm_da["year"].min(),
        split_year,
        gm_da["year"].max(),
        gm_da.attrs["units"],
        legend_loc=legend_loc,
    )
    if pre_industrial is not None:
        # The pre-industrial value is usually the lowest thing on the panel,
        # and often exactly zero, so without a little room under it
        # its marker is drawn half outside the axes.
        for ax in (ax_left, ax_right):
            y_min, y_max = ax.get_ylim()
            margin = PRE_INDUSTRIAL_MARKER_MARGIN * (y_max - y_min)
            ax.set_ylim(min(y_min, pre_industrial[1] - margin), y_max)


def plot_extension_pieces(  # noqa: PLR0913
    da: xr.DataArray,
    ax_left: matplotlib.axes.Axes,
    ax_right: matplotlib.axes.Axes,
    pieces: Mapping[str, np.typing.ArrayLike],
    split_year: int = 1900,
    legend_loc: str = "upper left",
) -> None:
    """
    Plot an extended timeseries, coloured by where each year's value came from

    A component which has been extended in time looks like one line,
    but only part of it is data: the rest is whatever assumption
    filled the years the data does not reach.
    Colouring the line by its source is the only thing which says so.

    Parameters
    ----------
    da
        Extended timeseries to plot, with a `year` dimension

    ax_left
        Left half of the broken axis

    ax_right
        Right half of the broken axis

    pieces
        Years which came from each source, labelled by source

        Every year in `da` must appear in exactly one of these,
        and they are drawn, and listed in the legend, in the order given.

    split_year
        Year the axis is broken at

    legend_loc
        Where to put the panel's legend

    Raises
    ------
    AssertionError
        Some year in `da` was not attributed to a source
    """
    pdf = da.to_pandas().rename("value").to_frame().reset_index()
    pdf["source"] = None
    for source, years in pieces.items():
        pdf.loc[pdf["year"].isin(np.asarray(years)), "source"] = source

    if pdf["source"].isnull().any(axis=None):
        unattributed = pdf.loc[pdf["source"].isnull(), "year"]
        msg = (
            "Every year must come from somewhere, "
            f"but {len(unattributed)} years were not attributed to a source, "
            f"the first of them {unattributed.iloc[0]}"
        )
        raise AssertionError(msg)

    # So the legend reads in the order the pieces were given,
    # which is the order the reader meets them in as they go along the line.
    pdf["source"] = pd.Categorical(pdf["source"], categories=list(pieces), ordered=True)

    for ax in (ax_left, ax_right):
        sns.scatterplot(
            pdf,
            x="year",
            y="value",
            hue="source",
            palette=source_sequence_palette(pieces),
            ax=ax,
            s=25,
            edgecolor=None,
            alpha=0.7,
        )

    add_break_lines_and_setup(
        ax_left,
        ax_right,
        int(da["year"].min()),
        split_year,
        int(da["year"].max()),
        da.attrs["units"],
        legend_loc=legend_loc,
    )


def plot_pc_timeseries_regression(  # noqa: PLR0913
    lat_grad_info: xr.Dataset,
    timeseries_data: xr.Dataset,
    timeseries_name: str,
    regression_info: dict[str, tuple[float, str]],
    ax: matplotlib.axes.Axes,
    x_unit: str,
    pcs_key: str = "principal-components",
    eof: int = 0,
    ur=openscm_units.unit_registry,
) -> None:
    """
    Plot the regression between a PC and a timeseries
    """
    pc_da = lat_grad_info[pcs_key].sel(eof=eof)

    common_years = np.intersect1d(pc_da["year"], timeseries_data["year"])

    timeseries_da = get_only_data_variable(timeseries_data.sel(year=common_years))
    timeseries_da.to_pandas()
    timeseries_da_units = timeseries_da.attrs["units"]
    conversion_factor = ur(timeseries_da_units).to(x_unit).m
    timeseries_values = timeseries_da.values * conversion_factor

    pc_units = pc_da.attrs["units"]
    pc_values = pc_da.sel(year=common_years).values

    ax.scatter(
        x=timeseries_values,
        y=pc_values,
        label="raw data",
        marker="x",
        s=30,
        color="tab:blue",
        alpha=0.7,
    )
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()

    regression_gradient = (
        ur.Quantity(regression_info["m"][0], regression_info["m"][1])
        .to(f"{pc_units} / ({x_unit})")
        .m
    )
    regression_y_int = (
        ur.Quantity(regression_info["c"][0], regression_info["c"][1]).to(pc_units).m
    )
    ax.axline(
        xy1=(0, regression_y_int),
        slope=regression_gradient,
        label="regression",
        linestyle="-",
        color="tab:orange",
        alpha=0.9,
    )
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    ax.set_xlabel(f"{timeseries_name} [{x_unit}]", fontsize="small")
    ax.set_ylabel(f"PC{eof} [{pc_units}]", fontsize="small")
    ax.tick_params(labelsize="small")

    add_compact_legend(ax, loc="best")


def plot_pcs_extended(  # noqa: PLR0913
    eof_pieces: xr.Dataset,
    ax_left: matplotlib.axes.Axes,
    ax_right: matplotlib.axes.Axes,
    pcs_key: str = "principal-components",
    split_year: int = 1950,
    # Should be numpy array of int, anyway
    pieces: dict[int, dict[str, list[int]]] | None = None,
    legend_loc: str = "upper left",
    max_ticks_per_half: int | None = None,
) -> None:
    """
    Plot extended latitudinal gradient PCs
    """
    pcs_da = eof_pieces[pcs_key]
    pcs_df = pcs_da.to_pandas().stack().rename("value").to_frame().reset_index()

    for eof, info in pieces.items():
        for source, years in info.items():
            pcs_df.loc[
                (pcs_df["eof"] == eof) & (pcs_df["year"].isin(years)), "source"
            ] = source

    if pcs_df["source"].isnull().any(axis=None):
        raise AssertionError

    for i, ax in enumerate((ax_left, ax_right)):
        sns.scatterplot(
            pcs_df,
            x="year",
            y="value",
            # hue="eof",
            # style="source",
            style="eof",
            hue="source",
            palette=extension_source_palette(pcs_df["source"].unique()),
            ax=ax,
            s=25,
            edgecolor=None,
            alpha=0.7,
        )

    add_break_lines_and_setup(
        ax_left,
        ax_right,
        pcs_da["year"].min(),
        split_year,
        pcs_da["year"].max(),
        pcs_da.attrs["units"],
        legend_loc=legend_loc,
        max_ticks_per_half=max_ticks_per_half,
    )

    legend = ax_left.get_legend()
    if legend is None:
        raise AssertionError

    handles = legend.legend_handles
    labels = [text.get_text() for text in legend.get_texts()]
    # Add fake handles and legends for other EOFs so 2 col splits as we want.
    for _ in range(len(pcs_df["source"].unique()) - len(pcs_df["eof"].unique())):
        handles.append(matplotlib.lines.Line2D([], [], color="none", label=""))
        labels.append("")

    title = legend.get_title().get_text()

    add_compact_legend(
        ax_left, handles=handles, labels=labels, title=title, ncols=2, loc=legend_loc
    )
    ax_left.get_legend().get_title().set_fontsize("x-small")


def plot_flying_carpet(
    native_resolution: xr.Dataset, ax: matplotlib.axes.Axes, tick_nearest: int = 5
) -> None:
    """
    Plot flying carpet
    """
    tmp = convert_year_month_to_time(get_only_data_variable(native_resolution))
    tmp = tmp.assign_coords(
        time=tmp["time"].dt.year + tmp["time"].dt.month / 12 - 1 / 24
    )
    mesh = tmp.plot.surface(
        x="time",
        y="lat",
        ax=ax,
        cmap=INTERPOLATION_COLOUR_MAP,
        levels=30,
        add_colorbar=False,
        # alpha=0.7,
    )
    ax.view_init(15, -135, 0)  # type: ignore
    # No label on the vertical axis. Matplotlib works out which way up to
    # write a 3D axis' label from the shape of the axes' box, and on this
    # panel it gets the vertical one upside down as often as not.
    # The colour bar carries what the label would have said.
    ax.set_zlabel("")  # type: ignore
    ax.set_ylabel(r"latitude [$^{\circ}$N]", fontsize="small")
    ax.set_yticks([-45, 45])
    ax.set_xlabel("time", fontsize="small")
    max_year = int(native_resolution["year"].max())
    other_tick = tick_nearest * np.floor(max_year / tick_nearest) - 5
    ax.set_xticks([other_tick, max_year])
    ax.tick_params(labelsize="small", pad=0.0)

    return mesh


def plot_monthly_means(
    pda: xr.Dataset,
    ax: matplotlib.axes.Axes,
) -> None:
    """
    Plot monthly means
    """
    tmp = convert_year_month_to_time(get_only_data_variable(pda))
    tmp = tmp.assign_coords(
        time=tmp["time"].dt.year + tmp["time"].dt.month / 12 - 1 / 24
    )
    pdf_l = []
    for region, rda in tmp.groupby("Region"):
        tmp_df = (
            rda.sel(Region=region).to_pandas().rename("value").to_frame().reset_index()
        )
        tmp_df["Region"] = region
        pdf_l.append(tmp_df)

    pdf = pd.concat(pdf_l)

    sns.scatterplot(
        pdf,
        x="time",
        y="value",
        hue="Region",
        palette=region_palette(pdf["Region"].unique()),
        ax=ax,
        s=15,
        alpha=0.7,
        edgecolor=None,
    )

    ax.set_ylabel(f"[{tmp.attrs['units']}]", fontsize="small")
    ax.set_xlabel("time", fontsize="small")
    ax.tick_params(labelsize="small")
    ax.set_xlim(pda["year"].min(), pda["year"].max())
    ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
    compact_existing_legend(ax, loc="best")


def plot_yearly_means(
    pda: xr.Dataset,
    ax_left: matplotlib.axes.Axes,
    ax_right: matplotlib.axes.Axes,
    split_year: int = 1850,
) -> None:
    """
    Plot global-mean derived from the observational network
    """
    tmp = get_only_data_variable(pda)
    pdf_l = []
    for region, rda in tmp.groupby("Region"):
        tmp_df = (
            rda.sel(Region=region).to_pandas().rename("value").to_frame().reset_index()
        )
        tmp_df["Region"] = region
        pdf_l.append(tmp_df)

    pdf = pd.concat(pdf_l)

    for i, ax in enumerate((ax_left, ax_right)):
        sns.scatterplot(
            pdf,
            x="year",
            y="value",
            hue="Region",
            palette=region_palette(pdf["Region"].unique()),
            ax=ax,
            s=15,
            alpha=0.7,
            edgecolor=None,
            # legend=False,
        )

    add_break_lines_and_setup(
        ax_left,
        ax_right,
        pda["year"].min(),
        split_year,
        pda["year"].max(),
        tmp.attrs["units"],
        # The record is flat and low until its last century, so the top of
        # the left half is the only part of this panel with room in it,
        # and centring keeps the legend off the axis' tick labels.
        legend_loc="upper center",
    )


def thin_ticks(ticks: np.typing.ArrayLike, max_ticks: int) -> np.typing.NDArray:
    """
    Take every nth tick, so that no more than `max_ticks` are left

    A colour bar which stands for a count wants a tick per count,
    right up until there are more counts than the colour bar has room for,
    at which point the labels run into each other
    and the colour bar says less than it would with a handful of ticks.

    Parameters
    ----------
    ticks
        Ticks we would show if there were room for all of them

    max_ticks
        The most ticks to leave

    Returns
    -------
        Every nth tick of `ticks`, for the smallest n which leaves
        no more than `max_ticks` of them
    """
    ticks = np.asarray(ticks)
    if ticks.size <= max_ticks:
        return ticks

    return ticks[:: int(np.ceil(ticks.size / max_ticks))]


def add_colour_bar(  # noqa: PLR0913
    fig: matplotlib.figure.Figure,
    mappable: matplotlib.cm.ScalarMappable,
    cax: matplotlib.axes.Axes,
    label: str,
    ticks: np.typing.ArrayLike | None = None,
    max_ticks: int = 6,
    label_on_top: bool = False,
    **kwargs: object,
) -> matplotlib.colorbar.Colorbar:
    """
    Add a colour bar to the figure

    Parameters
    ----------
    fig
        Figure to add the colour bar to

    mappable
        Mappable to draw the colour bar for

    cax
        Axes in which to draw the colour bar

        Where this sits is up to the figure's layout,
        see [local.historical_ghg_forcing_for_cmip7.layout][].

    label
        Label for the colour bar

    ticks
        Ticks to show on the colour bar

        If not supplied, the ticks are left to matplotlib.

    max_ticks
        The most ticks to actually show out of `ticks`

        See [thin_ticks][].

    label_on_top
        Put the label above the colour bar rather than beside it

        Useful where there is no width to spare beside the colour bar.

    **kwargs
        Passed on to `fig.colorbar`

    Returns
    -------
        The colour bar which was added
    """
    colour_bar = fig.colorbar(mappable, cax=cax, **kwargs)
    if ticks is not None:
        colour_bar.set_ticks(thin_ticks(ticks, max_ticks=max_ticks))

    if label_on_top:
        colour_bar.ax.set_title(label, fontsize="small", loc="left")
    else:
        colour_bar.set_label(label, fontsize="small")

    colour_bar.ax.tick_params(labelsize="small")

    return colour_bar
