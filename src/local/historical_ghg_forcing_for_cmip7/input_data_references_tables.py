"""
Generation of the tables which list the references of each gas' input data

The original run recorded the sources each gas depends on
in its `dependencies.db` (this is also what was written into the output files).
Rather than write these out by hand for every gas,
we read them from there, add the few sources the original run did not record
and write the tables here.
The build inlines them in the appendices in place of their tag.

The gas names are written the way the manuscript's source text writes them
(e.g. `C_2F_6`, `CFC12`), so the manuscript's replacements
format them exactly as they are formatted in the text around them.
"""

from __future__ import annotations

import sqlite3
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

from loguru import logger

from local.historical_ghg_forcing_for_cmip7.appendix_figures import (
    MANUSCRIPT_GAS_NAMES as CFC12_AND_C4F10_LIKE_MANUSCRIPT_GAS_NAMES,
)
from local.historical_ghg_forcing_for_cmip7.c4f10_like_methods_figure import (
    C4F10_LIKE_GASES,
)
from local.historical_ghg_forcing_for_cmip7.cfc12_like_methods_figure import (
    CFC12_LIKE_GASES,
)

DEPENDENCIES_DB = Path("data") / "processed" / "dependencies.db"
"""The original run's record of the sources each gas depends on

Relative to the directory which holds the original run's bundle.
Its `source` table holds each source's reference
and its `dependencies` table the sources of each gas.
"""

MANUSCRIPT_GAS_NAMES = {
    "co2": "CO_2",
    "ch4": "CH_4",
    "n2o": "N_2O",
    **CFC12_AND_C4F10_LIKE_MANUSCRIPT_GAS_NAMES,
    "c8f18": "C_8F_18",
}
"""How the manuscript's source text writes each gas' name"""

NOAA = "NOAA"
AGAGE = "AGAGE"
OTHER = "Other"
COLUMNS = (NOAA, AGAGE, OTHER)
"""The tables' columns of references, in order"""

COLUMN_HEADERS = {NOAA: "NOAA", AGAGE: "AGAGE", OTHER: "Other sources"}
"""Header of each column"""

COLUMN_WIDTHS = {NOAA: "3.2cm", AGAGE: "6.5cm", OTHER: "4.3cm"}
"""Width of each column

Together with the gas' column, these fill the width of the page.
"""

SOURCES_TO_SKIP = ("Nicholls et al., 2025 (in-prep)",)
"""Sources in the original run's record which don't go in the tables

The original run's placeholder for this manuscript.
"""

SOURCES_RECORDED_BUT_NOT_USED = {
    gas: ("Nicholls et al., 2020",) for gas in C4F10_LIKE_GASES
}
"""Sources the original run recorded for a gas, but which its output doesn't depend on

For the gases processed like C4F10, the original run's
`1404_c4f10-like_derive-latitudinal-gradient` notebook loads the RCMIP emissions
(so records them as a dependency) and regresses the latitudinal gradient against them,
but the regression isn't used: the code which would apply it is commented out
and the latitudinal gradient is instead held at zero before the Droste et al. data
and linearly extrapolated after it.
"""

NOAA_HATS_COMBINED_BIBKEYS = {
    "ccl4": "dutton_hats-ccl4_2022",
    "cfc11": "dutton_hats-cfc11_2022",
    "cfc113": "dutton_hats-cfc113_2022",
    "cfc12": "dutton_hats-cfc12_2022",
    "n2o": "dutton_hats-n2o_2022",
    "sf6": "dutton_hats-sf6_2022",
}
"""Bibtex key of NOAA's HATS combined product of each gas which has one"""

NOAA_HATS_FLASK_BIBKEY = "montzka_ods_1999"
"""Bibtex key of the reference the original run recorded for NOAA's HATS flask data"""

SOURCE_BIBKEYS: dict[str, tuple[str, ...]] = {
    "NOAA co2 surface-flask": ("lan_co2-flask_2024",),
    "NOAA co2 in-situ": ("thoning_co2-in-situ_2024",),
    "NOAA ch4 surface-flask": ("lan_ch4-flask_2024",),
    "NOAA ch4 in-situ": ("thoning_ch4-in-situ_2024",),
    # The reference AGAGE asks for, whichever of its networks the data is from
    "AGAGE": ("prinn_ale-gage-agage_2000",),
    "AGAGE ALE": ("prinn_ale-gage-agage_2000",),
    "AGAGE GAGE": ("prinn_ale-gage-agage_2000",),
    "AGAGE Arnold et al. 2013": ("arnold_nf3_2013",),
    "AGAGE Arnold et al. 2014": ("arnold_hfc4310mee_2014",),
    "AGAGE Arnold et al. 2018": ("arnold_cf4-nf3-east-asia_2018",),
    "AGAGE Babbin et al. 2020": ("babbin_n2o-south-pacific_2020",),
    "AGAGE Claxton et al. 2020": ("claxton_chlorocarbons_2020",),
    # Published online in 2018, in print in 2019
    "AGAGE Fang et al. 2018": ("fang_chloroform_2019",),
    "AGAGE Ganesan et al. 2020": ("ganesan_marine-n2o_2020",),
    "AGAGE Hossaini et al. 2019": ("hossaini_vsls-chlorine_2019",),
    "AGAGE Lunt et al. 2015": ("lunt_hfcs_2015",),
    "AGAGE Lunt et al. 2018": ("lunt_ccl4_2018",),
    "AGAGE Miller et al. 2010": ("miller_hfc23_2010",),
    "AGAGE Montzka et al. 2021": ("montzka_cfc11_2021",),
    "AGAGE Mühle et al. 2009": ("muhle_so2f2_2009",),
    "AGAGE Mühle et al. 2010": ("muhle_pfcs_2010",),
    "AGAGE Mühle et al. 2019": ("muhle_cc4f8_2019",),
    "AGAGE Nevison et al. 2004": ("nevison_n2o-cfcs-seasonal_2004",),
    "AGAGE Nevison et al. 2007": ("nevison_n2o-variability_2007",),
    "AGAGE O'Doherty et al. 2009": ("odoherty_hfc125_2009",),
    "AGAGE O'Doherty et al. 2014": ("odoherty_hfc143a-hfc32_2014",),
    "AGAGE Park et al. 2018": ("park_ccl4_2018",),
    "AGAGE Park et al. 2021": ("park_cfc11_2021",),
    "AGAGE Patra et al. 2021": ("patra_ch3ccl3_2021",),
    "AGAGE Prinn et al. 2018": ("prinn_agage_2018",),
    "AGAGE Rigby et al. 2008": ("rigby_methane-growth_2008",),
    "AGAGE Rigby et al. 2010": ("rigby_sf6_2010",),
    "AGAGE Rigby et al. 2013": ("rigby_lifetimes_2013",),
    "AGAGE Rigby et al. 2014": ("rigby_synthetic-ghgs_2014",),
    "AGAGE Rigby et al. 2017": ("rigby_methane-oxidation_2017",),
    "AGAGE Rigby et al. 2019": ("rigby_cfc11_2019",),
    "AGAGE Saikawa et al. 2012": ("saikawa_hcfc22_2012",),
    "AGAGE Saikawa et al. 2014": ("saikawa_n2o_2014",),
    "AGAGE Say et al. 2021": ("say_pfcs_2021",),
    "AGAGE Simmonds et al. 2016": ("simmonds_hfc152a_2016",),
    "AGAGE Simmonds et al. 2017": ("simmonds_hcfcs-hfcs_2017",),
    "AGAGE Simmonds et al. 2018": ("simmonds_hfc23-hcfc22_2018",),
    "AGAGE Simmonds et al. 2020": ("simmonds_sf6_2020",),
    "AGAGE Stanley et al. 2020": ("stanley_hfc23_2020",),
    "AGAGE Thompson et al. 2013": ("thompson_n2o-variability_2013",),
    "AGAGE Thompson et al. 2014": ("thompson_n2o-inversion_2014",),
    "AGAGE Thompson et al. 2014 b": ("thompson_transcom-n2o-part-2_2014",),
    "AGAGE Thompson et al. 2014 c": ("thompson_transcom-n2o-part-1_2014",),
    "AGAGE Trudinger et al. 2016": ("trudinger_pfcs_2016",),
    "AGAGE Vollmer et al. 2011": ("vollmer_hfcs_2011",),
    "AGAGE Vollmer et al. 2016": ("vollmer_halons_2016",),
    "AGAGE Vollmer et al. 2018": ("vollmer_cfcs_2018",),
    "AGAGE Vollmer et al. 2021": ("vollmer_hcfcs_2021",),
    "Adam et al., 2024": ("adam_hfc23_2024",),
    "Ahn et al., 2012": ("ahn_co2_2012",),
    "Azharuddin et al., 2024": ("azharuddin_n2o-holocene_2024",),
    "Bauska et al., 2015": ("bauska_co2_2015",),
    "Bernard et al., 2006": ("bernard_n2o-firn_2006",),
    "Droste et al., 2020": ("droste_pfcs_2020",),
    "EPICA": ("epica_edml-ch4_2006",),
    # The original run recorded the paper's reference with the DOI of its dataset
    # (https://doi.org/10.15784/601693). We cite the paper.
    "Ghosh et al., 2023": ("ghosh_n2o_2023",),
    # The paper (which is what the original run recorded) and the data
    "HadCRUT5": ("morice_hadcrut5_2021", "met-office_hadcrut5-data_2025"),
    "Ishijima et al., 2007": ("ishijima_n2o-firn_2007",),
    "King et al., 2024": ("king_co2_2024",),
    "Law Dome ice core": ("rubino_law-dome-dataset_2019",),
    "Meinshausen et al., 2017": ("meinshausen_historical-ghgs_2017",),
    "Meinshausen et al., 2020": ("meinshausen_ssp-ghgs_2020",),
    "Menking et al., 2025 (in-prep.)": ("menking_law-dome_2025",),
    # The dataset (which is what the original run recorded)
    # and the paper it is a supplement to
    "NEEM": ("rhodes-brook_neem-ch4-dataset_2019", "rhodes_neem-ch4_2013"),
    "Nicholls et al., 2020": ("nicholls_rcmip_2020",),
    "Park et al., 2012": ("park_n2o-isotopes_2012",),
    "Prokopiou et al., 2018": ("prokopiou_n2o-isotopes_2018",),
    # The original run called this 2006, but the reference it recorded is from 2003
    "Roeckmann et al., 2006": ("rockmann_n2o-isotopes_2003",),
    "Ryu et al., 2020": ("ryu_n2o_2020",),
    "Schilt et al., 2010": ("schilt_n2o_2010",),
    # The references the record's file header asks for
    # (the original run only recorded the first)
    "Scripps - Law Dome merged CO2 record": (
        "keeling_co2-exchanges_2001",
        "rubino_law-dome-dataset_2019",
    ),
    "Trudinger et al., 2016": ("trudinger_pfcs_2016",),
    "Velders et al., 2022": ("velders_hfcs_2022",),
    "WMO 2022 Ozone Assessment Ch. 7": ("wmo_ozone-ch7_2022",),
    "Western et al., 2024": ("western_hcfcs_2024",),
}
"""Bibtex keys of each source, by its name in the original run's record

NOAA's HATS data is not here, because its sources are named per gas,
see [get_source_bibkeys][].
"""

SOURCES_NOT_IN_DEPENDENCIES_DB: dict[str, dict[str, tuple[str, ...]]] = {
    # PRIMAP-hist emissions, which the latitudinal gradient is regressed against
    "co2": {OTHER: ("gutschow_primap-hist-dataset_2024", "gutschow_primap-hist_2016")},
    "ch4": {OTHER: ("gutschow_primap-hist-dataset_2024", "gutschow_primap-hist_2016")},
}
"""Bibtex keys of the sources which the original run used but did not record

By gas, then by the column they go in.
"""


@dataclass(frozen=True)
class GasGroup:
    """
    A group of gases whose references are shown together
    """

    name: str
    """Name of the group, used in the tables' labels"""

    description: str
    """How the group is described in the tables' captions"""

    gases: tuple[str, ...]
    """The gases in the group, in the order of their rows"""

    n_tables: int = 1
    """
    Number of tables to split the group over

    A table can't be split over pages, so each has to fit on one.
    """


GAS_GROUPS = (
    GasGroup(
        name="co2-ch4-n2o",
        description="CO_2, CH_4 and N_2O",
        gases=("co2", "ch4", "n2o"),
    ),
    GasGroup(
        name="cfc12-like",
        description="the gases processed like CFC12",
        gases=CFC12_LIKE_GASES,
        n_tables=2,
    ),
    GasGroup(
        name="c4f10-like-and-c8f18",
        description="the gases processed like C_4F_10 and C_8F_18",
        gases=(*C4F10_LIKE_GASES, "c8f18"),
    ),
)
"""The groups of gases, in the order of their tables"""


FIRST_AUTHOR_OVERRIDES = {
    # Daniel et al. (2022)
    "wmo_ozone-ch7_2022": "daniel",
}
"""First author of the references whose bibtex key doesn't start with it

Only used for sorting, see [sort_bibkeys][].
"""


def sort_bibkeys(bibkeys: Iterable[str]) -> tuple[str, ...]:
    """
    Sort bibtex keys by first author, then year

    This relies on our keys being written `<first author>_<topic>_<year>`
    (for those which aren't, see [FIRST_AUTHOR_OVERRIDES][]).

    Parameters
    ----------
    bibkeys
        Bibtex keys to sort

    Returns
    -------
    :
        The sorted keys
    """
    return tuple(
        sorted(
            bibkeys,
            key=lambda bibkey: (
                FIRST_AUTHOR_OVERRIDES.get(bibkey, bibkey.split("_")[0]),
                bibkey.split("_")[-1],
                bibkey,
            ),
        )
    )


def get_source_bibkeys(source: str) -> tuple[str, ...]:
    """
    Get the bibtex keys of a source

    Parameters
    ----------
    source
        The source's name in the original run's record

    Returns
    -------
    :
        The source's bibtex keys

    Raises
    ------
    KeyError
        We have no bibtex keys for `source`
    """
    if source.startswith("NOAA ") and source.endswith(" hats"):
        gas = source.removeprefix("NOAA ").removesuffix(" hats")

        return (NOAA_HATS_COMBINED_BIBKEYS.get(gas, NOAA_HATS_FLASK_BIBKEY),)

    return SOURCE_BIBKEYS[source]


def get_source_column(source: str) -> str:
    """
    Get the column a source's references go in

    Parameters
    ----------
    source
        The source's name in the original run's record

    Returns
    -------
    :
        The column, one of [COLUMNS][]
    """
    if source.startswith("NOAA "):
        return NOAA

    if source.startswith("AGAGE"):
        return AGAGE

    return OTHER


def get_bibkeys_by_column(gas: str, bundle_dir: Path) -> dict[str, tuple[str, ...]]:
    """
    Get the bibtex keys of a gas' input data, by the column they go in

    Parameters
    ----------
    gas
        Gas of interest

    bundle_dir
        Directory which holds the original run's bundle

    Returns
    -------
    :
        The bibtex keys in each of [COLUMNS][], sorted

        A reference only appears once.
        Where the original run recorded it both under AGAGE and on its own
        (because AGAGE asks for it to be cited and we also use its data directly),
        it goes in the other sources' column.

    Raises
    ------
    AssertionError
        The original run recorded no sources for `gas`
    """
    con = sqlite3.connect(bundle_dir / DEPENDENCIES_DB)
    try:
        sources = [
            row[0]
            for row in con.execute(
                "SELECT short_name FROM dependencies WHERE gas = ?", (gas,)
            )
        ]
    finally:
        con.close()

    if not sources:
        msg = f"No sources for {gas=} in {bundle_dir / DEPENDENCIES_DB}"
        raise AssertionError(msg)

    res: dict[str, set[str]] = {column: set() for column in COLUMNS}
    for source in sources:
        if source in SOURCES_TO_SKIP:
            continue

        if source in SOURCES_RECORDED_BUT_NOT_USED.get(gas, ()):
            continue

        res[get_source_column(source)].update(get_source_bibkeys(source))

    for column, bibkeys in SOURCES_NOT_IN_DEPENDENCIES_DB.get(gas, {}).items():
        res[column].update(bibkeys)

    res[AGAGE] -= res[OTHER]

    return {column: sort_bibkeys(res[column]) for column in COLUMNS}


def get_cite(bibkeys: tuple[str, ...]) -> str:
    """
    Get the latex which cites references

    Parameters
    ----------
    bibkeys
        Bibtex keys of the references

    Returns
    -------
    :
        The latex, a dash if there is nothing to cite
    """
    if not bibkeys:
        return "--"

    return rf"\citet{{{','.join(bibkeys)}}}"


def get_table(
    gases: tuple[str, ...],
    bibkeys_by_column: dict[str, dict[str, tuple[str, ...]]],
    caption: str,
    label: str,
    placement: str,
) -> str:
    """
    Get the latex of a table

    Parameters
    ----------
    gases
        The gases to show, in the order of their rows

    bibkeys_by_column
        The bibtex keys to cite for each gas, by the column they go in

        A column which is empty for every gas is left out of the table.

    caption
        The table's caption

    label
        The table's label

    placement
        Where latex is allowed to put the table (e.g. `"t"`)

    Returns
    -------
    :
        The latex of the table
    """
    columns = [
        column
        for column in COLUMNS
        if any(bibkeys_by_column[gas][column] for gas in gases)
    ]

    # The cells are ragged right paragraphs, so the citations wrap
    # (without stretched spaces), hence `\tabularnewline` rather than `\\`.
    rows = [
        " & ".join(
            (
                MANUSCRIPT_GAS_NAMES[gas],
                *(
                    rf"\raggedright {get_cite(bibkeys_by_column[gas][column])}"
                    for column in columns
                ),
            )
        )
        + r" \tabularnewline"
        for gas in gases
    ]

    column_spec = "".join(f"p{{{COLUMN_WIDTHS[column]}}}" for column in columns)

    return "\n".join(
        [
            rf"\begin{{table*}}[{placement}]",
            r"\caption{",
            f"    {caption}",
            "}",
            rf"\label{{{label}}}",
            r"\footnotesize",
            rf"\begin{{tabular}}{{l{column_spec}}}",
            r"\tophline",
            " & ".join(("Gas", *(COLUMN_HEADERS[c] for c in columns))) + r" \\",
            r"\middlehline",
            *(f"    {row}" for row in rows),
            r"\bottomhline",
            r"\end{tabular}",
            r"\end{table*}",
        ]
    )


def get_group_tables(
    group: GasGroup, bundle_dir: Path, first_placement: str = "t"
) -> list[str]:
    """
    Get the latex of the tables of a group of gases

    Parameters
    ----------
    group
        The group of gases

    bundle_dir
        Directory which holds the original run's bundle

    first_placement
        Where latex is allowed to put the group's first table

        The others can go at the top of a page.

    Returns
    -------
    :
        The latex of each of the group's tables
    """
    bibkeys_by_column = {
        gas: get_bibkeys_by_column(gas, bundle_dir) for gas in group.gases
    }

    # Split the gases as evenly as possible, keeping their order
    n_per_table = -(-len(group.gases) // group.n_tables)
    res = []
    for i in range(group.n_tables):
        gases = group.gases[i * n_per_table : (i + 1) * n_per_table]

        caption = f"References for the generation of {group.description}"
        label = f"tab:input-data-references-{group.name}"
        if group.n_tables > 1:
            caption = (
                f"{caption} ({MANUSCRIPT_GAS_NAMES[gases[0]]} "
                f"to {MANUSCRIPT_GAS_NAMES[gases[-1]]})"
            )
            label = f"{label}-{i + 1}"

        res.append(
            get_table(
                gases,
                bibkeys_by_column,
                caption=f"{caption}.",
                label=label,
                placement=first_placement if i == 0 else "t",
            )
        )

    return res


def generate_input_data_references_tables(
    outfile: Path,
    bundle_dir: Path,
    force_rerun: bool = False,
) -> Path:
    """
    Generate the tables which list the references of each gas' input data

    Parameters
    ----------
    outfile
        File in which to write the tables

    bundle_dir
        Directory which holds the original run's bundle

    force_rerun
        Re-generate the tables, even if the output file already exists

    Returns
    -------
    :
        `outfile`
    """
    if outfile.exists() and not force_rerun:
        logger.info(f"Using existing {outfile}")
        return outfile

    # The first table follows its appendix's heading,
    # so (like the first figure of the other appendices)
    # get latex to keep it after the heading, at the bottom of that page
    # or at the top of the next, rather than floating it above the heading.
    tables = [
        table
        for i, group in enumerate(GAS_GROUPS)
        for table in get_group_tables(
            group, bundle_dir, first_placement="hbt!" if i == 0 else "t"
        )
    ]

    outfile.parent.mkdir(exist_ok=True, parents=True)
    logger.info(f"Writing {outfile}")
    outfile.write_text("\n\n".join(tables) + "\n")

    return outfile
