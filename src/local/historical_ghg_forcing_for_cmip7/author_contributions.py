"""
Generation of the author contributions statement

Each author in the manuscript's metadata file lists their contributions as tags.
We turn these into one sentence per contribution,
"<initials of the authors who made it> <what the contribution was>.",
which the build inlines in the author contributions in place of its tag.
"""

from __future__ import annotations

import tomllib
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from loguru import logger

CONTRIBUTIONS = {
    "dataset-production": (
        "produced the dataset, including writing the underlying software"
    ),
    "dataset-contributed": (
        "contributed to improvements of production of the dataset "
        "(such as considering new methods and examining outputs)"
    ),
    "conceptualisation": "conceptualised the dataset",
    # "conceptualisation-support": "supported the conceptualisation",
    # "software": "developed the software used to produce the dataset",
    "input-data": "provided input data and advised on its use",
    "input-data-noaa-ccgg-flask": (
        "provided the CO_2 and CH_4 surface flask data from the "
        "NOAA Global Monitoring Laboratory (GML) "
        "Carbon Cycle Cooperative Global Air Sampling Network"
    ),
    "input-data-noaa-ccgg-in-situ-co2": (
        "provided the in situ CO_2 data from the NOAA GML observatories"
        # "at Barrow (Alaska), Mauna Loa (Hawaii), American Samoa and the South Pole"
    ),
    "input-data-noaa-ccgg-in-situ-ch4": (
        "provided the in situ CH_4 data from the NOAA GML observatories"
        # "at Barrow (Alaska) and Mauna Loa (Hawaii)"
    ),
    "input-data-noaa-hats-combined": (
        "provided the NOAA GML Halocarbons and other Atmospheric Trace Species (HATS) "
        "combined data products (CCl_4, CFC-11, CFC-113, CFC-12, N_2O and SF_6)"
    ),
    "input-data-noaa-hats-flask": "provided the NOAA GML HATS flask data",
    "input-data-agage-ale-gage": (
        "provided the data from AGAGE's predecessor networks, ALE and GAGE"
    ),
    "input-data-agage-archive": "curated the AGAGE data archive",
    "input-data-agage-cape-grim": (
        "provided the AGAGE data from Kennaook/Cape Grim (Tasmania, Australia)"
    ),
    "input-data-agage-mace-head": ("provided the AGAGE data from Mace Head (Ireland)"),
    "input-data-agage-ragged-point": (
        "provided the AGAGE data from Ragged Point (Barbados)"
    ),
    "input-data-agage-trinidad-head-samoa": (
        "provided the AGAGE data from Trinidad Head (California, USA) "
        "and Cape Matatula (American Samoa)"
    ),
    "input-data-agage-jungfraujoch": (
        "provided the AGAGE data from Jungfraujoch (Switzerland)"
    ),
    "input-data-agage-zeppelin": (
        r"provided the AGAGE data from Zeppelin (Ny-\AA{}lesund, Svalbard, Norway)"
    ),
    "input-data-agage-monte-cimone": (
        "provided the AGAGE data from Monte Cimone (Italy)"
    ),
    "input-data-agage-gosan": (
        "provided the AGAGE data from Gosan (Jeju Island, Republic of Korea)"
    ),
    "input-data-agage-tacolneston": ("provided the AGAGE data from Tacolneston (UK)"),
    "input-data-agage-taunus": (
        "provided the AGAGE data from Taunus Observatory (Germany)"
    ),
    "input-data-agage-hfc23": (
        "provided data for global-mean HFC-23 concentrations"
        # r"from the AGAGE 12-box model of \citet{adam_hfc23_2024}"
    ),
    "input-data-law-dome-menking": (
        "provided the spline fits to the CSIRO Law Dome ice core and firn "
        "CO_2 and N_2O data"
        # of \citet{menking_law-dome_2025}"
    ),
    "input-data-law-dome-rubino": (
        "provided the CSIRO Law Dome ice core 2000-year CO_2, CH_4 and N_2O "
        "data collection"
        # of \citet{rubino_law-dome-dataset_2019}"
    ),
    "input-data-scripps-co2": (
        "provided the Scripps CO_2 Program's merged Law Dome ice core, "
        r"Mauna Loa and South Pole CO_2 record"
        # \citep{keeling_co2-exchanges_2001}"
    ),
    "input-data-trudinger-pfcs": (
        r"provided CF_4, C_2F_6 and C_3F_8 histories"
        # of \citet{trudinger_pfcs_2016}"
    ),
    "input-data-droste-pfcs": (
        r"provided cC_4F_8, C_4F_10, C_5F_12, C_6F_14 and C_7F_16 data "
        "(from Cape Grim and Tacolneston air samples)"
        # r"of \citet{droste_pfcs_2020} "
    ),
    "input-data-wmo-2022": (
        "provided the mixing ratios of the ozone-depleting substances "
        r"from Chapter 7 of the WMO (2022) Scientific Assessment of Ozone Depletion"
        # r"\citep{wmo_ozone-ch7_2022}"
    ),
    "input-data-velders-hfcs": (
        r"provided HFC input"
        # of \citet{velders_hfcs_2022}"
    ),
    "input-data-western-hcfcs": (
        r"provided HCFC input"
        # of \citet{western_hcfcs_2024}"
    ),
    "input-emissions-preparation": (
        "contributed to preparation of the input emissions "
        "with which the simple climate model MAGICC was run"
    ),
    # "validation": "tested and validated the dataset",
    "figures-tables": "produced the figures and tables",
    # "figures-tables": "produced the figures and tables",
    # "figures-design": "designed the figures",
    "original-draft": "wrote the original draft",
    "review-and-editing": "contributed to review and editing",
    "funding-acquisition": "acquired funding",
}
"""What each contribution tag is written as

The tags are what authors can list under `contributions` in the metadata file.
The text follows the initials of the authors who made the contribution,
so it has to read correctly after one author, several authors and "All authors".
The sentences are written in the order the tags are given here.
"""

ALL_AUTHORS = "All authors"
"""What is written, instead of their initials, when every author made a contribution"""


def join_initials(initials: Sequence[str]) -> str:
    """
    Join authors' initials into a list which can start a sentence

    Parameters
    ----------
    initials
        The authors' initials

    Returns
    -------
    :
        The initials, e.g. "ZN", "ZN and MM", "ZN, MM and MP"
    """
    if len(initials) == 1:
        return initials[0]

    return f"{', '.join(initials[:-1])} and {initials[-1]}"


def get_authors_by_contribution(
    authors: Sequence[dict[str, Any]],
) -> dict[str, list[str]]:
    """
    Get the authors who made each contribution

    Parameters
    ----------
    authors
        The authors, as the metadata file gives them

        Each needs `contributions` (a list of the tags in [CONTRIBUTIONS][])
        and `contributions_initials` (how the author is written in the statement).

    Returns
    -------
    :
        The initials of the authors who made each contribution, in the authors' order

        The contributions are in the order of [CONTRIBUTIONS][].
        A contribution which no author made is left out.

    Raises
    ------
    ValueError
        An author has no contributions, no initials or the same initials as another,
        or lists a contribution which is not in [CONTRIBUTIONS][]
    """
    res: dict[str, list[str]] = {tag: [] for tag in CONTRIBUTIONS}
    seen_initials: dict[str, str] = {}
    for author in authors:
        name = f"{author['given_name']} {author['surname']}"

        for key in ("contributions", "contributions_initials"):
            if not author.get(key):
                msg = f"{name} has no `{key}` in the metadata file"
                raise ValueError(msg)

        initials = author["contributions_initials"]
        if initials in seen_initials:
            msg = (
                f"{name} and {seen_initials[initials]} "
                f"have the same `contributions_initials`: {initials!r}"
            )
            raise ValueError(msg)

        seen_initials[initials] = name

        for tag in author["contributions"]:
            if tag not in CONTRIBUTIONS:
                msg = (
                    f"{name} has an unknown contribution: {tag!r}. "
                    f"Expected one of: {', '.join(CONTRIBUTIONS)}"
                )
                raise ValueError(msg)

            if initials in res[tag]:
                msg = f"{name} lists the contribution {tag!r} more than once"
                raise ValueError(msg)

            res[tag].append(initials)

    return {tag: initials for tag, initials in res.items() if initials}


def get_author_contributions(authors: Sequence[dict[str, Any]]) -> str:
    """
    Get the author contributions statement

    Parameters
    ----------
    authors
        The authors, as the metadata file gives them
        (see [get_authors_by_contribution][])

    Returns
    -------
    :
        The statement, one sentence per contribution, one sentence per line
    """
    sentences = []
    for tag, initials in get_authors_by_contribution(authors).items():
        who = (
            ALL_AUTHORS
            if len(authors) > 1 and len(initials) == len(authors)
            else join_initials(initials)
        )
        sentences.append(f"{who} {CONTRIBUTIONS[tag]}.")

    return "\n".join(sentences)


def generate_author_contributions(outfile: Path, metadata_file: Path) -> Path:
    """
    Generate the author contributions statement

    This is always re-generated,
    so that changes to the metadata file are always picked up.

    Parameters
    ----------
    outfile
        File in which to write the statement

    metadata_file
        The manuscript's metadata file, which holds the authors and their contributions

    Returns
    -------
    :
        `outfile`
    """
    with open(metadata_file, "rb") as fh:
        metadata = tomllib.load(fh)

    outfile.parent.mkdir(exist_ok=True, parents=True)
    logger.info(f"Writing {outfile}")
    outfile.write_text(get_author_contributions(metadata["authors"]) + "\n")

    return outfile
