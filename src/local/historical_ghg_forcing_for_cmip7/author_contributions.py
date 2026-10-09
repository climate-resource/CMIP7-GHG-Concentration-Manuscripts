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
