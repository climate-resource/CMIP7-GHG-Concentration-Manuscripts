"""
Exporting of components as markdown
"""

from __future__ import annotations

from typing import Any

from pylatexenc.latex2text import LatexNodes2Text

latex_encoder = LatexNodes2Text()


def authors_to_markdown(metadata: dict[str, Any]) -> str:
    """
    Convert authors list to markdown

    (So we can drop straight in google docs, for example)
    """
    author_l = []
    affils_order = []
    for author in metadata["authors"]:
        author_affils_t = ",".join(author["affiliations"])
        author_t = latex_encoder.latex_to_text(
            f"{author['given_name']} {author['surname']}<sup>{author_affils_t}</sup>"
        )

        author_l.append(author_t)

        for affil in author["affiliations"]:
            if affil not in affils_order:
                affils_order.append(affil)

    affils_l = []
    for affil in sorted(affils_order):
        affil_t = f"<sup>{affil}</sup>{latex_encoder.latex_to_text(metadata['affiliations'][affil])}"  # noqa: E501
        affils_l.append(affil_t)

    res_authors = ", ".join(author_l)
    res_affiliations = "\n\n".join(affils_l)

    res = "\n\n\n".join((res_authors, res_affiliations))

    return res
