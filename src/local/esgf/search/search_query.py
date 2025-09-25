"""
General search query interface

The intent here is to provide a common interface for searching
and let this interface handle the details of translating to
CMIP era specific searches.
This is a deliberate abstraction away from ESGF.
This makes it harder for users to find information,
but easier to use as you don't need (in theory)
to understand all the mappings yourself
(you only need to learn one vocabulary, not multiple).
"""

from __future__ import annotations

from enum import StrEnum

from attrs import asdict, define

from local.esgf.search.input4MIPs import (
    MAPPING_FROM_GENERAL_TERMS as MAPPING_FROM_GENERAL_TERMS_INPUT4MIPS,
)
from local.esgf.search.query import query_esgf
from local.esgf.search.search_result import SearchResult


class KnownIndexNode(StrEnum):
    """
    Known index nodes for searching ESGF
    """

    ORNL = "https://esgf-node.ornl.gov/esgf-1-5-bridge"
    """
    Oakridge national lab
    """


@define
class SearchQuery:
    """
    ESGF search query

    Searches across facets use AND logic.
    Searches within facets use OR logic.
    If you want OR logic across facets or AND logic within facets
    then you will need multiple [SearchQuery][]'s.

    Note that this is an abstraction.
    The actual query must be transformed carefully
    as the names of different search terms have changed over time.
    """

    project: str | tuple[str, ...]
    """
    ESGF project(s) to search within
    """

    variable: str | tuple[str, ...] | None = None
    """
    Name of the variable(s) to search for
    """

    grid: str | tuple[str, ...] | None = None
    """
    Grid label(s) to search for
    """
    # This could get super messy as grid names have changed over time,
    # might have to introduce simplified names for grids
    # (or just don't use grid names in your search)

    time_sampling: str | tuple[str, ...] | None = None
    """
    Time sampling(s) to search for
    """

    cmip_era: str | tuple[str, ...] | None = None
    """
    CMIP era(s) in which to search for data
    """

    source_id: str | tuple[str, ...] | None = None
    """
    Source ID(s) to search for
    """

    def to_input4MIPs_terms(self) -> dict[str, str | tuple[str, ...]]:
        input4MIPs_terms = {
            MAPPING_FROM_GENERAL_TERMS_INPUT4MIPS[k]: v
            for k, v in asdict(self).items()
            if v is not None
        }

        return input4MIPs_terms

    def get_results(self, index_node: str | KnownIndexNode) -> tuple[SearchResult, ...]:
        projects = (self.project,) if isinstance(self.project, str) else self.project

        for project in projects:
            # TODO: consider moving this out
            if project == "input4MIPs":
                specific_search_terms = self.to_input4MIPs_terms()
            else:
                raise NotImplementedError(project)

            raw_query_response = query_esgf(
                endpoint=str(index_node),
                query_terms=specific_search_terms,
                # TODO: support these configuration options
                # format=format,
                # distrib=distrib,
                # limit=limit,
            )
            breakpoint()
