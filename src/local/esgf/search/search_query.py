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

from local.esgf.esgf_dataset_collection import ESGFDatasetCollection
from local.esgf.search.input4MIPs import (
    MAPPING_FROM_GENERAL_TERMS as MAPPING_FROM_GENERAL_TERMS_INPUT4MIPS,
)
from local.esgf.search.query import query_esgf


class KnownIndexNode(StrEnum):
    """
    Known index nodes for searching ESGF
    """

    # https://esgf.ceda.ac.uk/thredds/dodsC/esg_cmip6/input4MIPs/CMIP7/CMIP/CR/CR-CMIP-1-0-0/atmos/mon/c4f10/gr1z/v20250228/c4f10_input4MIPs_GHGConcentrations_CMIP_CR-CMIP-1-0-0_gr1z_000101-099912.nc
    CEDA = "https://esgf.ceda.ac.uk/esg-search/search"
    """
    Centre for Environmental Data Analysis
    """

    DKRZ = "https://esgf-data.dkrz.de/esg-search/search"
    """
    Deutches Klimarechenzentrum (German climate computing centre)
    """

    NCI = "https://esgf.nci.org.au/esg-search/search"
    """
    National Computing Infrastructure (Australia)
    """

    ORNL = "https://esgf-node.ornl.gov/esgf-1-5-bridge"
    """
    Oakridge national lab
    """


@define
class SearchQuery:
    """
    ESGF search query

    Searches with multiple values for a given search term
    are deliberately not allowed in this low-level interface.
    The reason is that the servers seem to not implement the logic correctly,
    and we don't want to pass that confusion on to users.
    If you want combination logic across search terms
    then you will need multiple [SearchQuery][]'s.

    Note that this is an abstraction.
    The actual query must be transformed carefully
    as the names of different search terms have changed over time.
    """

    project: str
    """
    ESGF project(s) to search within
    """

    variable: str | None = None
    """
    Name of the variable(s) to search for
    """

    grid: str | None = None
    """
    Grid label(s) to search for
    """
    # This could get super messy as grid names have changed over time,
    # might have to introduce simplified names for grids
    # (or just don't use grid names in your search)

    time_sampling: str | None = None
    """
    Time sampling(s) to search for
    """

    cmip_era: str | None = None
    """
    CMIP era(s) in which to search for data
    """

    source_id: str | None = None
    """
    Source ID(s) to search for
    """

    def to_esgf_seach_terms(self) -> dict[str, str]:
        if self.project == "input4MIPs":
            mapping = MAPPING_FROM_GENERAL_TERMS_INPUT4MIPS
        else:
            raise NotImplementedError(self.project)

        esgf_search_terms = {
            mapping[k]: v for k, v in asdict(self).items() if v is not None
        }

        return esgf_search_terms

    def get_results(
        self,
        index_node: str,
        distrib: bool = True,
        limit: int = 1_000,
    ) -> ESGFDatasetCollection:
        # TODO: docstring
        """
        Should the query be distributed?

        I.e. look at results both on the index node and other nodes?
        """
        """
        Maximum amount of results to retrieve
        """

        esgf_search_terms = self.to_esgf_seach_terms()

        esgf_dataset_collection = query_esgf(
            endpoint=str(index_node),
            query_terms=esgf_search_terms,
            distrib=distrib,
            limit=limit,
        )

        return esgf_dataset_collection
