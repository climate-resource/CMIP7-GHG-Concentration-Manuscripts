"""
General search result interface

The intent here is to provide a common interface for search results
and let this interface handle the details of translating from the results of
CMIP era specific searches.
This is a deliberate abstraction away from ESGF.
This makes it harder for users to find information,
but easier to use as you don't need (in theory)
to understand all the mappings yourself
(you only need to learn one vocabulary, not multiple).
"""

from __future__ import annotations

from attrs import define


@define
class SearchResult:
    """
    ESGF search result

    Note that this is an abstraction.
    The actual result is transformed carefully
    as the names of different search result terms have changed over time.
    """
