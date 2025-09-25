"""
Support for querying the input4MIPs project
"""

from __future__ import annotations

from attrs import define

MAPPING_FROM_GENERAL_TERMS = {
    "project": "project",
    "variable": "variable_id",
    "grid": "grid_label",
    "time_sampling": "frequency",
    "cmip_era": "mip_era",
    "source_id": "source_id",
}


# Not sure if this is needed
@define
class SearchQueryInput4MIPs:
    """
    ESGF search query specific to searching the input4MIPs project
    """
