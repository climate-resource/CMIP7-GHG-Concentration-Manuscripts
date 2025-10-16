"""
Support for querying the input4MIPs project
"""

from __future__ import annotations

# TODO: move this into mappings rather than search
MAPPING_FROM_GENERAL_TERMS = {
    "project": "project",
    "variable": "variable",
    "grid": "grid_label",
    "time_sampling": "frequency",
    "cmip_era": "mip_era",
    "source_id": "source_id",
    # Imperfect as versioning is done differently across different projects.
    "source_version": "source_version",
    "institution_id": "institution_id",
    "target_mip": "target_mip",
}
