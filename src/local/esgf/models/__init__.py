"""
Data models

Kept in their own module to allow special resolution of cross-refences etc.
that the user shouldn't have to worry about.

All these models are implemented
following the pattern recommended by sqlmodel and fastapi:
- https://sqlmodel.tiangolo.com/tutorial/fastapi/multiple-models
- https://sqlmodel.tiangolo.com/tutorial/fastapi/relationships/

This leads to some sub-classing, but this seems to be the least worst option
in order to be able to get data validation on the Python API side.
"""

from local.esgf.models.esgf_dataset import ESGFDataset, ESGFDatasetDB
from local.esgf.models.esgf_file import (
    ESGFFile,
    ESGFFileDB,
)
from local.esgf.models.esgf_file_access_url import (
    ESGFFileAccessURL,
    ESGFFileAccessURLDB,
)
from local.esgf.models.esgf_raw_metadata import (
    ESGFRawMetadata,
    ESGFRawMetadataDB,
)

# Ensure cross-links validate properly
ESGFFile.model_rebuild()
ESGFFileAccessURL.model_rebuild()
ESGFFileDB.model_rebuild()
ESGFFileAccessURLDB.model_rebuild()
ESGFRawMetadata.model_rebuild()
ESGFRawMetadataDB.model_rebuild()
ESGFDataset.model_rebuild()
ESGFDatasetDB.model_rebuild()

__all__ = [
    "ESGFDataset",
    "ESGFDatasetDB",
    "ESGFFile",
    "ESGFFileAccessURL",
    "ESGFFileAccessURLDB",
    "ESGFFileDB",
    "ESGFRawMetadata",
    "ESGFRawMetadataDB",
]
