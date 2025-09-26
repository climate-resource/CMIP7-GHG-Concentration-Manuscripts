"""
ESGF dataset models
"""
# Don't use this, it breaks sqlmodel
# from __future__ import annotations

from typing import TYPE_CHECKING, Optional

from pydantic import ConfigDict
from sqlmodel import Field, Relationship, SQLModel

if TYPE_CHECKING:
    from local.esgf.models.esgf_file import ESGFFile
    from local.esgf.models.esgf_raw_metadata import ESGFRawMetadata, ESGFRawMetadataDB


class ESGFDatasetBase(SQLModel):
    """
    Dataset on ESGF
    """

    project: str
    """
    ESGF project
    """

    variable: str
    """
    Name of the variable
    """

    grid: str
    """
    Grid label
    """

    time_sampling: str
    """
    Time sampling
    """

    cmip_era: str
    """
    CMIP era
    """

    source_id: str
    """
    Source ID
    """


class ESGFDatasetDB(ESGFDatasetBase, table=True):
    """
    Dataset on ESGF

    Database model
    """

    id: int | None = Field(default=None, primary_key=True)

    esgf_files: list["ESGFFile"] = Relationship(back_populates="esgf_dataset")
    """
    Files that make up this dataset
    """

    esgf_raw_metadata: Optional["ESGFRawMetadataDB"] = Relationship(
        back_populates="esgf_dataset"
    )
    """
    Raw metadata associated with this dataset
    """


class ESGFDatasetNoLinks(ESGFDatasetBase):
    """
    Dataset on ESGF

    Data model without any links.
    As a user, you will probably prefer [ESGFDataset][].
    """

    model_config = ConfigDict(extra="forbid")


class ESGFDataset(ESGFDatasetNoLinks):
    """
    Dataset on ESGF

    Data model i.e. the one that users will likely want to interact with
    """

    esgf_files: list["ESGFFile"] = Field(default_factory=list)
    """
    Files that make up this dataset
    """

    esgf_raw_metadata: Optional["ESGFRawMetadata"] = None
    """
    Raw metadata associated with this dataset
    """
