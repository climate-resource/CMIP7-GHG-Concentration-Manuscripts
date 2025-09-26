"""
Model of a dataset on ESGF
"""

# Don't use this, it breaks sqlmodel
# from __future__ import annotations

from typing import Union

from sqlmodel import Field, Relationship, SQLModel


class ESGFDataset(SQLModel, table=True):
    """
    Representation of a dataset on ESGF

    Not a perfect mapping as some names are altered
    (on purpose, to make representation simpler).
    """

    id: int | None = Field(default=None, primary_key=True)

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

    esgf_files: list["ESGFFile"] = Relationship(back_populates="esgf_dataset")

    esgf_raw_metadata: Union["ESGFRawMetadata", None] = Relationship(
        back_populates="esgf_dataset"
    )


class ESGFFile(SQLModel, table=True):
    """
    File on ESGF
    """

    id: int | None = Field(default=None, primary_key=True)

    # Download locations or similar

    # Could add in time range here

    esgf_dataset_id: int = Field(foreign_key="esgfdataset.id")
    esgf_dataset: ESGFDataset = Relationship(back_populates="esgf_files")
    """
    Dataset to which this file belongs
    """

    esgf_file_access_urls: list["ESGFFileAccessURL"] = Relationship(
        back_populates="esgf_file"
    )
    """
    Access URLs available for this file
    """


class ESGFFileAccessURL(SQLModel, table=True):
    """
    File access URL from ESGF
    """

    id: int | None = Field(default=None, primary_key=True)

    url: str
    """
    URL
    """

    mime_type: str
    """
    Media type, see e.g. https://developer.mozilla.org/en-US/docs/Web/HTTP/Guides/MIME_types
    """

    service_name: str
    """
    Name of the service according to ESGF
    """

    esgf_file_id: int = Field(foreign_key="esgffile.id")
    esgf_file: ESGFFile = Relationship(back_populates="esgf_file_access_urls")
    """
    File to which this access URL belongs
    """


class ESGFRawMetadata(SQLModel, table=True):
    """
    Raw metadata with names as used on ESGF
    """

    id: int | None = Field(default=None, primary_key=True)

    variable: str
    experiment: str

    esgf_dataset_id: int = Field(foreign_key="esgfdataset.id")
    esgf_dataset: ESGFDataset = Relationship(back_populates="esgf_raw_metadata")
    """
    Dataset to which this raw metadata applies
    """
