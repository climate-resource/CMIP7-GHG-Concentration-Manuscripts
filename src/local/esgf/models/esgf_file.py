"""
ESGF file models
"""
# Don't use this, it breaks sqlmodel
# from __future__ import annotations

from typing import TYPE_CHECKING

from pydantic import ConfigDict
from sqlmodel import Field, Relationship, SQLModel

if TYPE_CHECKING:
    from local.esgf.models.esgf_file_access_url import (
        ESGFFileAccessURL,
        ESGFFileAccessURLDB,
    )


class ESGFFileBase(SQLModel):
    """
    File on ESGF
    """


class ESGFFileDB(ESGFFileBase, table=True):
    """
    File on ESGF

    Database model
    """

    id: int | None = Field(default=None, primary_key=True)

    # esgf_dataset_id: int = Field(foreign_key="esgfdataset.id")
    # esgf_dataset: ESGFDataset = Relationship(back_populates="esgf_files")
    # """
    # Dataset to which this file belongs
    # """

    esgf_file_access_urls: list["ESGFFileAccessURLDB"] = Relationship(
        back_populates="esgf_file"
    )
    """
    Access URLs available for this file
    """


class ESGFFileNoLinks(ESGFFileBase):
    """
    File on ESGF

    Data model without any links.
    As a user, you will probably prefer [ESGFFile][].
    """

    model_config = ConfigDict(extra="forbid")


class ESGFFile(ESGFFileNoLinks):
    """
    File on ESGF

    Data model i.e. the one that users will likely want to interact with
    """

    esgf_file_access_urls: list["ESGFFileAccessURL"] = Field(default_factory=list)
    """
    Access URLs available for this file
    """
