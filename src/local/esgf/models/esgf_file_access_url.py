"""
ESGF file access URL models
"""
# Don't use this, it breaks sqlmodel
# from __future__ import annotations

from typing import TYPE_CHECKING, Optional

from pydantic import ConfigDict
from sqlmodel import Field, Relationship, SQLModel

if TYPE_CHECKING:
    from local.esgf.models.esgf_file import ESGFFileDB


class ESGFFileAccessURLBase(SQLModel):
    """
    File access URL from ESGF
    """

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


class ESGFFileAccessURLDB(ESGFFileAccessURLBase, table=True):
    """
    File access URL from ESGF

    Database model
    """

    id: int | None = Field(default=None, primary_key=True)

    # Would move esgf_file_id onto base if we had a 'create' class
    esgf_file_id: int | None = Field(default=None, foreign_key="esgffiledb.id")
    esgf_file: Optional["ESGFFileDB"] = Relationship(
        back_populates="esgf_file_access_urls"
    )
    """
    File to which this access URL belongs
    """


class ESGFFileAccessURL(ESGFFileAccessURLBase):
    """
    File access URL from ESGF

    Data model i.e. the one that users will likely want to interact with
    """

    model_config = ConfigDict(extra="forbid")

    def to_db_model(self) -> ESGFFileAccessURLDB:
        # Can't just use model_validate because of cross-references
        db_model_init_kwargs = {
            model_field: getattr(self, model_field)
            for model_field in ESGFFileAccessURLBase.model_fields
        }
        res = ESGFFileAccessURLDB(**db_model_init_kwargs)

        return res
