"""
ESGF file models

When we say "local ESGF file",
what we mean is local (i.e. downloaded) file
that corresponds to a file available via ESGF.
"""
# Don't use this, it breaks sqlmodel
# from __future__ import annotations

from typing import TYPE_CHECKING

from pydantic import ConfigDict
from sqlmodel import Field, Relationship, SQLModel

if TYPE_CHECKING:
    from local.esgf.models.esgf_file import ESGFFileDB


class ESGFFileLocalBase(SQLModel):
    """
    Local file corresponding to an ESGF file
    """

    path: str
    """
    Path to the local dataset
    """


class ESGFFileLocalDB(ESGFFileLocalBase, table=True):
    """
    Local file corresponding to an ESGF file

    Database model
    """

    id: int | None = Field(default=None, primary_key=True)

    esgf_file: "ESGFFileDB" = Relationship(back_populates="esgf_file_local")
    """
    ESGF file to which this local file corresponds
    """


class ESGFFileLocalNoLinks(ESGFFileLocalBase):
    """
    Local file corresponding to an ESGF file

    Data model without any links.
    As a user, you will probably prefer [ESGFFileLocal][].
    """

    model_config = ConfigDict(extra="forbid")


class ESGFFileLocal(ESGFFileLocalNoLinks):
    """
    Local file corresponding to an ESGF file

    Data model i.e. the one that users will likely want to interact with
    """

    # No esgf_dataset_local to avoid circularity

    # No esgf_file to avoid circularity

    def to_db_model(self) -> ESGFFileLocalDB:
        """
        Convert to the database model

        Returns
        -------
        :
            Database model equivalent of `self`
        """
        # Can't just use model_validate because of cross-references
        db_model_init_kwargs = {
            model_field: getattr(self, model_field)
            for model_field in ESGFFileLocalBase.model_fields
        }
        res = ESGFFileLocalDB(**db_model_init_kwargs)

        return res
