"""
Local ESGF dataset models

When we say "local ESGF dataset",
what we mean is a local (i.e. downloaded) dataset
that corresponds to a dataset record available via ESGF.
"""
# Don't use this, it breaks sqlmodel
# from __future__ import annotations

from typing import TYPE_CHECKING

from pydantic import ConfigDict
from sqlmodel import Field, Relationship, SQLModel

if TYPE_CHECKING:
    from local.esgf.models.esgf_dataset import ESGFDatasetDB
    from local.esgf.models.esgf_file_local import ESGFFileLocal, ESGFFileLocalDB


class ESGFDatasetLocalBase(SQLModel):
    """
    Local dataset corresponding to an ESGF dataset
    """


# TODO : think about uniqueness constraints
class ESGFDatasetLocalDB(ESGFDatasetLocalBase, table=True):
    """
    Local dataset corresponding to an ESGF dataset

    Database model
    """

    id: int | None = Field(default=None, primary_key=True)

    esgf_dataset_id: int = Field(foreign_key="esgfdatasetdb.id")
    esgf_dataset: list["ESGFDatasetDB"] = Relationship(
        back_populates="esgf_dataset_local"
    )
    """
    ESGF dataset to which this local dataset corresponds
    """

    esgf_files_local: list["ESGFFileLocalDB"] = Relationship(
        back_populates="esgf_dataset_local"
    )
    """
    Local files that make up this dataset
    """


class ESGFDatasetLocalNoLinks(ESGFDatasetLocalBase):
    """
    Local dataset corresponding to an ESGF dataset

    Data model without any links.
    As a user, you will probably prefer [ESGFDatasetLocal][].
    """

    model_config = ConfigDict(extra="forbid")


class ESGFDatasetLocal(ESGFDatasetLocalNoLinks):
    """
    Local dataset corresponding to an ESGF dataset

    Data model i.e. the one that users will likely want to interact with
    """

    esgf_files_local: list["ESGFFileLocal"] = Field(default_factory=list)
    """
    Local files that make up this dataset
    """

    def to_db_model(self) -> ESGFDatasetLocalDB:
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
            for model_field in ESGFDatasetLocalBase.model_fields
        }
        db_model_init_kwargs["esgf_files_local"] = [
            v.to_db_model() for v in self.esgf_files_local
        ]
        res = ESGFDatasetLocalDB(**db_model_init_kwargs)

        return res
