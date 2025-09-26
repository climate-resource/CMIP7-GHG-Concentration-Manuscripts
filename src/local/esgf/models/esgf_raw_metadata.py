"""
ESGF raw metadata models
"""
# Don't use this, it breaks sqlmodel
# from __future__ import annotations

from typing import TYPE_CHECKING, Optional

from pydantic import ConfigDict
from sqlmodel import Field, Relationship, SQLModel

if TYPE_CHECKING:
    from local.esgf.models import ESGFDatasetDB

MAPPING_TO_GENERAL_TERMS = {
    "project": "project",
    "variable_id": "variable",
    "grid_label": "grid",
    "frequency": "time_sampling",
    "mip_era": "cmip_era",
    "source_id": "source_id",
}
"""
Mapping from raw metadata terms to general terms
"""


class ESGFRawMetadataBase(SQLModel):
    """
    ESGF raw metadata
    """

    frequency: str | None = None
    """
    Frequency label of the ESGF record

    A bit of a misnomer, it's actually time sampling.
    Frequency should have units 1 / time,
    but the values here represet values of time e.g. month, year, 3 hr.
    """

    grid_label: str | None = None
    """
    Grid label of the ESGF record
    """

    mip_era: str | None = None
    """
    MIP era to which the record belongs
    """

    project: str
    """
    Project to which the record belongs
    """

    source_id: str | None
    """
    Source ID of the record

    A bit of a fuzzy concept.
    Think of it as a unique ID for the source of the data,
    where source is loosely defined
    (not always institute, not always model, can include some version information).
    """

    variable_id: str
    """
    ID of the variable represented by the ESGF record
    """


class ESGFRawMetadataDB(ESGFRawMetadataBase, table=True):
    """
    ESGF raw metadata

    Database model
    """

    id: int | None = Field(default=None, primary_key=True)

    # Would move esgf_dataset_id onto base if we had a 'create' class
    esgf_dataset_id: int | None = Field(default=None, foreign_key="esgfdatasetdb.id")
    esgf_dataset: Optional["ESGFDatasetDB"] = Relationship(
        back_populates="esgf_raw_metadata"
    )
    """
    Dataset to which this raw metadata belongs
    """


class ESGFRawMetadata(ESGFRawMetadataBase):
    """
    ESGF raw metadata

    Data model i.e. the one that users will likely want to interact with
    """

    model_config = ConfigDict(extra="forbid")

    def to_db_model(self) -> ESGFRawMetadataDB:
        # Can't just use model_validate because of cross-references
        db_model_init_kwargs = {
            model_field: getattr(self, model_field)
            for model_field in ESGFRawMetadataBase.model_fields
        }
        res = ESGFRawMetadataDB(**db_model_init_kwargs)

        return res
