"""
ESGF dataset models
"""
# Don't use this, it breaks sqlmodel
# from __future__ import annotations

from collections.abc import Iterable
from typing import TYPE_CHECKING, Any, Optional, TypeVar

from pydantic import ConfigDict
from sqlmodel import Field, Relationship, SQLModel

from local.esgf.models.esgf_file import ESGFFile, ESGFFileDB, to_esgf_files
from local.esgf.models.esgf_raw_metadata import (
    MAPPING_TO_GENERAL_TERMS as MAPPING_TO_GENERAL_TERMS_FROM_RAW,
)
from local.esgf.models.esgf_raw_metadata import ESGFRawMetadata

if TYPE_CHECKING:
    from local.esgf.models.esgf_raw_metadata import ESGFRawMetadataDB

T = TypeVar("T")


class ESGFDatasetBase(SQLModel):
    """
    Dataset on ESGF
    """

    # No version concept yet as version
    # is highly variable across projects
    # depending on how they use source id
    # vs. drs version
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


# TODO : think about uniqueness constraints
class ESGFDatasetDB(ESGFDatasetBase, table=True):
    """
    Dataset on ESGF

    Database model
    """

    id: int | None = Field(default=None, primary_key=True)

    esgf_files: list["ESGFFileDB"] = Relationship(back_populates="esgf_dataset")
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


def coerce_to_single_value(inv: T | list[T]) -> T:
    """
    Coerce the input to a single value

    Parameters
    ----------
    inv
        Input value

    Returns
    -------
    :
        `inv` converted to a single value

    Raises
    ------
    AssertionError
        `inv` has a length greater than one
    """
    if isinstance(inv, list):
        if len(inv) != 1:
            raise AssertionError(inv)

        return inv[0]

    return inv


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

    @classmethod
    def from_esgf_file_records(
        cls, esgf_file_records: Iterable[dict[str, Any]]
    ) -> "ESGFDataset":
        init_kwargs = {}
        init_kwargs_raw_metadata = {}

        all_record_keys = set(k for r in esgf_file_records for k in r)
        for key in all_record_keys:
            if key not in MAPPING_TO_GENERAL_TERMS_FROM_RAW:
                continue

            general_term = MAPPING_TO_GENERAL_TERMS_FROM_RAW[key]
            esgf_record_values = []
            for esgf_file_record in esgf_file_records:
                try:
                    esgf_record_values.append(
                        coerce_to_single_value(esgf_file_record[key])
                    )
                except AssertionError as exc:
                    msg = (
                        f"For the following record, {key=} has more than one value. "
                        f"Values for {key}: {esgf_file_record[key]}. "
                        f"{esgf_file_record=}"
                    )
                    raise AssertionError(msg) from exc

            if len(set(esgf_record_values)) != 1:
                msg = (
                    f"The component records do not agree on the value of {key}. "
                    f"Values found: {esgf_record_values}. "
                    f"Raw responses from ESGF {esgf_file_records}"
                )
                raise AssertionError(msg)

            init_kwargs_raw_metadata[key] = esgf_record_values[0]
            init_kwargs[general_term] = esgf_record_values[0]

        init_kwargs["esgf_raw_metadata"] = ESGFRawMetadata.model_validate(
            init_kwargs_raw_metadata
        )
        esgf_files = to_esgf_files(esgf_file_records)
        init_kwargs["esgf_files"] = esgf_files
        res = cls.model_validate(init_kwargs)

        return res

    def to_db_model(self) -> ESGFDatasetDB:
        # Can't just use model_validate because of cross-references
        db_model_init_kwargs = {
            model_field: getattr(self, model_field)
            for model_field in ESGFDatasetBase.model_fields
        }
        db_model_init_kwargs["esgf_files"] = [v.to_db_model() for v in self.esgf_files]
        db_model_init_kwargs["esgf_raw_metadata"] = (
            self.esgf_raw_metadata.to_db_model()
            if self.esgf_raw_metadata is not None
            else None
        )
        res = ESGFDatasetDB(**db_model_init_kwargs)

        return res
