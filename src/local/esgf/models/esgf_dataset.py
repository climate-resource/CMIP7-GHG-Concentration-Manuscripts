"""
ESGF dataset models
"""
# Don't use this, it breaks sqlmodel
# from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path
from typing import TYPE_CHECKING, Any, Optional, TypeVar

from pydantic import ConfigDict
from sqlmodel import Field, Relationship, SQLModel

from local.esgf.download import download_files_parallel_progress
from local.esgf.models.esgf_dataset_local import ESGFDatasetLocal, ESGFDatasetLocalDB
from local.esgf.models.esgf_file import ESGFFile, ESGFFileDB, to_esgf_files
from local.esgf.models.esgf_file_local import ESGFFileLocal
from local.esgf.models.esgf_raw_metadata import (
    MAPPING_TO_GENERAL_TERMS as MAPPING_TO_GENERAL_TERMS_FROM_RAW,
)
from local.esgf.models.esgf_raw_metadata import ESGFRawMetadata

if TYPE_CHECKING:
    from local.esgf.models.esgf_dataset_local import ESGFDatasetLocalDB
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

    esgf_dataset_local: Optional["ESGFDatasetLocalDB"] = Relationship(
        back_populates="esgf_dataset"
    )
    """
    Local version of this dataset
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

    esgf_dataset_local: Optional["ESGFDatasetLocal"] = None
    """
    Local version of this dataset
    """

    @classmethod
    def from_esgf_file_records(
        cls, esgf_file_records: Iterable[dict[str, Any]]
    ) -> "ESGFDataset":
        """
        Initialise from file records retrieved via the ESGF search API

        Parameters
        ----------
        esgf_file_records
            File records from which to initialise

            These should all apply to the same dataset.

        Returns
        -------
        :
            Initialised dataset
        """
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

    def ensure_all_files_available_locally(self) -> tuple[Path, ...]:
        """
        Ensure that all files in the dataset are available locally

        Returns
        -------
        :
            Local paths
        """
        local_files_available = []
        to_download = []
        for esgf_file in self.esgf_files:
            if esgf_file.esgf_file_local is None:
                msg = f"No local path specified for {esgf_file=}"
                raise AssertionError(msg)

            lp = Path(esgf_file.esgf_file_local.path)
            if lp.exists():
                local_files_available.append(lp)
            else:
                to_download.append(esgf_file)

        local_files_downloaded = (
            download_files_parallel_progress(to_download) if to_download else []
        )

        local_files = tuple(
            (
                *local_files_available,
                *local_files_downloaded,
            )
        )

        return local_files

    def set_local_files_root_dir(self, local_files_root_dir: Path) -> ESGFDatasetDB:
        """
        Set the local files' root directory

        This sets the local file paths if they're not already set,
        or updates them if they are set

        Parameters
        ----------
        local_files_root_dir
            Root directory to use for local files

        Returns
        -------
        :
            Updated instance of `self`
        """
        esgf_files_local = []
        for esgf_file in self.esgf_files:
            local_file = ESGFFileLocal(
                path=str(local_files_root_dir / esgf_file.path_esgf)
            )
            esgf_files_local.append(local_file)
            esgf_file.esgf_file_local = local_file

        esgf_dataset_local = ESGFDatasetLocal(esgf_files_local=esgf_files_local)
        self.esgf_dataset_local = esgf_dataset_local

        return self

    def to_db_model(self) -> ESGFDatasetDB:
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
            for model_field in ESGFDatasetBase.model_fields
        }
        db_model_init_kwargs["esgf_files"] = [v.to_db_model() for v in self.esgf_files]
        db_model_init_kwargs["esgf_raw_metadata"] = (
            self.esgf_raw_metadata.to_db_model()
            if self.esgf_raw_metadata is not None
            else None
        )
        # TODO: Hmmm I don't like having to break the pattern here.
        # Think about whether the links are correct
        # (do we need ESGFDatasetLocal? Or we can we just use local files?)
        esgf_files_local = [
            v.esgf_file_local
            for v in db_model_init_kwargs["esgf_files"]
            if v.esgf_file_local is not None
        ]
        if esgf_files_local:
            db_model_init_kwargs["esgf_dataset_local"] = ESGFDatasetLocalDB(
                esgf_files_local=esgf_files_local
            )

        res = ESGFDatasetDB(**db_model_init_kwargs)

        return res
