"""
ESGF file models
"""
# Don't use this, it breaks sqlmodel
# from __future__ import annotations

from collections.abc import Iterable
from typing import TYPE_CHECKING, Any, Optional

from pydantic import ConfigDict
from sqlmodel import Field, Relationship, SQLModel

from local.esgf.models.esgf_file_access_url import ESGFFileAccessURL
from local.esgf.models.esgf_file_local import ESGFFileLocal

if TYPE_CHECKING:
    from local.esgf.models.esgf_dataset import ESGFDatasetDB
    from local.esgf.models.esgf_file_access_url import ESGFFileAccessURLDB
    from local.esgf.models.esgf_file_local import ESGFFileLocalDB


class ESGFFileBase(SQLModel):
    """
    File on ESGF
    """

    path_esgf: str
    """
    Filepath according to ESGF
    """


class ESGFFileDB(ESGFFileBase, table=True):
    """
    File on ESGF

    Database model
    """

    id: int | None = Field(default=None, primary_key=True)

    esgf_dataset_id: int = Field(foreign_key="esgfdatasetdb.id")
    esgf_dataset: "ESGFDatasetDB" = Relationship(back_populates="esgf_files")
    """
    Dataset to which this file belongs
    """

    esgf_file_access_urls: list["ESGFFileAccessURLDB"] = Relationship(
        back_populates="esgf_file"
    )
    """
    Access URLs available for this file
    """

    esgf_file_local_id: int | None = Field(foreign_key="esgffilelocaldb.id")
    esgf_file_local: "ESGFFileLocalDB" = Relationship(back_populates="esgf_file")
    """
    Local version of this file
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

    # No esgf_dataset to avoid circularity issues

    esgf_file_access_urls: list["ESGFFileAccessURL"] = Field(default_factory=list)
    """
    Access URLs available for this file
    """

    esgf_file_local: Optional["ESGFFileLocal"] = None
    """
    Local version of this file
    """

    def to_db_model(self) -> ESGFFileDB:
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
            for model_field in ESGFFileBase.model_fields
        }
        db_model_init_kwargs["esgf_file_access_urls"] = [
            v.to_db_model() for v in self.esgf_file_access_urls
        ]
        db_model_init_kwargs["esgf_file_local"] = (
            self.esgf_file_local.to_db_model()
            if self.esgf_file_local is not None
            else None
        )
        res = ESGFFileDB(**db_model_init_kwargs)

        return res


def to_esgf_files(esgf_file_records: Iterable[dict[str, Any]]) -> tuple[ESGFFile, ...]:
    """
    Convert ESGF file records to [ESGFFile][]'s

    Parameters
    ----------
    esgf_file_records
        ESGF file records


    Returns
    -------
    :
        [ESGFFile][]'s
    """
    file_access_urls_grouped = {}
    file_paths_by_id = {}
    for result_d in esgf_file_records:
        # Want the ID without the node
        our_file_id = result_d["instance_id"]
        if our_file_id not in file_access_urls_grouped:
            file_access_urls_grouped[our_file_id] = []

        # Best way to get file path I can think of,
        # I really hope this convention doesn't break...
        file_path = our_file_id.replace(".", "/").replace("/nc", ".nc")
        if our_file_id not in file_paths_by_id:
            file_paths_by_id[our_file_id] = file_path
        elif file_path != file_paths_by_id[our_file_id]:
            msg = (
                f"{file_path} != {file_paths_by_id[our_file_id]}. {esgf_file_records=}"
            )
            raise AssertionError(msg)

        for access_url in result_d["url"]:
            # https://esgf.github.io/esg-search/ESGF_Search_RESTful_API.html#access-urls
            url, mime_type, service_name = access_url.split("|")
            esgf_file_access_url = ESGFFileAccessURL(
                url=url,
                mime_type=mime_type,
                service_name=service_name,
            )
            file_access_urls_grouped[our_file_id].append(esgf_file_access_url)

    esgf_files = tuple(
        ESGFFile(
            path_esgf=file_paths_by_id[file_id],
            esgf_file_access_urls=file_access_urls,
            esgf_file_local=None,
        )
        for file_id, file_access_urls in file_access_urls_grouped.items()
    )

    return esgf_files
