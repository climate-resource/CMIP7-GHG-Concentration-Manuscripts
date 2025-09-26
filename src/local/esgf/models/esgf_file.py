"""
ESGF file models
"""
# Don't use this, it breaks sqlmodel
# from __future__ import annotations

from collections.abc import Iterable
from typing import TYPE_CHECKING, Any

from pydantic import ConfigDict
from sqlmodel import Field, Relationship, SQLModel

from local.esgf.models.esgf_file_access_url import ESGFFileAccessURL

if TYPE_CHECKING:
    from local.esgf.models.esgf_dataset import ESGFDatasetDB
    from local.esgf.models.esgf_file_access_url import ESGFFileAccessURLDB


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

    def to_db_model(self) -> ESGFFileDB:
        # Can't just use model_validate because of cross-references
        db_model_init_kwargs = {
            model_field: getattr(self, model_field)
            for model_field in ESGFFileBase.model_fields
        }
        db_model_init_kwargs["esgf_file_access_urls"] = [
            v.to_db_model() for v in self.esgf_file_access_urls
        ]
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
    for result_d in esgf_file_records:
        # Get rid of the node before creating the file ID
        our_file_id = result_d["id"].split("|")[0]
        if our_file_id not in file_access_urls_grouped:
            file_access_urls_grouped[our_file_id] = []

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
        ESGFFile(esgf_file_access_urls=file_access_urls)
        for file_access_urls in file_access_urls_grouped.values()
    )

    return esgf_files
