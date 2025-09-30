"""
Collection of ESGF datasets

This does not get used in the database, hence does not appear there
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from attrs import define

if TYPE_CHECKING:
    from local.esgf.models import ESGFDataset


@define
class ESGFDatasetCollection:
    """
    Collection of [ESGFDataset][local.esgf.models.ESGFDataset]
    """

    esgf_datasets: tuple[ESGFDataset, ...]
    """
    ESGF datasets in the collection
    """

    # def to_df(self) -> pd.DataFrame:
    #     raise NotImplementedError
