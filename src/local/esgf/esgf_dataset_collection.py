"""
Collection of ESGF datasets

This does not get used in the database, hence does not appear there
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from attrs import define

from local.esgf.download import download_files_parallel_progress

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

    def ensure_all_files_available_locally(self) -> tuple[Path, ...]:
        """
        Ensure that all files in all the datasets are available locally

        Returns
        -------
        :
            Local paths
        """
        local_files = []
        to_download = []
        for esgf_dataset in self.esgf_datasets:
            local_files_split = esgf_dataset.split_local_files_by_availability()

            local_files.extend(local_files_split["available_locally"])
            to_download.extend(local_files_split["not_available_locally"])

        local_files = tuple(
            (*local_files, *download_files_parallel_progress(to_download))
        )

        return local_files

    def set_local_files_root_dir(
        self, local_files_root_dir: Path
    ) -> ESGFDatasetCollection:
        """
        Set the local files' root directory

        This sets the local file paths
        for all values of `self.esgf_datasets` if they're not already set,
        or updates them if they are set.

        Parameters
        ----------
        local_files_root_dir
            Root directory to use for local files

        Returns
        -------
        :
            Updated instance of `self`
        """
        for ed in self.esgf_datasets:
            ed.set_local_files_root_dir(local_files_root_dir)

        return self
