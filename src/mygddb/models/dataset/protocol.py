"""
Protocol for the dataset model definition
"""

from __future__ import annotations

from typing import Any, Protocol


class DatasetLike(Protocol):
    """
    Dataset-like

    Probably too general to actually be useful.
    I have written this more to try and help my own mental model.
    """

    metadata: Any
    """
    Metadata associated with the dataset
    """

    def ensure_available_locally(self) -> None:
        """
        Ensure that the dataset is available locally
        """

    # Helpful concepts but unhelpful because they can't be made concreete
    def load_data(self) -> None:
        """
        Load data
        """
        # Makes more sense to do this as
        # `DataseLoader().load_dataset(dataset: DatasetLike)`
        # so the details of how to load the data
        # ()
