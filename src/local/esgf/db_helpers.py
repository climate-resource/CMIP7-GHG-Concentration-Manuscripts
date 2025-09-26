"""
Database helpers
"""

from pathlib import Path

import sqlalchemy.engine.base
from sqlmodel import SQLModel, create_engine

# Have to import so the tables will be created and therefore can be registered
from local.esgf.models import (
    ESGFDatasetDB,  # noqa: F401
    ESGFFileAccessURLDB,  # noqa: F401
    ESGFFileDB,  # noqa: F401
    ESGFRawMetadataDB,  # noqa: F401
)


def create_all_tables(engine: sqlalchemy.engine.base.Engine) -> None:
    """
    Create all tables based on the models that have been registered with [sqlmodel][]

    Parameters
    ----------
    engine
        Engine to use for interacting with the database
    """
    SQLModel.metadata.create_all(engine)


def get_sqlite_engine(
    sqlite_db_path: Path,
) -> sqlalchemy.engine.base.Engine:
    """
    Get SQLite engine

    Parameters
    ----------
    sqlite_db_path
        Path to the SQLite database

    Returns
    -------
    :
        Engine to use for interacting with the database
    """
    sqlite_url = f"sqlite:///{sqlite_db_path.as_posix()}"
    engine = create_engine(
        sqlite_url,
        # echo=True,
    )

    return engine
