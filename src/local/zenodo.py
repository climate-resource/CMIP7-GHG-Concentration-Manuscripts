"""
Helpers for retrieving data from Zenodo

Thin wrappers around [openscm_zenodo][], plus the tarball extraction
we need because Zenodo records generally archive directory trees as tarballs.
"""

from __future__ import annotations

import tarfile
from collections.abc import Collection
from pathlib import Path

from loguru import logger
from openscm_zenodo import ZenodoClient


def download_files(
    record_id: str,
    filenames: Collection[str],
    dest_dir: Path,
) -> list[Path]:
    """
    Download files from a Zenodo record

    Files which have already been downloaded are not downloaded again
    (`openscm_zenodo` compares the MD5 of any existing file
    with the MD5 recorded on the Zenodo record),
    so calling this repeatedly is cheap.

    No Zenodo token is required, as long as the record is public.

    Parameters
    ----------
    record_id
        ID of the Zenodo record from which to download

    filenames
        Names of the files to download, as they appear on the record

    dest_dir
        Directory into which to download the files

    Returns
    -------
        Path to each downloaded file, in the same order as `filenames`
    """
    dest_dir.mkdir(exist_ok=True, parents=True)

    logger.info(
        f"Ensuring {len(filenames)} file(s) "
        f"from Zenodo record {record_id} are in {dest_dir}"
    )
    with ZenodoClient() as client:
        return client.download_files(record_id, dest_dir, filenames=filenames)


def extract_tar_gz(
    tarball: Path,
    dest_dir: Path,
    strip_components: int = 0,
) -> None:
    """
    Extract a gzipped tarball

    Parameters
    ----------
    tarball
        Tarball to extract

    dest_dir
        Directory into which to extract the tarball's contents

    strip_components
        Number of leading path components to strip from each member's name

        This is the equivalent of `tar`'s `--strip-components`.
        It is needed because Zenodo tarballs commonly carry
        a prefix which describes where the files sat in the original run.
    """
    dest_dir.mkdir(exist_ok=True, parents=True)

    logger.info(f"Extracting {tarball} into {dest_dir}")
    with tarfile.open(tarball, "r:gz") as th:
        for member in th.getmembers():
            stripped = Path(*Path(member.name).parts[strip_components:])
            if not stripped.parts:
                # Member is (part of) the stripped prefix, nothing to extract
                continue

            member.name = str(stripped)
            # `filter="data"` refuses members which would escape `dest_dir`.
            # It is the default from Python 3.14, we set it explicitly
            # so the behaviour doesn't depend on the Python version.
            th.extract(member, dest_dir, filter="data")
