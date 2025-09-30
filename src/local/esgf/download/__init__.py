"""
Support for downloading from ESGF
"""

from __future__ import annotations

import concurrent.futures
import shutil
import tempfile
import threading
from collections.abc import Iterable
from pathlib import Path
from typing import TYPE_CHECKING

import httpx
import tqdm.auto

if TYPE_CHECKING:
    from local.esgf.models import ESGFFile


# TODO: move to using a database to manage this
# so we can deal with re-queing, downloads that fail halfway through,
# prioritising nodes etc.
# Also makes it easier to decouple parallel and progress bar configuration
def download_files_parallel_progress(
    esgf_files: Iterable[ESGFFile], n_processes: int = 3
) -> tuple[Path, ...]:
    thread_positions = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=n_processes) as pool:
        futures = []
        for esgf_file in esgf_files:
            # TODO: split out node handling rather than hard-coding DKZR
            access_url = [
                v
                for v in esgf_file.esgf_file_access_urls
                if v.service_name == "HTTPServer"
                # and "dkrz" in v.url
                and "globus" not in v.url
            ]
            if len(access_url) < 1:
                # Interesting, can have more than one HTTPServer on the same node
                # TODO: deal with this as part of URL sorting
                raise NotImplementedError(esgf_file.esgf_file_access_urls)

            access_url = access_url[0]

            # # Very cool trick to get the header alone, thanks Bouwe
            # # https://github.com/Climate-REF/climate-ref/issues/212#issuecomment-3347484978
            # import netCDF4
            # with netCDF4.Dataset(f"{access_url.url}#bytes") as ds:
            #     ds.product

            future = pool.submit(
                download_file_parallel_progress_helper,
                url=access_url.url,
                size=esgf_file.size,
                out_path=Path(esgf_file.esgf_file_local.path),
                thread_positions=thread_positions,
            )
            futures.append(future)

        iterator_results = tqdm.auto.tqdm(
            concurrent.futures.as_completed(futures),
            desc="Datasets",
            total=len(futures),
            position=0,
        )

        res_l = [future.result() for future in iterator_results]

    return tuple(res_l)


def download_file_parallel_progress_helper(
    url: str,
    size: int,
    out_path: Path,
    thread_positions: dict[int, int],
    progress: bool = False,
    max_desc_width: int = 100,
) -> Path:
    # TODO: fix up hard-coded widths and chunk sizes
    thread_id = threading.get_ident()
    if thread_id not in thread_positions:
        thread_positions[thread_id] = len(thread_positions) + 1

    # No-one knows why this is needed, but it is in jupyter notebooks
    print(end=" ")

    if len(url) > max_desc_width:
        desc = f"{url[:30]}...{url[-70:]}"
    else:
        desc = url

    with (
        tempfile.TemporaryDirectory() as td,
        httpx.stream("GET", url, follow_redirects=True) as request,
    ):
        tmpf = Path(td) / out_path.name
        with (
            open(tmpf, "wb") as fh,
            tqdm.auto.tqdm(
                desc=desc,
                total=size,
                miniters=1,
                unit="B",
                unit_scale=True,
                unit_divisor=2**10,
                position=thread_positions[thread_id],
                leave=False,
            ) as pbar,
        ):
            for chunk in request.iter_bytes(chunk_size=2**17):
                fh.write(chunk)
                pbar.update(len(chunk))

        # If we got to here, can move the tempfile to our destination file
        out_path.parent.mkdir(exist_ok=True, parents=True)
        shutil.move(tmpf, out_path)

    return out_path
