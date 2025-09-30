"""
Support for downloading from ESGF
"""

from __future__ import annotations

import concurrent.futures
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
                if v.service_name == "HTTPServer" and "dkrz" in v.url
            ]
            if len(access_url) != 1:
                raise NotImplementedError
            access_url = access_url[0]

            # # Very cool trick to get the header alone, thanks Bouwe
            # # https://github.com/Climate-REF/climate-ref/issues/212#issuecomment-3347484978
            # import netCDF4
            # with netCDF4.Dataset(f"{access_url.url}#bytes") as ds:
            #     ds.product

            future = pool.submit(
                download_file_parallel_progress_helper,
                access_url.url,
                out_path=Path(esgf_file.esgf_file_local.path),
                thread_positions=thread_positions,
            )
            futures.append(future)

        iterator_results = concurrent.futures.as_completed(futures)

        res_l = [future.result() for future in iterator_results]

    return tuple(res_l)


def download_file_parallel_progress_helper(
    url: str,
    out_path: Path,
    thread_positions: dict[int, int],
    progress: bool = False,
    max_desc_width: int = 100,
) -> Path:
    # TODO: fix up hard-coded widths and chunk sizes
    thread_id = threading.get_ident()
    if thread_id not in thread_positions:
        thread_positions[thread_id] = len(thread_positions)

    # No-one knows why this is needed, but it is in jupyter notebooks
    print(end=" ")

    if len(url) > max_desc_width:
        desc = f"{url[:30]}...{url[-70:]}"
    else:
        desc = url

    out_path.parent.mkdir(exist_ok=True, parents=True)
    with open(out_path, "wb") as fh, httpx.stream("GET", url) as request:
        with tqdm.auto.tqdm(
            desc=desc,
            # TODO: think about this, might be a better way for ESGF
            total=int(request.headers.get("content-length"), 0),
            miniters=1,
            unit="B",
            unit_scale=True,
            unit_divisor=2**10,
            # unit_divisor=2**20,
            position=thread_positions[thread_id],
            leave=False,
        ) as pbar:
            for chunk in request.iter_bytes(chunk_size=2**12):
                fh.write(chunk)
                pbar.update(len(chunk))

    return out_path
