import concurrent.futures
import threading
from collections.abc import Iterable
from pathlib import Path

import httpx
import tqdm.auto


def download_file_parallel_helper(
    url: str,
    out_path: Path,
    thread_positions: dict[int, int],
    progress: bool = False,
) -> Path:
    thread_id = threading.get_ident()
    if thread_id not in thread_positions:
        thread_positions[thread_id] = len(thread_positions)

    # with open("thread-info.txt", "a") as fh:
    #     fh.write(f"Writing from {thread_id=}")
    #     fh.write(f"\t{thread_positions=}")
    #     fh.write(f"\t{thread_positions[thread_id]=}")
    #     fh.write("\n")

    if progress:
        # try:
        #     from tqdm.auto import tqdm
        # except ImportError as exc:
        #     raise MissingOptionalDependencyError(
        #         "dist(..., progress=True)", requirement="tdqm"
        #     ) from exc

        # No-one knows why this is needed, but it is in jupyter notebooks
        print(end=" ")

        if len(url) > 100:
            desc = f"{url[:30]}...{url[-70:]}"
        else:
            desc = url

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

    else:
        raise NotImplementedError(progress)

    return out_path


def download(
    # URLS, out path would be better interface
    urls: Iterable[str],
    n_processes: int = 3,
    progress: bool = True,
    root_out_path: Path = Path("."),
) -> tuple[Path, ...]:
    # with open("thread-info.txt", "w") as fh:
    #     fh.write("Starting\n")

    thread_positions = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=n_processes) as pool:
        futures = [
            pool.submit(
                download_file_parallel_helper,
                url,
                # Obviously silly
                out_path=url.split("/")[-1],
                progress=progress,
                thread_positions=thread_positions,
            )
            for i, url in enumerate(urls)
        ]

        iterator_results = concurrent.futures.as_completed(futures)

        res_l = [future.result() for future in iterator_results]


to_download = [
    "https://esgf1.dkrz.de/thredds/fileServer/input4mips/input4MIPs/CMIP7/CMIP/CR/CR-CMIP-1-0-0/atmos/yr/co2/gm/v20250228/co2_input4MIPs_GHGConcentrations_CMIP_CR-CMIP-1-0-0_gm_0001-0999.nc",
    "https://esgf1.dkrz.de/thredds/fileServer/input4mips/input4MIPs/CMIP7/CMIP/CR/CR-CMIP-1-0-0/atmos/yr/co2/gm/v20250228/co2_input4MIPs_GHGConcentrations_CMIP_CR-CMIP-1-0-0_gm_1000-1749.nc",
    "https://esgf1.dkrz.de/thredds/fileServer/input4mips/input4MIPs/CMIP7/CMIP/CR/CR-CMIP-1-0-0/atmos/yr/co2/gm/v20250228/co2_input4MIPs_GHGConcentrations_CMIP_CR-CMIP-1-0-0_gm_1750-2022.nc",
    "https://esgf1.dkrz.de/thredds/fileServer/input4mips/input4MIPs/CMIP7/CMIP/CR/CR-CMIP-1-0-0/atmos/yr/co2/gm/v20250228/co2_input4MIPs_GHGConcentrations_CMIP_CR-CMIP-1-0-0_gm_0001-0999.nc",
    "https://esgf1.dkrz.de/thredds/fileServer/input4mips/input4MIPs/CMIP7/CMIP/CR/CR-CMIP-1-0-0/atmos/yr/co2/gm/v20250228/co2_input4MIPs_GHGConcentrations_CMIP_CR-CMIP-1-0-0_gm_1000-1749.nc",
    "https://esgf1.dkrz.de/thredds/fileServer/input4mips/input4MIPs/CMIP7/CMIP/CR/CR-CMIP-1-0-0/atmos/yr/co2/gm/v20250228/co2_input4MIPs_GHGConcentrations_CMIP_CR-CMIP-1-0-0_gm_1750-2022.nc",
    "https://esgf1.dkrz.de/thredds/fileServer/input4mips/input4MIPs/CMIP7/CMIP/CR/CR-CMIP-1-0-0/atmos/yr/co2/gm/v20250228/co2_input4MIPs_GHGConcentrations_CMIP_CR-CMIP-1-0-0_gm_0001-0999.nc",
    "https://esgf1.dkrz.de/thredds/fileServer/input4mips/input4MIPs/CMIP7/CMIP/CR/CR-CMIP-1-0-0/atmos/yr/co2/gm/v20250228/co2_input4MIPs_GHGConcentrations_CMIP_CR-CMIP-1-0-0_gm_1000-1749.nc",
    "https://esgf1.dkrz.de/thredds/fileServer/input4mips/input4MIPs/CMIP7/CMIP/CR/CR-CMIP-1-0-0/atmos/yr/co2/gm/v20250228/co2_input4MIPs_GHGConcentrations_CMIP_CR-CMIP-1-0-0_gm_1750-2022.nc",
]

download(
    to_download[::-1],
)
