"""
Query ESGF
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, TypeVar

import httpx
from httpx_retries import Retry, RetryTransport
from loguru import logger

from local.esgf.esgf_dataset_collection import ESGFDatasetCollection
from local.esgf.models import ESGFDataset, ESGFFileAccessURL

T = TypeVar("T")


@dataclass(frozen=True)
class SearchRetryConfig:
    """
    Retry configuration for ESGF search requests

    By default, this retries the transient exceptions handled by
    [httpx-retries][], including [httpx.RemoteProtocolError][].
    """

    total: int = 5
    """
    Maximum number of retries after the first request
    """

    backoff_factor: float = 0.5
    """
    Exponential backoff factor between retries
    """

    max_backoff_wait: float = 30.0
    """
    Maximum time to sleep between retry attempts
    """


DEFAULT_SEARCH_RETRY_CONFIG = SearchRetryConfig()


# TODO: think about how to organise this mapping etc.


def query_esgf(
    endpoint: str,
    query_terms: dict[str, tuple[str, ...]],
    distrib: bool = True,
    limit: int = 1000,
    retry_config: SearchRetryConfig | None = DEFAULT_SEARCH_RETRY_CONFIG,
) -> ESGFDatasetCollection:
    """
    Query ESGF and parse the file results into datasets
    """
    raw_response = query_esgf_files(
        endpoint=endpoint,
        query_terms=query_terms,
        distrib=distrib,
        limit=limit,
        retry_config=retry_config,
    )

    esgf_dataset_collection = parse_raw_esgf_search_result(
        raw_response.raise_for_status().json()
    )

    return esgf_dataset_collection


def query_esgf_files(
    endpoint: str,
    query_terms: dict[str, tuple[str, ...]],
    distrib: bool = True,
    limit: int = 1000,
    retry_config: SearchRetryConfig | None = DEFAULT_SEARCH_RETRY_CONFIG,
    # TODO: support other configuration options?
    # Much longer list here:
    # https://esgf.github.io/esg-search/ESGF_Search_RESTful_API.html#the-esgf-search-restful-api
) -> httpx.Response:
    """
    Query ESGF for file records
    """
    # Query files
    # Then process them later into results
    # that are grouped by dataset etc.
    format: str = "application/solr+json"
    result_type: str = "File"

    params = {
        **query_terms,
        "format": format,
        "distrib": distrib,
        "limit": limit,
        "type": result_type,
    }
    logger.debug(f"Querying {endpoint} with {params=}")
    if retry_config is None:
        client = httpx.Client()
    else:
        retry = Retry(
            total=retry_config.total,
            backoff_factor=retry_config.backoff_factor,
            max_backoff_wait=retry_config.max_backoff_wait,
        )
        transport = RetryTransport(retry=retry)
        client = httpx.Client(transport=transport)

    try:
        response = client.get(endpoint, params=params)
    except httpx.HTTPError as exc:
        msg = f"Error raised while trying to access {endpoint} with {params=}"
        raise AssertionError(msg) from exc
    finally:
        client.close()

    logger.debug(f"Query URL: {response.url}")
    try:
        response.raise_for_status()
    except httpx.HTTPStatusError as exc:
        msg = f"Error raised while trying to access {response.url}"
        raise AssertionError(msg) from exc

    return response


def get_single_value(ind: dict[Any, Any], key: Any) -> T:
    """
    Get a single value from a potentially list-valued ESGF field
    """
    res = ind[key]
    if isinstance(res, list):
        if len(res) != 1:
            msg = f"More than one value for {key=}, values={res}"
            raise AssertionError(msg)

        return res[0]

    return res


def parse_raw_esgf_search_result(
    raw_search_json: dict[str, Any],
) -> ESGFDatasetCollection:
    """
    Parse a raw ESGF search response into a dataset collection
    """
    # work out which dataset each file belongs to
    # parse everything into dataset/ESGFDataset objects

    num_results = raw_search_json["response"]["numFound"]
    max_supported_results_without_scrolling = 10_000
    if num_results > max_supported_results_without_scrolling:
        raise NotImplementedError

    # This is by far the hardest part.
    # We go over all the retrieved file records.
    # We figure out the dataset they belong to.
    # We also need to associate the file access options with each file.
    # Once we have done this, we can go back over the grouped files
    # to create our dataset instances.
    dataset_file_ids = {}
    file_access_urls = {}
    results_by_file_id = {}
    for result_d in raw_search_json["response"]["docs"]:
        # Get rid of the node before creating the dataset ID
        our_dataset_id = result_d["dataset_id"].split("|")[0]
        # Get rid of the node before creating the file ID
        our_file_id = result_d["id"].split("|")[0]

        if our_file_id not in results_by_file_id:
            results_by_file_id[our_file_id] = []

        results_by_file_id[our_file_id].append(result_d)

        for access_url in result_d["url"]:
            # https://esgf.github.io/esg-search/ESGF_Search_RESTful_API.html#access-urls
            url, mime_type, service_name = access_url.split("|")
            esgf_file_access_url = ESGFFileAccessURL(
                url=url,
                mime_type=mime_type,
                service_name=service_name,
            )
            if our_file_id not in file_access_urls:
                file_access_urls[our_file_id] = []

            file_access_urls[our_file_id].append(esgf_file_access_url)

        if our_dataset_id not in dataset_file_ids:
            dataset_file_ids[our_dataset_id] = []

        if our_file_id not in dataset_file_ids[our_dataset_id]:
            dataset_file_ids[our_dataset_id].append(our_file_id)

    esgf_datasets_l = []
    for dataset_id, file_ids in dataset_file_ids.items():
        esgf_file_records = [
            esgf_result
            for id in file_ids
            # Each file can be listed more than one place due to duplication over nodes
            for esgf_result in results_by_file_id[id]
        ]
        esgf_dataset = ESGFDataset.from_esgf_file_records(esgf_file_records)

        esgf_datasets_l.append(esgf_dataset)

    res = ESGFDatasetCollection(esgf_datasets=tuple(esgf_datasets_l))

    return res
