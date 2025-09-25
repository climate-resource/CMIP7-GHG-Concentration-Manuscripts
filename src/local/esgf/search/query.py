"""
Query ESGF
"""

from __future__ import annotations

from typing import Any

import httpx
from loguru import logger


def query_esgf(
    endpoint: str,
    query_terms: dict[str, tuple[str, ...]],
    format: str = "application/solr+json",
    type: str = "File",
    distrib: bool = True,
    limit: int = 1000,
) -> dict[str, Any]:
    params = {
        **query_terms,
        "format": format,
        "distrib": distrib,
        "limit": limit,
        "type": type,
    }
    logger.debug(f"Querying {endpoint} with {params=}")
    response = httpx.get(endpoint, params=params)
    logger.debug(f"Query URL: {response.url}")

    res = response.raise_for_status().json()

    num_results = res["response"]["numFound"]
    # TODO: put this check somewhere else
    max_supported_results_without_scrolling = 10_000
    if num_results > max_supported_results_without_scrolling:
        raise NotImplementedError

    logger.debug(f"{num_results} {'result' if num_results == 1 else 'results'}")

    return res
