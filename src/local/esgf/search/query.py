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
    distrib: bool = True,
    limit: int = 1000,
) -> dict[str, Any]:
    params = {
        **query_terms,
        "format": format,
        "distrib": distrib,
        "limit": limit,
    }
    logger.debug(f"Querying {endpoint} with {params=}")
    response = httpx.get(endpoint, params=params)

    res = response.raise_for_status().json()

    num_results = res["response"]["numFound"]
    # TODO: put this somewhere else
    max_supported_results_without_scrolling = 10_000
    if num_results > max_supported_results_without_scrolling:
        raise NotImplementedError

    logger.debug(f"{num_results} {'result' if num_results == 1 else 'results'}")
    # Actual entries to parse are in res["response"]["docs"]
