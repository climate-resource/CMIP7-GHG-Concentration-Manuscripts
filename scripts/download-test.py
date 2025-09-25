from local.data_loading import fetch_and_load_ghg_file
from local.esgf.search.search_query import KnownIndexNode

fetch_and_load_ghg_file(
    ghg="co2",
    grid="gm",
    time_sampling="yr",
    cmip_era="CMIP7",
    source_id="CR-CMIP-1-0-0",
    # cmip_era="CMIP6Plus",
    # source_id="CR-CMIP-0-4-0",
    # index_node=KnownIndexNode.ORNL,
    # index_node=KnownIndexNode.DKRZ,
    # index_node=KnownIndexNode.CEDA,
    index_node=KnownIndexNode.NCI,
)
