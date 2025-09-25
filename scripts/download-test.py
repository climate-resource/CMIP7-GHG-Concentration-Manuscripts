from local.data_loading import fetch_and_load_ghg_file

fetch_and_load_ghg_file(
    ghg="co2",
    grid="gm",
    time_sampling="yr",
    cmip_era="CMIP7",
    source_id="CR-CMIP-1-0-0",
)
