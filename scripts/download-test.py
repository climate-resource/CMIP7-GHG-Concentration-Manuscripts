"""
Test of downloading and saving results
"""

from pathlib import Path

from local.data_loading import fetch_and_load_ghg_dataset
from local.esgf.db_helpers import create_all_tables, get_sqlite_engine
from local.esgf.search.search_query import KnownIndexNode

REPO_ROOT = Path(__file__).parents[1]
local_data_root_dir = REPO_ROOT / "data" / "raw" / "esgf"
local_data_root_dir.mkdir(exist_ok=True, parents=True)
sqlite_file = REPO_ROOT / "download-test-database.db"
# Obviously we wouldn't delete the database every time
# in production, but while experimenting it's handy
# to always start with a clean slate.
if sqlite_file.exists():
    sqlite_file.unlink()

engine = get_sqlite_engine(sqlite_file)
create_all_tables(engine)


ds = fetch_and_load_ghg_dataset(
    local_data_root_dir=local_data_root_dir,
    ghg="hfc23",
    grid="gnz",
    time_sampling="mon",
    cmip_era="CMIP7",
    source_id="CR-CMIP-1-0-0",
    # cmip_era="CMIP6Plus",
    # source_id="CR-CMIP-0-4-0",
    index_node=KnownIndexNode.DKRZ,
    # index_node=KnownIndexNode.ORNL,
    engine=engine,
)
print(ds)
