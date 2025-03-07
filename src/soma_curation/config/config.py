import tiledbsoma as soma

from pydantic import BaseModel, ConfigDict, computed_field
from functools import lru_cache
from typing import Optional

from ..mtx_collection.mtx_collection import MtxCollection
from ..schema.objects import DatabaseSchema
from ..schema.load import load_schema


DEFAULT_TILEDB_CONFIG = {
    "py.max_incomplete_retries": 100,
    "py.init_buffer_bytes": 536870912,
    "soma.init_buffer_bytes": 536870912,
    "sm.consolidation.buffer_size": 1073741824,
    "sm.consolidation.step_max_frags": 100,
    "sm.consolidation.step_min_frags": 3,
}


# TODO: convert to singleton to speed up (if necessary)
def SOMA_TileDB_Context() -> soma.options.SOMATileDBContext:
    """
    Create and return a SOMA TileDB context with default configuration.

    Returns:
    - soma.options.SOMATileDBContext
        The configured SOMA TileDB context.
    """
    return soma.options.SOMATileDBContext(tiledb_config=DEFAULT_TILEDB_CONFIG, timestamp=None)


class PipelineConfig(BaseModel):
    """
    A pipeline configuration object.
    """

    atlas_name: str
    h5ad_storage_dir: str
    raw_storage_dir: str
    atlas_storage_dir: str

    db_schema_uri: Optional[str] = None

    processes: int = 4

    model_config = ConfigDict(arbitrary_types_allowed=True)

    @computed_field(repr=False)
    @property
    def db_schema(self) -> DatabaseSchema:
        if not self.db_schema_uri:
            # Optionally: return the default schema if no URI provided
            return load_schema(None)
        return load_schema(self.db_schema_uri)

    @computed_field
    @property
    def mtx_collection(self) -> MtxCollection:
        return MtxCollection(storage_directory=self.raw_storage_dir, db_schema=self.db_schema)


# TODO: figure out if LRU cache works in multiprocessing
@lru_cache(maxsize=1)
def get_pipeline_config(
    atlas_name: str,
    raw_storage_dir: str,
    h5ad_storage_dir: str,
    atlas_storage_dir: str,
    processes: int = 4,
    db_schema_uri: Optional[str] = None,
) -> PipelineConfig:
    return PipelineConfig(
        atlas_name=atlas_name,
        raw_storage_dir=raw_storage_dir,
        h5ad_storage_dir=h5ad_storage_dir,
        atlas_storage_dir=atlas_storage_dir,
        processes=processes,
        db_schema_uri=db_schema_uri,
    )
