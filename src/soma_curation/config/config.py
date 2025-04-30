import tiledbsoma as soma

from pydantic import BaseModel, ConfigDict, computed_field
from functools import lru_cache
from typing import Optional, Union
from enum import Enum

from ..collection import MtxCollection, H5adCollection
from ..schema import DatabaseSchema, load_schema


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


class RawCollectionType(str, Enum):
    MTX = "mtx"
    H5AD = "h5ad"


class PipelineConfig(BaseModel):
    """
    A pipeline configuration object.
    """

    atlas_name: str
    h5ad_storage_dir: str
    raw_storage_dir: str
    atlas_storage_dir: str
    log_dir: str
    raw_collection_type: RawCollectionType = RawCollectionType.MTX

    filenames_pickle: Optional[str] = "filenames.pkl"
    registration_mapping_pickle: Optional[str] = "rm.pkl"

    db_schema_uri: Optional[str] = None
    include_studies: Optional[list[str]] = None

    processes: int = 4

    model_config = ConfigDict(arbitrary_types_allowed=True)

    @computed_field(repr=False)
    @property
    def db_schema(self) -> DatabaseSchema:
        return load_schema(self.db_schema_uri)

    @computed_field
    @property
    def collection(self) -> Union[MtxCollection, H5adCollection]:
        if self.raw_collection_type == RawCollectionType.MTX:
            return MtxCollection(
                storage_directory=self.raw_storage_dir,
                db_schema=self.db_schema,
                include=self.include_studies
            )
        else:
            return H5adCollection(
                storage_directory=self.raw_storage_dir,
                include=self.include_studies
            )


# TODO: figure out if LRU cache works in multiprocessing
@lru_cache(maxsize=1)
def get_pipeline_config(
    atlas_name: str,
    raw_storage_dir: str,
    h5ad_storage_dir: str,
    atlas_storage_dir: str,
    processes: int = 4,
    db_schema_uri: Optional[str] = None,
    raw_collection_type: RawCollectionType = RawCollectionType.MTX,
) -> PipelineConfig:
    return PipelineConfig(
        atlas_name=atlas_name,
        raw_storage_dir=raw_storage_dir,
        h5ad_storage_dir=h5ad_storage_dir,
        atlas_storage_dir=atlas_storage_dir,
        processes=processes,
        db_schema_uri=db_schema_uri,
        raw_collection_type=raw_collection_type,
    )
