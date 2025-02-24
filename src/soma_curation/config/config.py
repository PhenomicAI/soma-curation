from pydantic import BaseModel
from functools import lru_cache
from pydantic import computed_field

from ..mtx_collection.mtx_collection import MtxCollection
from ..schema import get_schema, DatabaseSchema
from ..atlas.crud import ExpandedPath


class PipelineConfig(BaseModel):
    """
    A pipeline configuration object.
    """

    atlas_name: str
    h5ad_storage_dir: str
    raw_storage_dir: str
    atlas_storage_dir: str
    processes: int = 4

    _db_schema: DatabaseSchema = None
    _mtx_collection: MtxCollection = None

    @property
    def db_schema(self) -> DatabaseSchema:
        if self._db_schema is None:
            self._db_schema = get_schema()
        return self._db_schema

    @property
    def mtx_collection(self) -> MtxCollection:
        if self._mtx_collection is None:
            self._mtx_collection = MtxCollection(
                storage_directory=self.raw_storage_dir, db_schema=self.db_schema, include=None
            )
        return self._mtx_collection


@lru_cache(maxsize=1)
def get_pipeline_config(
    atlas_name: str, raw_storage_dir: str, h5ad_storage_dir: str, atlas_storage_dir: str, processes: int = 4
) -> PipelineConfig:
    return PipelineConfig(
        atlas_name=atlas_name,
        raw_storage_dir=raw_storage_dir,
        h5ad_storage_dir=h5ad_storage_dir,
        atlas_storage_dir=atlas_storage_dir,
        processes=processes,
    )
