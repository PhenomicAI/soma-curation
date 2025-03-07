import pyarrow as pa
import shutil
import tiledbsoma as soma

from datetime import datetime, timezone
from pydantic import BaseModel, Field, computed_field, ConfigDict, AfterValidator
from pathlib import Path
from typing import Union, Generator
from contextlib import contextmanager
from typing_extensions import Annotated

from ..schema.objects import DatabaseSchema
from ..utils.git_utils import get_git_commit_sha
from ..sc_logging import logger
from ..config.config import SOMA_TileDB_Context


def expand_paths(value: Union[str, Path]) -> Path:
    if isinstance(value, str):
        value = Path(value)
    return value.expanduser()


ExpandedPath = Annotated[str, AfterValidator(expand_paths)]


class AtlasManager(BaseModel):
    """
    AtlasManager class is designed to manage the creation, deletion, and schema management for single-cell RNA sequencing atlases.

    Attributes:
    - atlas_name: str
        The name of the atlas.
    - globals_: DatabaseSchema
        The schema defining the structure and validation rules for the data.
    - storage_directory: Path
        The directory where the atlas is stored.
    - context: soma.SOMATileDBContext
        The context for TileDB-SOMA operations, with a default factory function.

    Methods:
    - exists: Checks if the atlas already exists.
    - create: Creates a new atlas.
    - delete: Deletes an existing atlas.
    """

    atlas_name: str
    storage_directory: ExpandedPath
    globals_: DatabaseSchema = Field(repr=False)
    context: soma.SOMATileDBContext = Field(default_factory=SOMA_TileDB_Context, repr=False)

    model_config = ConfigDict(arbitrary_types_allowed=True)

    @computed_field(repr=False)
    @property
    def experiment_path(self) -> Path:
        """
        Compute the path to the experiment directory.

        Returns:
        - Path
            The path to the experiment directory.
        """
        return self.storage_directory / self.atlas_name

    @computed_field(repr=True)
    @property
    def version(self) -> Path:
        """
        Compute the path to the experiment directory.

        Returns:
        - Path
            The path to the experiment directory.
        """
        if self.exists():
            with self.open(mode="r") as exp:
                if "pai_soma_object_version" in exp.metadata:
                    return exp.metadata["pai_soma_object_version"]
                else:
                    return "v1.0"
        else:
            return "v1.0"

    @contextmanager
    def open(self, **kwargs) -> Generator[soma.Experiment, None, None]:
        yield soma.Experiment.open(self.experiment_path.as_posix(), context=self.context, **kwargs)

    def exists(self) -> bool:
        """
        Check if the atlas already exists.

        Returns:
        - bool
            True if the atlas exists, False otherwise.
        """
        if self.experiment_path.exists():
            try:
                self.open(mode="r")
                return True
            except soma._exception.DoesNotExistError:
                return False
            except Exception as e:
                logger.critical("Unexpected error occured", e)
                return False
        return False

    def create(self) -> None:
        """
        Create a new atlas.

        This method initializes the experiment directory, sets up metadata, and creates necessary collections and dataframes.
        """
        if self.exists():
            logger.warning(
                f"Atlas {self.atlas_name} in directory {self.storage_directory} exists, skipping creation..."
            )
            return None
        else:
            if not self.storage_directory.expanduser().exists():
                self.storage_directory.expanduser().mkdir()
        logger.info(f"Creating atlas {self.atlas_name} in directory {self.storage_directory}...")

        with soma.Experiment.create(self.experiment_path.as_posix(), context=self.context) as experiment:
            experiment.metadata["created_on"] = datetime.now(tz=timezone.utc).isoformat(timespec="seconds")
            experiment.metadata["pai_schema_version"] = self.globals_.PAI_SCHEMA_VERSION

            experiment.metadata["pai_soma_object_version"] = self.version

            sha = get_git_commit_sha()
            experiment.metadata["git_commit_sha"] = sha
            experiment.metadata["atlas_name"] = self.atlas_name

            # create `obs`
            obs_schema = pa.schema(
                self.globals_.PAI_OBS_TERM_COLUMNS,
                metadata={
                    k: "nullable"
                    for k in (
                        [x[0] for x in self.globals_.PAI_OBS_SAMPLE_COLUMNS]
                        + [x[0] for x in self.globals_.PAI_OBS_CELL_COLUMNS]
                        + [x[0] for x in self.globals_.PAI_OBS_COMPUTED_COLUMNS]
                    )
                },
            )
            experiment.add_new_dataframe(
                "obs",
                schema=obs_schema,
                index_column_names=[x[0] for x in self.globals_.PAI_OBS_INDEX_COLUMNS],
                platform_config=self.globals_.PAI_OBS_PLATFORM_CONFIG,
            )

            # create `ms`
            measurements = experiment.add_new_collection("ms")

            # Create RNA measurements
            rna_measurement = measurements.add_new_collection(self.globals_.MEASUREMENT_RNA_NAME, soma.Measurement)

            # create empty `obsm` collection in the measurement
            rna_measurement.add_new_collection("obsm")

            # create `var` in the measurement
            var_schema = pa.schema(self.globals_.PAI_VAR_TERM_COLUMNS)
            var = rna_measurement.add_new_dataframe(
                "var",
                schema=var_schema,
                index_column_names=[x[0] for x in self.globals_.PAI_VAR_INDEX_COLUMNS],
                platform_config=self.globals_.PAI_VAR_PLATFORM_CONFIG,
                domain=[[0, len(self.globals_.VAR_DF) - 1]],
            )

            # TODO: clean this little issue here
            assert set([x[0] for x in self.globals_.PAI_VAR_COLUMNS]) == set(
                ["gene", "ens"]
            ), "Make sure that your var df aligns"
            table = pa.Table.from_pandas(self.globals_.VAR_DF, preserve_index=False)
            var.write(table)

            # create `X` in the measurement
            X_collection = rna_measurement.add_new_collection("X")
            for layer_name, dtype in self.globals_.PAI_X_LAYERS:
                platform_config = self.globals_.PAI_X_LAYERS_PLATFORM_CONFIG[layer_name]
                if layer_name.startswith("row"):
                    logger.info(f"Converting row-reads tile width to {self.globals_.NUM_GENES}")
                    platform_config["tiledb"]["create"]["dims"]["soma_dim_1"]["tile"] = self.globals_.NUM_GENES

                X_collection.add_new_sparse_ndarray(
                    layer_name,
                    type=dtype,
                    shape=(None, self.globals_.NUM_GENES),
                    platform_config=platform_config,
                )

            logger.info(f"Converting presence tile width to {self.globals_.NUM_GENES}")
            platform_config = self.globals_.PAI_PRESENCE_PLATFORM_CONFIG
            platform_config["tiledb"]["create"]["dims"]["soma_dim_1"]["tile"] = self.globals_.NUM_GENES
            rna_measurement.add_new_sparse_ndarray(
                self.globals_.PAI_PRESENCE_MATRIX_NAME,
                type=self.globals_.PAI_PRESENCE_LAYER,
                shape=(None, self.globals_.NUM_GENES),
                platform_config=platform_config,
            )

    def delete(self) -> None:
        """
        Delete an existing atlas.

        This method removes the directory and all its contents.
        """
        if self.exists():
            logger.info(f"Deleting atlas {self.atlas_name} in directory {self.storage_directory}...")
            shutil.rmtree(self.experiment_path)
        else:
            logger.info(
                f"Atlas {self.atlas_name} does not exist in directory {self.storage_directory}, skipping deletion..."
            )
