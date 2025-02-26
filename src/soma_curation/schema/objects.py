import pyarrow as pa
import pandas as pd
import importlib

from pydantic import BaseModel, ConfigDict, Field, computed_field
from functools import cached_property
from typing import Dict, Any, List, Literal, Set


class ValidationSchema(BaseModel):
    REQUIRED_OBS_COLUMNS: List[str]
    REQUIRED_VAR_COLUMNS: List[str]
    GENE_INTERSECTION_THRESHOLD_FRAC: float


class DatabaseSchema(BaseModel):
    """
    DatabaseSchema class defines the schema used for structuring and validating single-cell RNA sequencing data in TileDB-SOMA format.

    Attributes:
    - PAI_SCHEMA_VERSION: str
        The version of the schema.
    - MEASUREMENT_RNA_NAME: str
        The name of the RNA measurement.
    - PAI_PRESENCE_MATRIX_NAME: str
        The name of the presence matrix.
    - PAI_OBS_INDEX_COLUMNS: Dict[str, pa.DataType]
        Dictionary of observation index columns and their data types.
    - PAI_OBS_SAMPLE_COLUMNS: Dict[str, pa.DataType]
        Dictionary of observation sample columns and their data types.
    - PAI_OBS_CELL_COLUMNS: Dict[str, pa.DataType]
        Dictionary of observation cell columns and their data types.
    - PAI_OBS_COMPUTED_COLUMNS: Dict[str, pa.DataType]
        Dictionary of computed columns and their data types.
    - PAI_OBS_TERM_COLUMNS: Dict[str, pa.DataType]
        Dictionary of observation term columns and their data types.
    - PAI_VAR_INDEX_COLUMNS: Dict[str, pa.DataType]
        Dictionary of variable index columns and their data types.
    - PAI_VAR_COLUMNS: Dict[str, pa.DataType]
        Dictionary of variable columns and their data types.
    - PAI_VAR_TERM_COLUMNS: Dict[str, pa.DataType]
        Dictionary of variable term columns and their data types.
    - PAI_X_LAYERS: Dict[str, pa.DataType]
        Dictionary of X layers and their data types.
    - PAI_OBSM_INDEX_COLUMN: Dict[str, pa.DataType]
        Dictionary of observation matrix index columns and their data types.
    - PAI_PRESENCE_LAYER: pa.DataType
        Data type of the presence layer.
    - PAI_OBS_PLATFORM_CONFIG: Dict[str, Dict[str, Dict[str, Any]]]
        Platform configuration for observation columns.
    - PAI_VAR_PLATFORM_CONFIG: Dict[str, Dict[str, Dict[str, Any]]]
        Platform configuration for variable columns.
    - PAI_X_LAYERS_PLATFORM_CONFIG: Dict[str, Dict[str, Dict[str, Dict[str, Any]]]]
        Platform configuration for X layers.
    - PAI_PRESENCE_PLATFORM_CONFIG: Dict[str, Dict[str, Any]]
        Platform configuration for presence layers.
    - COMPUTED_COLUMN_FUNCTIONS: Dict[str, Any]
        Dictionary of computed column functions.
    - VALIDATION_SCHEMA: ValidationSchema
        Validation schema object.

    Methods:
    - fetch_and_functions: Fetch and set functions for computed columns.
    - apply_filters: Apply filters to observation and variable columns based on their data types.

    Properties:
    - CORE_GENES: set
        Set of core genes loaded from the core gene set file.
    """

    # Names of things
    PAI_SCHEMA_VERSION: str = Field(repr=True)
    MEASUREMENT_RNA_NAME: str = Field(repr=False)
    PAI_PRESENCE_MATRIX_NAME: str = Field(repr=False)

    # Core gene path
    CORE_GENE_SET_PATH: str = Field(repr=False)

    # Columns and types
    PAI_OBS_INDEX_COLUMNS: Dict[str, pa.DataType] = Field(repr=False)
    PAI_OBS_SAMPLE_COLUMNS: Dict[str, pa.DataType] = Field(repr=False)
    PAI_OBS_CELL_COLUMNS: Dict[str, pa.DataType] = Field(repr=False)
    PAI_OBS_COMPUTED_COLUMNS: Dict[str, pa.DataType] = Field(repr=False)
    PAI_OBS_TERM_COLUMNS: Dict[str, pa.DataType] = Field(repr=False)
    PAI_VAR_INDEX_COLUMNS: Dict[str, pa.DataType] = Field(repr=False)
    PAI_VAR_COLUMNS: Dict[str, pa.DataType] = Field(repr=False)
    PAI_VAR_TERM_COLUMNS: Dict[str, pa.DataType] = Field(repr=False)
    PAI_X_LAYERS: Dict[str, pa.DataType] = Field(repr=False)
    PAI_OBSM_INDEX_COLUMN: Dict[str, pa.DataType] = Field(repr=False)
    PAI_PRESENCE_LAYER: pa.DataType = Field(repr=False)

    # Platform configs for creation
    PAI_OBS_PLATFORM_CONFIG: Dict[str, Dict[str, Dict[str, Any]]] = Field(repr=False)
    PAI_VAR_PLATFORM_CONFIG: Dict[str, Dict[str, Dict[str, Any]]] = Field(repr=False)
    PAI_X_LAYERS_PLATFORM_CONFIG: Dict[str, Dict[str, Dict[str, Dict[str, Any]]]] = Field(repr=False)
    PAI_PRESENCE_PLATFORM_CONFIG: Dict[str, Dict[str, Any]] = Field(repr=False)

    COMPUTED_COLUMN_FUNCTIONS: Dict[str, Any] = Field(default={}, repr=False)

    VALIDATION_SCHEMA: ValidationSchema = Field(repr=False)

    model_config = ConfigDict(arbitrary_types_allowed=True)

    def model_post_init(self, ctx):
        """
        Post-initialization hook to apply filters and fetch functions for computed columns.

        Parameters:
        - ctx: Context
            The context in which the model is initialized.
        """
        [self.apply_filters(l) for l in ["obs", "var"]]
        self.fetch_and_functions()

    def fetch_and_functions(self):
        """
        Fetch and set functions for computed columns from the specified module.
        """
        module = importlib.import_module("soma_curation.dataset.standardize", package=__package__)

        for col in self.PAI_OBS_COMPUTED_COLUMNS:
            func = getattr(module, self.COMPUTED_COLUMN_FUNCTIONS[col])

            self.COMPUTED_COLUMN_FUNCTIONS[col] = func

    def apply_filters(self, layer: Literal["obs", "var"]):
        """
        Apply filters to observation and variable columns based on their data types.

        Parameters:
        - layer: Literal["obs", "var"]
            The layer to which filters are applied.

        Raises:
        - ValueError
            If the dimensional attributes are not integers.
        """
        layer_mapping = {
            "obs": {
                "index": self.PAI_OBS_INDEX_COLUMNS,
                "platform_config": self.PAI_OBS_PLATFORM_CONFIG,
                "term_columns": self.PAI_OBS_TERM_COLUMNS,
            },
            "var": {
                "index": self.PAI_OBS_INDEX_COLUMNS,
                "platform_config": self.PAI_VAR_PLATFORM_CONFIG,
                "term_columns": self.PAI_VAR_TERM_COLUMNS,
            },
        }

        layer_index: Dict[str, pa.DataType] = layer_mapping[layer]["index"]
        layer_term_columns: Dict[str, pa.DataType] = layer_mapping[layer]["term_columns"]
        layer_platform_config: Dict[str, Dict[str, Dict[str, Any]]] = layer_mapping[layer]["platform_config"]

        for column, dtype in layer_term_columns.items():
            if column in layer_index:
                dim_attr_key = "dims"
                if dim_attr_key not in layer_platform_config["tiledb"]["create"]:
                    layer_platform_config["tiledb"]["create"][dim_attr_key] = {}

                if pa.types.is_integer(dtype):
                    layer_platform_config["tiledb"]["create"][dim_attr_key][column] = {
                        "filters": ["DoubleDeltaFilter", {"_type": "ZstdFilter", "level": 9}]
                    }
                else:
                    raise ValueError("Set dimensional attributes to integers. Other data types are not efficient")
            else:
                dim_attr_key = "attrs"
                if dim_attr_key not in layer_platform_config["tiledb"]["create"]:
                    layer_platform_config["tiledb"]["create"][dim_attr_key] = {}

                # Numeric
                if pa.types.is_integer(dtype) or pa.types.is_floating(dtype):
                    layer_platform_config["tiledb"]["create"][dim_attr_key][column] = {
                        "filters": ["ByteShuffleFilter", {"_type": "ZstdFilter", "level": 9}]
                    }
                # Dictionary
                elif pa.types.is_dictionary(dtype):
                    layer_platform_config["tiledb"]["create"][dim_attr_key][column] = {
                        "filters": [{"_type": "ZstdFilter", "level": 9}]
                    }
                else:
                    layer_platform_config["tiledb"]["create"][dim_attr_key][column] = {
                        "filters": [{"_type": "ZstdFilter", "level": 9}]
                    }

    @computed_field(repr=False)
    @cached_property
    def CORE_GENES(self) -> Set[str]:
        """
        Load and return the core genes from the core gene set file.

        Returns:
        - set
            Set of core genes.
        """
        return set(pd.read_csv(self.CORE_GENE_SET_PATH, sep="\t", header=None)[0])

    @computed_field(repr=False)
    @cached_property
    def SORTED_CORE_GENES(self) -> List[str]:
        """
        Return the sorted core genes.

        Returns:
        - List[str]
            Set of sorted core genes.
        """
        return sorted(self.CORE_GENES)

    @computed_field(repr=False)
    @cached_property
    def NUM_GENES(self) -> int:
        """
        Load and return the num genes.

        Returns:
        - int
            Number of genes.
        """
        return len(self.SORTED_CORE_GENES)

    @computed_field(repr=False)
    @cached_property
    def VAR_DF(self) -> pd.DataFrame:
        """
        Load and return the var df based on the core genes.

        Returns:
        - pd.DataFrame
            Dataframe of of core genes, soma_joinid and ensembl.
        """
        var_df = pd.DataFrame(
            {
                "gene": self.SORTED_CORE_GENES,
                "soma_joinid": list(range(self.NUM_GENES)),
                "ens": self.SORTED_CORE_GENES,
            }
        )

        return var_df


def convert_types(d: dict):
    """
    Recursively convert string data types (e.g. 'int64', 'large_string')
    to PyArrow data types in-place.
    """
    pyarrow_mapping = {
        "large_string": pa.large_string(),
        "categorical__large_string": pa.dictionary(pa.int32(), pa.large_string(), ordered=False),
        "string": pa.string(),
        "int64": pa.int64(),
        "int32": pa.int32(),
        "uint8": pa.uint8(),
        "uint32": pa.uint32(),
        "float32": pa.float32(),
        "bool_": pa.bool_(),
    }
    for key, val in d.items():
        if isinstance(val, dict):
            convert_types(val)  # recurse
        elif isinstance(val, str):
            if val in pyarrow_mapping:
                d[key] = pyarrow_mapping[val]
