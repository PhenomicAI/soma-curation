import pyarrow as pa
import yaml
import tiledbsoma as soma
import tiledb
import pandas as pd
import importlib

from functools import cached_property
from pydantic import BaseModel, ConfigDict, validate_call, Field, computed_field
from typing import Dict, Any, Literal, List, Set


DEFAULT_TILEDB_CONFIG = {
    "py.max_incomplete_retries": 100,
    "py.init_buffer_bytes": 536870912,
    "soma.init_buffer_bytes": 536870912,
    "sm.consolidation.buffer_size": 1073741824,
    "sm.consolidation.step_max_frags": 100,
    "sm.consolidation.step_min_frags": 3,
}


def SOMA_TileDB_Context() -> soma.options.SOMATileDBContext:
    """
    Create and return a SOMA TileDB context with default configuration.

    Returns:
    - soma.options.SOMATileDBContext
        The configured SOMA TileDB context.
    """
    return soma.options.SOMATileDBContext(tiledb_ctx=TileDB_Ctx(), timestamp=None)


def TileDB_Ctx() -> tiledb.Ctx:
    """
    Create and return a TileDB context with default configuration.

    Returns:
    - tiledb.Ctx
        The configured TileDB context.
    """
    return tiledb.Ctx(DEFAULT_TILEDB_CONFIG)


# your YAML configuration
sample_global_config = """
# Names of arrays in the PAI-SOMA schema
PAI_SCHEMA_VERSION: 1.0.0
MEASUREMENT_RNA_NAME: RNA
PAI_PRESENCE_MATRIX_NAME: feature_presence_matrix

# Core genes path
CORE_GENE_SET_PATH: 

# Columns and data types
PAI_OBS_INDEX_COLUMNS: &PAI_OBS_INDEX_COLUMNS
    soma_joinid: int64

PAI_OBS_SAMPLE_COLUMNS: &PAI_OBS_SAMPLE_COLUMNS
    sample_name: categorical__large_string
    scrnaseq_protocol: categorical__large_string
    study_name: categorical__large_string

PAI_OBS_CELL_COLUMNS: &PAI_OBS_CELL_COLUMNS
    barcode: large_string
    cell_type: categorical__large_string

PAI_OBS_COMPUTED_COLUMNS: &PAI_OBS_COMPUTED_COLUMNS
    nnz: uint32
    umi_counts: uint32
    pct_mito: float32
    pct_ribo: float32

COMPUTED_COLUMN_FUNCTIONS:
    nnz: compute_nnz
    umi_counts: compute_umi_counts
    pct_mito: compute_pct_mito
    pct_ribo: compute_pct_ribo
    
PAI_OBS_TERM_COLUMNS:
    <<: *PAI_OBS_INDEX_COLUMNS
    <<: *PAI_OBS_CELL_COLUMNS
    <<: *PAI_OBS_SAMPLE_COLUMNS
    <<: *PAI_OBS_COMPUTED_COLUMNS

PAI_VAR_INDEX_COLUMNS: &PAI_VAR_INDEX_COLUMNS
    soma_joinid: int64

PAI_VAR_COLUMNS: &PAI_VAR_COLUMNS
    gene: large_string
    ens: large_string
    
PAI_VAR_TERM_COLUMNS:
    <<: *PAI_VAR_INDEX_COLUMNS
    <<: *PAI_VAR_COLUMNS

PAI_X_LAYERS:
    col_raw: uint32
    col_norm: float32
    row_raw: uint32
    row_norm: float32

PAI_OBSM_INDEX_COLUMN:
    soma_joinid: int64

PAI_OBSM_LAYERS:
    embeddings: float32    
    umap: float32

PAI_PRESENCE_LAYER: uint8

# Platform config for TileDB arrays
PAI_OBS_PLATFORM_CONFIG:
    tiledb:
        create:
            capacity: 16384
            tile_order: row-major
            cell_order: row-major
            offsets_filters: [DoubleDeltaFilter, {_type: ZstdFilter, level: 9}]
            allows_duplicates: False

PAI_VAR_PLATFORM_CONFIG:
    tiledb:
        create:
            capacity: 131072
            offsets_filters: [DoubleDeltaFilter, {_type: ZstdFilter, level: 9}]
            allows_duplicates: False

PAI_X_LAYERS_PLATFORM_CONFIG:
    col_raw:
        tiledb:
            create:
                capacity: 131072
                dims:
                    soma_dim_0: 
                        tile: 262144
                        filters: [ByteShuffleFilter, {_type: ZstdFilter, level: 9}]
                    soma_dim_1:
                        tile: 1
                        filters: [ByteShuffleFilter, {_type: ZstdFilter, level: 9}]
                attrs:
                    soma_data:
                        filters: [{_type: ZstdFilter, level: 5}]
                cell_order: col-major
                tile_order: col-major
                allows_duplicates: False
    col_norm:
        tiledb:
            create: 
                capacity: 131072
                dims:
                    soma_dim_0:
                        tile: 262144
                        filters: [ByteShuffleFilter, {_type: ZstdFilter, level: 9}]
                    soma_dim_1:
                        tile: 1
                        filters: [ByteShuffleFilter, {_type: ZstdFilter, level: 9}]
                attrs:
                    soma_data:
                        filters: [{_type: ZstdFilter, level: 5}]
                cell_order: col-major
                tile_order: col-major
                allows_duplicates: False
    row_raw:
        tiledb:
            create:
                capacity: 131072
                dims: 
                    soma_dim_0:
                        tile: 1
                        filters: [ByteShuffleFilter, {_type: ZstdFilter, level: 9}]
                    soma_dim_1:
                        tile: 35804
                        filters: [ByteShuffleFilter, {_type: ZstdFilter, level: 9}]
                attrs:
                    soma_data:
                        filters: [{_type: ZstdFilter, level: 5}]
                cell_order: row-major
                tile_order: row-major
                allows_duplicates: False
    row_norm:
        tiledb:
            create:
                capacity: 131072
                dims:
                    soma_dim_0:
                        tile: 1
                        filters: [ByteShuffleFilter, {_type: ZstdFilter, level: 9}]
                    soma_dim_1:
                        tile: 35804
                        filters: [ByteShuffleFilter, {_type: ZstdFilter, level: 9}]
                attrs:
                    soma_data:
                        filters: [{_type: ZstdFilter, level: 5}]
                cell_order: row-major
                tile_order: row-major
                allows_duplicates: False

PAI_PRESENCE_PLATFORM_CONFIG:
    tiledb:
        create:
            capacity: 131072
            dims:
                soma_dim_0:
                    tile: 1
                    filters: [ByteShuffleFilter, {_type: ZstdFilter, level: 9}]
                soma_dim_1:
                    tile: 35804
                    filters: [ByteShuffleFilter, {_type: ZstdFilter, level: 9}]
            cell_order: row-major
            tile_order: row-major
            allows_duplicates: False
"""

sample_validation_schema = """
REQUIRED_OBS_COLUMNS: 
    - barcode 
    - sample_name 
    - study_name 
REQUIRED_VAR_COLUMNS: 
    - gene
GENE_INTERSECTION_THRESHOLD_FRAC: 0.5
"""


class ValidationSchema(BaseModel):
    """
    ValidationSchema class defines the schema used for validating single-cell RNA sequencing data.

    Attributes:
    - REQUIRED_OBS_COLUMNS: List[str]
        List of required columns in the .obs attribute.
    - REQUIRED_VAR_COLUMNS: List[str]
        List of required columns in the .var attribute.
    - GENE_INTERSECTION_THRESHOLD_FRAC: float
        The threshold fraction for gene intersection.
    """

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
        for col in self.PAI_OBS_COMPUTED_COLUMNS:
            module = importlib.import_module("data_curation.soma.dataset.standardize")
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


def convert_types(dict_: str) -> Dict[str, Any]:
    """
    Convert string data types to PyArrow data types in a dictionary.

    Parameters:
    - dict_: str
        The dictionary containing data type strings.

    Returns:
    - Dict[str, Any]
        The dictionary with converted PyArrow data types.
    """
    pyarrow_mapping = {
        "large_string": pa.large_string(),
        "categorical__large_string": pa.dictionary(pa.int32(), pa.large_string(), ordered=0),
        "string": pa.string(),
        "int64": pa.int64(),
        "int32": pa.int32(),
        "uint8": pa.uint8(),
        "uint32": pa.uint32(),
        "float32": pa.float32(),
        "bool_": pa.bool_(),
    }

    for k, v in dict_.items():
        if isinstance(v, dict):
            convert_types(dict_[k])
        elif isinstance(v, str):
            if v in pyarrow_mapping.keys():
                dict_[k] = pyarrow_mapping[v]

    return None


@validate_call
def read_yaml(yaml_string: str) -> Dict[str, Any]:
    """
    Read and parse a YAML string into a dictionary.

    Parameters:
    - yaml_string: str
        The YAML string to parse.

    Returns:
    - Dict[str, Any]
        The parsed dictionary.
    """
    dict_: Dict[str, Any] = yaml.safe_load(yaml_string)
    return dict_


@validate_call
def read_globals_yaml(yaml_string: str) -> Dict[str, Any]:
    """
    Read and parse a YAML string into a dictionary and convert data types.

    Parameters:
    - yaml_string: str
        The YAML string to parse.

    Returns:
    - Dict[str, Any]
        The parsed and converted dictionary.
    """
    dict_ = read_yaml(yaml_string)
    convert_types(dict_)
    return dict_


@validate_call
def get_validation_schema() -> ValidationSchema:
    """
    Get the validation schema.

    Returns:
    - ValidationSchema
        The validation schema object.
    """
    return ValidationSchema(**sample_validation_schema)


@validate_call
def get_schema() -> DatabaseSchema:
    """
    Get the database schema.

    Returns:
    - DatabaseSchema
        The database schema object.
    """
    dict_ = read_globals_yaml(sample_global_config)
    validation_schema = get_validation_schema()
    dict_["VALIDATION_SCHEMA"] = validation_schema
    return DatabaseSchema(**dict_)
