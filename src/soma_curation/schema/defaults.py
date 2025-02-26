from typing import Any, Dict, List
import pyarrow as pa
import yaml
import pandas as pd
import importlib
from pydantic import BaseModel, Field
from functools import cached_property
from typing_extensions import Literal
from cloudpathlib import AnyPath
from pathlib import Path

# Default database schema
DEFAULT_DATABASE_SCHEMA_DICT = {
    "PAI_SCHEMA_VERSION": "1.0.0",
    "MEASUREMENT_RNA_NAME": "RNA",
    "PAI_PRESENCE_MATRIX_NAME": "feature_presence_matrix",
    "CORE_GENE_SET_PATH": None,
    "PAI_OBS_INDEX_COLUMNS": {"soma_joinid": "int64"},
    "PAI_OBS_SAMPLE_COLUMNS": {
        "sample_name": "categorical__large_string",
        "scrnaseq_protocol": "categorical__large_string",
        "study_name": "categorical__large_string",
    },
    "PAI_OBS_CELL_COLUMNS": {
        "barcode": "large_string",
        "cell_type": "categorical__large_string",
    },
    "PAI_OBS_COMPUTED_COLUMNS": {
        "nnz": "uint32",
        "umi_counts": "uint32",
        "pct_mito": "float32",
        "pct_ribo": "float32",
    },
    "COMPUTED_COLUMN_FUNCTIONS": {
        "nnz": "compute_nnz",
        "umi_counts": "compute_umi_counts",
        "pct_mito": "compute_pct_mito",
        "pct_ribo": "compute_pct_ribo",
    },
    "PAI_OBS_TERM_COLUMNS": {
        "soma_joinid": "int64",
        "barcode": "large_string",
        "cell_type": "categorical__large_string",
        "sample_name": "categorical__large_string",
        "scrnaseq_protocol": "categorical__large_string",
        "study_name": "categorical__large_string",
        "nnz": "uint32",
        "umi_counts": "uint32",
        "pct_mito": "float32",
        "pct_ribo": "float32",
    },
    "PAI_VAR_INDEX_COLUMNS": {"soma_joinid": "int64"},
    "PAI_VAR_COLUMNS": {
        "gene": "large_string",
        "ens": "large_string",
    },
    "PAI_VAR_TERM_COLUMNS": {
        "soma_joinid": "int64",
        "gene": "large_string",
        "ens": "large_string",
    },
    "PAI_X_LAYERS": {"col_raw": "uint32", "col_norm": "float32", "row_raw": "uint32", "row_norm": "float32"},
    "PAI_OBSM_INDEX_COLUMN": {"soma_joinid": "int64"},
    # Just an example – you can expand if you want
    "PAI_OBSM_LAYERS": {
        "embeddings": "float32",
        "umap": "float32",
    },
    "PAI_PRESENCE_LAYER": "uint8",
    # Platform configs for creation
    "PAI_OBS_PLATFORM_CONFIG": {
        "tiledb": {
            "create": {
                "capacity": 16384,
                "tile_order": "row-major",
                "cell_order": "row-major",
                "offsets_filters": ["DoubleDeltaFilter", {"_type": "ZstdFilter", "level": 9}],
                "allows_duplicates": False,
            }
        }
    },
    "PAI_VAR_PLATFORM_CONFIG": {
        "tiledb": {
            "create": {
                "capacity": 131072,
                "offsets_filters": ["DoubleDeltaFilter", {"_type": "ZstdFilter", "level": 9}],
                "allows_duplicates": False,
            }
        }
    },
    "PAI_X_LAYERS_PLATFORM_CONFIG": {
        "col_raw": {
            "tiledb": {
                "create": {
                    "capacity": 131072,
                    "dims": {
                        "soma_dim_0": {
                            "tile": 262144,
                            "filters": ["ByteShuffleFilter", {"_type": "ZstdFilter", "level": 9}],
                        },
                        "soma_dim_1": {
                            "tile": 1,
                            "filters": ["ByteShuffleFilter", {"_type": "ZstdFilter", "level": 9}],
                        },
                    },
                    "attrs": {"soma_data": {"filters": [{"_type": "ZstdFilter", "level": 5}]}},
                    "cell_order": "col-major",
                    "tile_order": "col-major",
                    "allows_duplicates": False,
                }
            }
        },
        "col_norm": {
            "tiledb": {
                "create": {
                    "capacity": 131072,
                    "dims": {
                        "soma_dim_0": {
                            "tile": 262144,
                            "filters": ["ByteShuffleFilter", {"_type": "ZstdFilter", "level": 9}],
                        },
                        "soma_dim_1": {
                            "tile": 1,
                            "filters": ["ByteShuffleFilter", {"_type": "ZstdFilter", "level": 9}],
                        },
                    },
                    "attrs": {"soma_data": {"filters": [{"_type": "ZstdFilter", "level": 5}]}},
                    "cell_order": "col-major",
                    "tile_order": "col-major",
                    "allows_duplicates": False,
                }
            }
        },
        "row_raw": {
            "tiledb": {
                "create": {
                    "capacity": 131072,
                    "dims": {
                        "soma_dim_0": {
                            "tile": 1,
                            "filters": ["ByteShuffleFilter", {"_type": "ZstdFilter", "level": 9}],
                        },
                        "soma_dim_1": {
                            "tile": 35804,
                            "filters": ["ByteShuffleFilter", {"_type": "ZstdFilter", "level": 9}],
                        },
                    },
                    "attrs": {"soma_data": {"filters": [{"_type": "ZstdFilter", "level": 5}]}},
                    "cell_order": "row-major",
                    "tile_order": "row-major",
                    "allows_duplicates": False,
                }
            }
        },
        "row_norm": {
            "tiledb": {
                "create": {
                    "capacity": 131072,
                    "dims": {
                        "soma_dim_0": {
                            "tile": 1,
                            "filters": ["ByteShuffleFilter", {"_type": "ZstdFilter", "level": 9}],
                        },
                        "soma_dim_1": {
                            "tile": 35804,
                            "filters": ["ByteShuffleFilter", {"_type": "ZstdFilter", "level": 9}],
                        },
                    },
                    "attrs": {"soma_data": {"filters": [{"_type": "ZstdFilter", "level": 5}]}},
                    "cell_order": "row-major",
                    "tile_order": "row-major",
                    "allows_duplicates": False,
                }
            }
        },
    },
    "PAI_PRESENCE_PLATFORM_CONFIG": {
        "tiledb": {
            "create": {
                "capacity": 131072,
                "dims": {
                    "soma_dim_0": {
                        "tile": 1,
                        "filters": ["ByteShuffleFilter", {"_type": "ZstdFilter", "level": 9}],
                    },
                    "soma_dim_1": {
                        "tile": 35804,
                        "filters": ["ByteShuffleFilter", {"_type": "ZstdFilter", "level": 9}],
                    },
                },
                "cell_order": "row-major",
                "tile_order": "row-major",
                "allows_duplicates": False,
            }
        }
    },
}

# Default validation schema
DEFAULT_VALIDATION_DICT = {
    "REQUIRED_OBS_COLUMNS": ["barcode", "sample_name", "study_name"],
    "REQUIRED_VAR_COLUMNS": ["gene"],
    "GENE_INTERSECTION_THRESHOLD_FRAC": 0.1,
}
