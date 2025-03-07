import yaml
import importlib

from cloudpathlib import AnyPath
from typing import Optional, Dict, Any

from .objects import DatabaseSchema, ValidationSchema, convert_types
from .defaults import DEFAULT_DATABASE_SCHEMA_DICT
from ..sc_logging import logger


def deep_merge_dict(orig: Dict[str, Any], override: Dict[str, Any]) -> Dict[str, Any]:
    """
    Recursively merges 'override' into 'orig'.
      - If both orig[k] and override[k] are dicts, merge them.
      - Otherwise override orig[k] with override[k].
    """
    for k, v in override.items():
        if k in orig and isinstance(orig[k], dict) and isinstance(v, dict):
            deep_merge_dict(orig[k], v)
        else:
            orig[k] = v
    return orig


def load_schema(db_config_uri: Optional[str] = None) -> DatabaseSchema:
    """
    Load a minimal or full DatabaseSchema YAML from a local file, S3, GCS, or Azure,
    merge it with your defaults, and return a DatabaseSchema object.

    If `db_config_uri` is None or empty, no user overrides are read, and the default
    database schema is used.

    If `validation_config_uri` is provided, merges that YAML with the default
    validation dictionary. Otherwise uses the default validation schema.
    """

    # 1. Load the user’s “database” config if db_config_uri is provided and non-empty
    if db_config_uri:
        db_path = AnyPath(db_config_uri)
        with db_path.open("r") as f:
            user_db_dict = yaml.safe_load(f)
        if not isinstance(user_db_dict, dict):
            raise ValueError(f"Invalid DB schema at {db_config_uri}")
    else:
        # No user DB config => no overrides
        user_db_dict = {}

    # 2. Merge user’s partial dictionary with your default DB schema
    merged_db_dict = deep_merge_dict(DEFAULT_DATABASE_SCHEMA_DICT.copy(), user_db_dict)

    # 3. Build ValidationSchema object
    merged_db_dict["VALIDATION_SCHEMA"] = ValidationSchema(**DEFAULT_DATABASE_SCHEMA_DICT["VALIDATION_SCHEMA"])

    # 4. Convert string types to pyarrow dtypes
    convert_types(merged_db_dict)

    # 5. Default core gene set path if not set
    if merged_db_dict["CORE_GENE_SET_PATH"] is None:
        logger.warning("CORE_GENE_SET_PATH is null or missing. Using dummy_core_geneset.tsv.gz instead.")
        core_gene_set_path = importlib.resources.files("soma_curation.constants").joinpath("dummy_core_geneset.tsv.gz")
        merged_db_dict["CORE_GENE_SET_PATH"] = str(core_gene_set_path)

    # 6. Instantiate your DatabaseSchema
    return DatabaseSchema(**merged_db_dict)
