import yaml
import importlib

from cloudpathlib import AnyPath
from typing import Optional

from .objects import DatabaseSchema, ValidationSchema, convert_types
from .defaults import DEFAULT_DATABASE_SCHEMA_DICT, DEFAULT_VALIDATION_DICT
from ..sc_logging import logger


def merge_dicts(user: dict, default: dict) -> dict:
    """
    Recursively merge keys from 'default' into 'user' if they don't exist in 'user'.
    If both 'user' and 'default' have a dict under the same key,
    merge them recursively.
    """
    for k, v in default.items():
        if k not in user:
            user[k] = v
        else:
            if isinstance(v, dict) and isinstance(user[k], dict):
                merge_dicts(user[k], v)
    return user


def load_schema(
    db_config_uri: Optional[str] = None,
    validation_config_uri: Optional[str] = None,
) -> DatabaseSchema:
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
    merged_db_dict = merge_dicts(user_db_dict, DEFAULT_DATABASE_SCHEMA_DICT.copy())

    # 3. Load user’s partial validation config if provided
    if validation_config_uri:
        val_path = AnyPath(validation_config_uri)
        with val_path.open("r") as vf:
            user_val_dict = yaml.safe_load(vf)
        if not isinstance(user_val_dict, dict):
            raise ValueError(f"Invalid validation schema at {validation_config_uri}")

        merged_val_dict = merge_dicts(user_val_dict, DEFAULT_VALIDATION_DICT.copy())
    else:
        # Use default validation if not provided
        merged_val_dict = DEFAULT_VALIDATION_DICT.copy()

    # 4. Build ValidationSchema object
    validation_obj = ValidationSchema(**merged_val_dict)

    # 5. Insert that into final dictionary
    merged_db_dict["VALIDATION_SCHEMA"] = validation_obj

    # 6. Convert string types to pyarrow dtypes
    convert_types(merged_db_dict)

    # 7. Default core gene set path if not set
    if merged_db_dict["CORE_GENE_SET_PATH"] is None:
        logger.warning("CORE_GENE_SET_PATH is null or missing. Using dummy_core_geneset.tsv.gz instead.")
        core_gene_set_path = importlib.resources.files("soma_curation.constants").joinpath("dummy_core_geneset.tsv.gz")
        merged_db_dict["CORE_GENE_SET_PATH"] = str(core_gene_set_path)

    # 8. Instantiate your DatabaseSchema
    return DatabaseSchema(**merged_db_dict)
