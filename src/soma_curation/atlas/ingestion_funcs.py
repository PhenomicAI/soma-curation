import tiledbsoma.io

from typing import List, Literal

from ..sc_logging import logger
from ..schema import SOMA_TileDB_Context


def _create_registration_mapping(
    experiment_uri: str,
    filenames: List[str],
    measurement_name: Literal["RNA"] = "RNA",
    obs_field_name: str = "barcode",
    var_field_name: str = "gene",
    context: tiledbsoma.SOMATileDBContext = SOMA_TileDB_Context(),
) -> tiledbsoma.io.ExperimentAmbientLabelMapping:

    rm = tiledbsoma.io.register_h5ads(
        experiment_uri=experiment_uri,
        h5ad_file_names=filenames,
        measurement_name=measurement_name,
        obs_field_name=obs_field_name,
        var_field_name=var_field_name,
        context=context,
    )

    return rm


def _ingest_h5ad_worker(
    experiment_uri: str,
    h5ad_path: str,
    registration_mapping: tiledbsoma.io.ExperimentAmbientLabelMapping,
    measurement_name: Literal["RNA"] = "RNA",
    obs_field_name: str = "barcode",
    var_field_name: str = "gene",
    x_layer_name: str = "row_raw",
    raw_x_layer_name: str = "row_raw",
    context: tiledbsoma.SOMATileDBContext = SOMA_TileDB_Context(),
) -> bool:
    try:
        logger.info(f"Worker ingesting file: {h5ad_path}")

        tiledbsoma.io.resize_experiment(
            uri=experiment_uri, nobs=registration_mapping.get_obs_shape(), nvars=registration_mapping.get_var_shapes()
        )

        tiledbsoma.io.from_h5ad(
            experiment_uri=experiment_uri,
            input_path=h5ad_path,
            measurement_name=measurement_name,
            registration_mapping=registration_mapping,
            obs_id_name=obs_field_name,
            var_id_name=var_field_name,
            X_layer_name=x_layer_name,
            raw_X_layer_name=raw_x_layer_name,
            context=context,
        )
        logger.info(f"Finished ingesting file: {h5ad_path}")

        return True
    except Exception as ex:
        logger.error(f"Failed to ingest {h5ad_path}, error: {ex}")
        raise Exception(f"{ex}")
