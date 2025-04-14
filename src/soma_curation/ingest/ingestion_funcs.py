import tiledbsoma.io
import tiledbsoma as soma

from typing import List, Literal, Union
from cloudpathlib import AnyPath
from pathlib import Path
from tiledbsoma.io import ExperimentAmbientLabelMapping
import pyarrow as pa

from ..sc_logging import logger
from ..dataset.anndataset import AnnDataset
from ..config.config import PipelineConfig, SOMA_TileDB_Context
from ..collection import MtxCollection, H5adCollection


def create_registration_mapping(
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


def resize_experiment(
    experiment_uri: str,
    registration_mapping: tiledbsoma.io.ExperimentAmbientLabelMapping,
    context: tiledbsoma.SOMATileDBContext = SOMA_TileDB_Context(),
) -> tiledbsoma.io.ExperimentAmbientLabelMapping:
    try:
        tiledbsoma.io.resize_experiment(
            uri=experiment_uri,
            nobs=registration_mapping.get_obs_shape(),
            nvars=registration_mapping.get_var_shapes(),
            context=context,
        )

        return True
    except Exception:
        logger.warning("Resize experiment failed")
        raise ValueError("Resize experiment failed")


def ingest_h5ad_soma(path: str, experiment_path: str, rm: ExperimentAmbientLabelMapping):
    """
    Ingest a single H5AD file into the existing experiment using
    the precomputed registration mapping.
    """
    try:
        logger.info(f"[ingest_h5ad] Ingesting {path}")
        _ingest_h5ad_soma(
            experiment_uri=experiment_path,
            h5ad_path=path,
            registration_mapping=rm,
        )
        logger.info(f"Successfully ingested {path}")
        return path  # Return the filename on success
    except Exception as e:
        logger.error(f"Failed to ingest {path}: {e}")
        raise


def convert_and_std_mtx_to_h5ad(study_name: str, sample_name: str, pc: PipelineConfig):
    """
    Function to convert a single (study, sample) to an H5AD file.
    Returns the path to the resulting H5AD file on success.
    Raises an exception if anything fails.
    """
    try:
        logger.info(f"[convert_to_h5ad] Processing study='{study_name}', sample='{sample_name}'")
        adata = pc.collection.get_anndata(study_name=study_name, sample_name=sample_name)
    except Exception as e:
        logger.error(
            f"Error fetching and assembling study='{study_name}', sample='{sample_name}' from raw storage: {e}"
        )
        raise
    try:
        anndataset = AnnDataset(artifact=adata, db_schema=pc.db_schema)
    except Exception as e:
        logger.error(f"Error converting study='{study_name}', sample='{sample_name}' to AnnDataset object: {e}")
        raise

    try:
        anndataset.standardize()
    except Exception as e:
        logger.error(f"Error standardizing study='{study_name}', sample='{sample_name}' to H5AD: {e}")
        raise

    try:
        # Construct the output filename
        filename = AnyPath(pc.h5ad_storage_dir) / f"{study_name}-{sample_name}.h5ad"
        anndataset.write(filename)
        logger.info(f"Successfully converted '{study_name}'-'{sample_name}' -> '{filename}'")
        return filename.as_posix() if isinstance(filename, Path) else filename.as_uri()
    except Exception as e:
        logger.error(f"Error writing study={study_name}, sample={sample_name} to H5AD: {e}")
        raise

# TODO: This is a duplicate of the function above. Might need to refactor.
def convert_and_std_h5ad_to_h5ad(h5ad_path: str, pc: PipelineConfig):
    """
    Function to convert a single H5AD file to a standardized H5AD file.
    Returns the path to the resulting standardized H5AD file on success.
    """
    try:
        logger.info(f"[convert_to_h5ad] Processing H5AD file: '{h5ad_path}'")
        adata = pc.collection.get_anndata(h5ad_path)
    except Exception as e:
        logger.error(f"Error fetching H5AD file '{h5ad_path}' from storage: {e}")
        raise

    try:
        anndataset = AnnDataset(artifact=adata, db_schema=pc.db_schema)
    except Exception as e:
        logger.error(f"Error converting '{h5ad_path}' to AnnDataset object: {e}")
        raise

    try:
        anndataset.standardize()
    except Exception as e:
        logger.error(f"Error standardizing '{h5ad_path}' to H5AD: {e}")
        raise

    try:
        # Construct the output filename - keep the original name but in the new directory
        original_name = Path(h5ad_path).name
        filename = AnyPath(pc.h5ad_storage_dir) / original_name
        anndataset.write(filename)
        logger.info(f"Successfully standardized '{h5ad_path}' -> '{filename}'")
        return filename.as_posix() if isinstance(filename, Path) else filename.as_uri()
    except Exception as e:
        logger.error(f"Error writing standardized H5AD for '{h5ad_path}': {e}")
        raise


def _ingest_h5ad_soma(
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

def compute_presence_matrix(
    sample_idx: int,
    sample_name: str,
    study_name: str,
    collection: Union[MtxCollection, H5adCollection],
    global_var_list: List[str],
    pai_presence_matrix_name: str,
    exp_uri: str,
):
    try:
        logger.info(f"Computing presence matrix for {sample_name}-{study_name}")
        coo_pm_matrix = collection.presence_matrix(
            study_name=study_name, sample_name=sample_name, global_var_list=global_var_list
        )
        pa_table = pa.Table.from_pydict(
            {
                "soma_dim_0": coo_pm_matrix.row + sample_idx,
                "soma_dim_1": coo_pm_matrix.col,
                "soma_data": coo_pm_matrix.data,
            }
        )
        with soma.Experiment.open(exp_uri, mode="w") as exp:
            exp.ms["RNA"][pai_presence_matrix_name].write(pa_table)
        logger.info(f"Successfully computed presence matrix for {sample_name} {study_name}")
        return sample_idx
    except Exception as e:
        print(e)
        logger.info(f"Error computing presence matrix for {sample_name} {study_name}: {e}")
        return None