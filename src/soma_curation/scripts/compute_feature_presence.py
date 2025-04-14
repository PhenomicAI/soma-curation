# TODO: Investigate dense arrays for presence matrix # atlas/crud.py
import argparse
import tiledbsoma as soma
import multiprocessing
import pandas as pd

from typing import Tuple

from cloudpathlib import AnyPath

from ..collection import MtxCollection, H5adCollection
from ..schema import load_schema, DatabaseSchema
from ..executor.executors import MultiprocessingExecutor
from ..sc_logging import logger, _set_level, init_worker_logging
from ..ingest.ingestion_funcs import compute_presence_matrix

if multiprocessing.get_start_method(True) != "spawn":
    multiprocessing.set_start_method("spawn", True)


def determine_sample_df_to_process(experiment_uri: str, schema: DatabaseSchema) -> Tuple[pd.DataFrame, int]:
    """
    Determine which samples need to be processed for feature presence matrix computation and the length of the presence matrix to resize to.

    Args:
        experiment_uri (str): URI of the experiment
        schema: Database schema containing PAI_PRESENCE_MATRIX_NAME

    Returns:
        pd.DataFrame: DataFrame containing sample_name, study_name, and sample_idx for samples that need processing
        int: The length of the presence matrix to resize to
    """
    exp = soma.Experiment.open(experiment_uri)
    sample_names_df = (
        exp.obs.read(column_names=["sample_name", "study_name"])
        .concat()
        .to_pandas()
        .drop_duplicates()
        .reset_index(drop=True)
    )
    sample_names_df["sample_idx"] = sample_names_df["sample_name"].cat.codes
    resize_length = sample_names_df["sample_idx"].max() + 1
    logger.info(
        f"Found {len(sample_names_df)} samples in SOMA experiment, going to resize presence matrix to {resize_length} rows"
    )
    last_sample_idx = exp.ms["RNA"][schema.PAI_PRESENCE_MATRIX_NAME].non_empty_domain()[0][1]
    exp.close()

    if last_sample_idx == 0:
        logger.info("No samples exist in the presence matrix, retaining all samples")
        sample_names_df = sample_names_df
    else:
        logger.info(f"Considering all sample idxs from {last_sample_idx + 1} to {resize_length - 1}")
        sample_names_df = sample_names_df[~sample_names_df["sample_idx"].isin(range(last_sample_idx + 1))]
    return sample_names_df, resize_length


def main():
    _set_level(level=10, add_file_handler=True, log_dir=AnyPath("logs"), log_file="pipeline.log")
    parser = argparse.ArgumentParser(description="Run feature presence computation.")
    parser.add_argument("--n-processes", type=int, default=4, help="URI of the atlas.")
    parser.add_argument("--exp-uri", type=str, default="test", help="URI of the atlas.")
    parser.add_argument("--raw-storage-dir", type=str, default="human", help="Directory to raw storage.")
    parser.add_argument("--db-schema-fp", type=str, default=None, help="Filepath for database schema.")
    parser.add_argument(
        "--raw-collection-type",
        type=str,
        choices=["mtx", "h5ad"],
        default="mtx",
        help="Type of collection to process (mtx or h5ad)",
    )
    args = parser.parse_args()

    schema = load_schema(args.db_schema_fp)

    if args.raw_collection_type == "mtx":
        collection = MtxCollection(storage_directory=args.raw_storage_dir, db_schema=schema)
    else:
        collection = H5adCollection(storage_directory=args.raw_storage_dir)

    global_var_list = schema.SORTED_CORE_GENES

    logger.info("Determining samples to process...")
    samples_to_generate_df, resize_length = determine_sample_df_to_process(args.exp_uri, schema)

    # Determine what to resize the presence matrix to
    presence_matrix_shape = (resize_length, len(global_var_list))
    logger.info(f"Resizing presence matrix to {presence_matrix_shape}")
    with soma.Experiment.open(args.exp_uri, mode="w") as exp:
        exp.ms["RNA"][schema.PAI_PRESENCE_MATRIX_NAME].resize(presence_matrix_shape)
    logger.info("Resized presence matrix")

    tasks_for_ingestion = [
        (sidx, sample_name, study_name, collection, global_var_list, schema.PAI_PRESENCE_MATRIX_NAME, args.exp_uri)
        for sidx, sample_name, study_name in zip(
            samples_to_generate_df["sample_idx"],
            samples_to_generate_df["sample_name"],
            samples_to_generate_df["study_name"],
        )
    ]
    logger.info(f"Running {len(tasks_for_ingestion)} tasks in parallel")
    mp_executor = MultiprocessingExecutor(
        processes=args.n_processes,
        init_worker_logging=init_worker_logging,
        init_args=(10, AnyPath("logs").as_posix(), "pipeline.log"),
    )
    ingest_result = mp_executor.run(tasks_for_ingestion, compute_presence_matrix)
    logger.info(
        f"Ingestion complete. {ingest_result.num_successes} successes, " f"{ingest_result.num_failures} failures."
    )

    # If needed, you can decide how to handle the failures
    if ingest_result.num_failures > 0:
        logger.warning(f"{ingest_result.num_failures} H5AD ingestion tasks failed.")


if __name__ == "__main__":
    main()
