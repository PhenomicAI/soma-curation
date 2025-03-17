import argparse
import sys
import pickle
import multiprocessing
from pathlib import Path

from ..sc_logging import logger, _set_level, init_worker_logging
from ..ingest.ingestion_funcs import (
    create_registration_mapping,
    ingest_h5ad_soma,
    resize_experiment,
    convert_and_std_mtx_to_h5ad,
)
from ..atlas.crud import AtlasManager
from ..executor.executors import MultiprocessingExecutor
from ..config.config import PipelineConfig

if multiprocessing.get_start_method(True) != "spawn":
    multiprocessing.set_start_method("spawn", True)


def main():
    parser = argparse.ArgumentParser(description="Run ingestion pipeline in parallel.")
    parser.add_argument(
        "--processes", type=int, default=4, help="Number of worker processes to use for parallel tasks."
    )
    parser.add_argument("--atlas-name", type=str, default="test", help="Name of the atlas.")
    parser.add_argument("--raw-storage-dir", type=str, default="human", help="Directory to raw storage.")
    parser.add_argument("--h5ad-storage-dir", type=str, default="h5ads", help="Directory to write H5AD files.")
    parser.add_argument("--atlas-storage-dir", type=str, default="./test", help="Directory to store the atlas.")
    parser.add_argument("--db-schema-fp", type=str, default=None, help="Filepath for database schema.")
    parser.add_argument("--log-dir", type=str, default="./logs", help="Directory to store logs.")
    args = parser.parse_args()

    # 1) Create a pipeline config
    pc = PipelineConfig(
        atlas_name=args.atlas_name,
        raw_storage_dir=args.raw_storage_dir,
        h5ad_storage_dir=args.h5ad_storage_dir,
        atlas_storage_dir=args.atlas_storage_dir,
        processes=args.processes,
        db_schema_uri=args.db_schema_fp,
        log_dir=args.log_dir,
    )

    # 2) Set up logging. We create a logs directory inside the atlas directory.
    log_dir = Path(pc.log_dir)
    log_dir.mkdir(parents=True, exist_ok=True)
    _set_level(level=10, add_file_handler=True, log_dir=log_dir, log_file=f"{pc.atlas_name}.log")

    logger.info("Starting pipeline execution...")
    logger.info(pc)

    # 5) Initialize AtlasManager
    am = AtlasManager(atlas_name=pc.atlas_name, storage_directory=pc.atlas_storage_dir, db_schema=pc.db_schema)

    if am.exists():
        raise ValueError(f"Atlas '{args.atlas_name}' already exists, aborting.")
    am.create()

    # ---------------------------------------------------------------------
    # PARALLEL STEP (1): Convert each study-sample into H5AD
    # ---------------------------------------------------------------------
    # Collect tasks: each item is (study, sample) so we can handle them
    # in parallel using the executor below.
    tasks_to_convert = []
    for study in pc.mtx_collection.list_studies():
        for sample in pc.mtx_collection.list_samples(study_name=study):
            tasks_to_convert.append((study, sample, pc))
            logger.info(f"Adding {(study, sample)} to the multiprocessing queue.")

    # Create a multiprocessing executor with the specified number of processes
    mp_executor = MultiprocessingExecutor(
        processes=args.processes,
        init_worker_logging=init_worker_logging,
        init_args=(10, log_dir.as_posix(), f"{pc.atlas_name}.log"),
    )
    logger.info("Starting parallel conversion to H5AD files...")
    convert_result = mp_executor.run(tasks_to_convert, convert_and_std_mtx_to_h5ad)

    # Gather successful conversions for the next step
    filenames = convert_result.successes
    logger.info(
        f"Parallel conversion complete. {convert_result.num_successes} successes, "
        f"{convert_result.num_failures} failures."
    )

    if not filenames:
        logger.error("No files were successfully converted; exiting early.")
        logger.error("Deleting atlas and exiting.")
        am.delete()
        sys.exit(1)

    # ---------------------------------------------------------------------
    # SERIAL STEP: Create registration mapping
    # ---------------------------------------------------------------------
    logger.info("Creating registration mapping (serial step)...")
    rm = create_registration_mapping(experiment_uri=str(am.experiment_path), filenames=filenames)
    with open("rm.pkl", "wb") as f:
        pickle.dump(rm, f)
    logger.info("Registration mapping created and saved to rm.pkl.")

    # ---------------------------------------------------------------------
    # SERIAL STEP: Resize the experiment
    # ---------------------------------------------------------------------
    logger.info("Resizing experiment (serial step)...")
    with open("rm.pkl", "rb") as f:
        rm = pickle.load(f)
    resize_experiment(str(am.experiment_path), registration_mapping=rm)
    logger.info("Experiment resized successfully.")

    # ---------------------------------------------------------------------
    # PARALLEL STEP (2): Ingest H5AD files into the SOMA experiment
    # ---------------------------------------------------------------------

    logger.info("Starting parallel ingestion of H5AD files into experiment...")
    tasks_for_ingestion = [(fname, str(am.experiment_path), rm) for fname in filenames]
    ingest_result = mp_executor.run(tasks_for_ingestion, ingest_h5ad_soma)
    logger.info(
        f"Ingestion complete. {ingest_result.num_successes} successes, " f"{ingest_result.num_failures} failures."
    )

    # If needed, you can decide how to handle the failures
    if ingest_result.num_failures > 0:
        logger.warning(f"{ingest_result.num_failures} H5AD ingestion tasks failed.")

    logger.info("All pipeline steps completed. Exiting.")
