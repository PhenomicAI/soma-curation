import argparse
import sys
import pickle
import multiprocessing
from pathlib import Path
from cloudpathlib import AnyPath

from ..sc_logging import logger, _set_level, init_worker_logging
from ..ingest.ingestion_funcs import (
    create_registration_mapping,
    ingest_h5ad_soma,
    resize_experiment,
    convert_and_std_mtx_to_h5ad,
    convert_and_std_h5ad_to_h5ad,
)
from ..atlas.crud import AtlasManager
from ..executor.executors import MultiprocessingExecutor
from ..config.config import PipelineConfig, RawCollectionType

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
    parser.add_argument(
        "--filenames-pkl",
        type=str,
        default=None,
        help="Pickle path for storing/loading successful H5AD filenames. "
        "If provided and the file exists, Step 1 is skipped.",
    )
    parser.add_argument(
        "--registration-mapping-pkl",
        type=str,
        default=None,
        help="Pickle path for storing/loading registration mapping. "
        "If provided and the file exists, creating the registration mapping is skipped.",
    )
    parser.add_argument(
        "--raw-collection-type",
        type=str,
        choices=["mtx", "h5ad"],
        default="mtx",
        help="Type of collection to process (mtx or h5ad)",
    )
    parser.add_argument(
        "--include-studies",
        type=str,
        nargs="+",
        default=None,
        help="List of studies to include. If not provided, all studies will be processed.",
    )
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
        filenames_pickle=args.filenames_pkl or "filenames.pkl",
        registration_mapping_pickle=args.registration_mapping_pkl or "rm.pkl",
        raw_collection_type=RawCollectionType(args.raw_collection_type),
        include_studies=args.include_studies,
    )

    # 2) Set up logging. We create a logs directory inside the atlas directory.
    log_dir = Path(pc.log_dir)
    log_dir.mkdir(parents=True, exist_ok=True)
    _set_level(level=10, add_file_handler=True, log_dir=log_dir, log_file=f"{pc.atlas_name}.log")

    logger.info("Starting pipeline execution...")
    logger.info(pc)

    # 5) Initialize AtlasManager
    am = AtlasManager(atlas_name=pc.atlas_name, storage_directory=pc.atlas_storage_dir, db_schema=pc.db_schema)

    if not am.exists():
        am.create()

    # ---------------------------------------------------------------------
    # STEP (1): Convert each study-sample into H5AD
    # ---------------------------------------------------------------------
    filenames = []
    filenames_pkl = AnyPath(pc.filenames_pickle)

    if filenames_pkl.is_file():
        logger.info(f"Found existing filenames pickle at {str(filenames_pkl)}. Skipping H5AD conversion.")
        with open(filenames_pkl, "rb") as f:
            filenames = pickle.load(f)
    else:
        tasks_to_convert = []

        if pc.raw_collection_type == RawCollectionType.MTX:
            for study in pc.collection.list_studies():
                for sample in pc.collection.list_samples(study_name=study):
                    tasks_to_convert.append((study, sample, pc))
                    logger.info(f"Adding MTX {(study, sample)} to the multiprocessing queue.")
            conversion_func = convert_and_std_mtx_to_h5ad
        else:
            for h5ad_file in pc.collection.list_h5ad_files():
                tasks_to_convert.append((h5ad_file, pc))
                logger.info(f"Adding H5AD file {h5ad_file} to the multiprocessing queue.")
            conversion_func = convert_and_std_h5ad_to_h5ad

        # Create a multiprocessing executor with the specified number of processes
        mp_executor = MultiprocessingExecutor(
            processes=pc.processes,
            init_worker_logging=init_worker_logging,
            init_args=(10, log_dir.as_posix(), f"{pc.atlas_name}.log"),
        )
        logger.info("Starting parallel conversion to standardized H5AD files...")
        convert_result = mp_executor.run(tasks_to_convert, conversion_func)

        # Gather successful conversions for the next step
        filenames = convert_result.successes
        logger.info(
            f"Parallel conversion complete. {convert_result.num_successes} successes, "
            f"{convert_result.num_failures} failures."
        )

        if not filenames:
            logger.error("No files were successfully converted; exiting early.")
            sys.exit(1)

        # Save the successfully converted filenames to pickle
        with filenames_pkl.open("wb") as f:
            pickle.dump(filenames, f)
        logger.info(f"Filenames saved to {str(filenames_pkl)}.")

    # ---------------------------------------------------------------------
    # STEP (2): Create registration mapping or load from file
    # ---------------------------------------------------------------------
    registration_mapping_pkl = AnyPath(pc.registration_mapping_pickle)
    if registration_mapping_pkl.is_file():
        logger.info(f"Found existing registration mapping at {str(registration_mapping_pkl)}, skipping creation.")
        with registration_mapping_pkl.open("rb") as f:
            rm = pickle.load(f)
    else:
        logger.info("Creating registration mapping (serial step)...")
        rm = create_registration_mapping(experiment_uri=str(am.experiment_path), filenames=filenames)
        with registration_mapping_pkl.open("wb") as f:
            pickle.dump(rm, f)
        logger.info(f"Registration mapping created and saved to {str(registration_mapping_pkl)}.")

    # ---------------------------------------------------------------------
    # STEP (3): Resize the experiment
    # ---------------------------------------------------------------------
    logger.info("Resizing experiment (serial step)...")
    rm.prepare_experiment(str(am.experiment_path))
    logger.info("Experiment resized successfully.")

    # ---------------------------------------------------------------------
    # STEP (4): Ingest H5AD files into the SOMA experiment
    # ---------------------------------------------------------------------

    logger.info("Starting parallel ingestion of H5AD files into experiment...")
    tasks_for_ingestion = [(fname, str(am.experiment_path), rm) for fname in filenames]
    mp_executor = MultiprocessingExecutor(
        processes=pc.processes,
        init_worker_logging=init_worker_logging,
        init_args=(10, log_dir.as_posix(), f"{pc.atlas_name}.log"),
    )
    ingest_result = mp_executor.run(tasks_for_ingestion, ingest_h5ad_soma)
    logger.info(
        f"Ingestion complete. {ingest_result.num_successes} successes, " f"{ingest_result.num_failures} failures."
    )

    # If needed, you can decide how to handle the failures
    if ingest_result.num_failures > 0:
        logger.warning(f"{ingest_result.num_failures} H5AD ingestion tasks failed.")

    logger.info("All pipeline steps completed. Exiting.")
