import argparse
import sys
import pickle
from pathlib import Path


from ..schema import get_schema
from ..sc_logging import logger, _set_level
from ..constants import create_dummy_structure
from ..ingest.ingestion_funcs import (
    _create_registration_mapping,
    _ingest_h5ad_worker,
    _resize_experiment,
)
from ..atlas.crud import AtlasManager
from ..dataset.anndataset import AnnDataset
from ..mtx_collection.mtx_collection import MtxCollection
from ..executor.executors import MultiprocessingExecutor


def main():
    parser = argparse.ArgumentParser(description="Run ingestion pipeline in parallel.")
    parser.add_argument(
        "--processes", type=int, default=4, help="Number of worker processes to use for parallel tasks."
    )
    parser.add_argument("--atlas-name", type=str, default="test", help="Name of the atlas.")
    parser.add_argument("--raw-storage-dir", type=str, default="human", help="Path to raw storage directory.")
    parser.add_argument("--h5ad-storage-dir", type=str, default="h5ads", help="Directory to write H5AD files.")
    parser.add_argument("--atlas-storage-dir", type=str, default="./test", help="Local directory to store the atlas.")
    parser.add_argument("--organism", type=str, default="human", help="Organism name (if needed).")
    args = parser.parse_args()

    # 1) Create a dummy structure for the raw data (if needed).
    create_dummy_structure(args.raw_storage_dir)

    # 2) Set up logging. We create a logs directory inside the atlas directory.
    log_dir = Path(args.atlas_storage_dir) / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    _set_level(level=20, add_file_handler=True, log_dir=log_dir, log_file=f"{args.atlas_name}.log")

    logger.info("Starting pipeline execution...")

    # 3) Get the database schema
    db_schema = get_schema()

    # 4) Initialize MtxCollection
    mtx_collection = MtxCollection(storage_directory=args.raw_storage_dir, db_schema=db_schema, include=None)

    # 5) Initialize AtlasManager
    am = AtlasManager(atlas_name=args.atlas_name, storage_directory=args.atlas_storage_dir, globals_=db_schema)

    if am.exists():
        raise ValueError(f"Atlas '{args.atlas_name}' already exists, aborting.")
    am.create()

    # ---------------------------------------------------------------------
    # PARALLEL STEP (1): Convert each study-sample into H5AD
    # ---------------------------------------------------------------------
    # Collect tasks: each item is (study, sample) so we can handle them
    # in parallel using the executor below.
    tasks_to_convert = []
    for study in mtx_collection.list_studies():
        for sample in mtx_collection.list_samples(study_name=study):
            tasks_to_convert.append((study, sample))

    def convert_to_h5ad(task):
        """
        Function to convert a single (study, sample) to an H5AD file.
        Returns the path to the resulting H5AD file on success.
        Raises an exception if anything fails.
        """
        study, sample = task
        try:
            logger.info(f"[convert_to_h5ad] Processing study={study}, sample={sample}")
            adata = mtx_collection.get_anndata(study_name=study, sample_name=sample)
            anndataset = AnnDataset(artifact=adata, database_schema=db_schema)
            anndataset.standardize()
            # Construct the output filename
            filename = Path(args.h5ad_storage_dir) / f"{study}-{sample}.h5ad"
            anndataset.write(filename)
            logger.info(f"Successfully converted {study}-{sample} -> {filename}")
            return filename.as_posix()
        except Exception as e:
            logger.error(f"Error converting study={study}, sample={sample} to H5AD: {e}")
            # Raise so that it goes to the executor's failure list
            raise

    # Create a multiprocessing executor with the specified number of processes
    mp_executor = MultiprocessingExecutor(processes=args.processes)
    logger.info("Starting parallel conversion to H5AD files...")
    convert_result = mp_executor.run(tasks_to_convert, convert_to_h5ad)

    # Gather successful conversions for the next step
    filenames = convert_result.successes
    logger.info(
        f"Parallel conversion complete. {convert_result.num_successes} successes, "
        f"{convert_result.num_failures} failures."
    )

    if not filenames:
        logger.error("No files were successfully converted; exiting early.")
        sys.exit(1)

    # ---------------------------------------------------------------------
    # SERIAL STEP: Create registration mapping
    # ---------------------------------------------------------------------
    logger.info("Creating registration mapping (serial step)...")
    rm = _create_registration_mapping(experiment_uri=am.experiment_path.as_posix(), filenames=filenames)
    with open("rm.pkl", "wb") as f:
        pickle.dump(rm, f)
    logger.info("Registration mapping created and saved to rm.pkl.")

    # ---------------------------------------------------------------------
    # SERIAL STEP: Resize the experiment
    # ---------------------------------------------------------------------
    logger.info("Resizing experiment (serial step)...")
    with open("rm.pkl", "rb") as f:
        rm = pickle.load(f)
    _resize_experiment(am.experiment_path.as_posix(), registration_mapping=rm)
    logger.info("Experiment resized successfully.")

    # ---------------------------------------------------------------------
    # PARALLEL STEP (2): Ingest H5AD files into the SOMA experiment
    # ---------------------------------------------------------------------
    def ingest_h5ad(path):
        """
        Ingest a single H5AD file into the existing experiment using
        the precomputed registration mapping.
        """
        try:
            logger.info(f"[ingest_h5ad] Ingesting {path}")
            _ingest_h5ad_worker(
                experiment_uri=am.experiment_path.as_posix(),
                h5ad_path=path,
                registration_mapping=rm,
            )
            logger.info(f"Successfully ingested {path}")
            return path  # Return the filename on success
        except Exception as e:
            logger.error(f"Failed to ingest {path}: {e}")
            raise

    logger.info("Starting parallel ingestion of H5AD files into experiment...")
    ingest_result = mp_executor.run(filenames, ingest_h5ad)
    logger.info(
        f"Ingestion complete. {ingest_result.num_successes} successes, " f"{ingest_result.num_failures} failures."
    )

    # If needed, you can decide how to handle the failures
    if ingest_result.num_failures > 0:
        logger.warning(f"{ingest_result.num_failures} H5AD ingestion tasks failed.")

    logger.info("All pipeline steps completed. Exiting.")


if __name__ == "__main__":
    main()
