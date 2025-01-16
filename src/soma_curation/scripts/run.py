import argparse
import sys
import pickle

from pathlib import Path

from schema import get_schema
from atlas.crud import AtlasManager
from dataset.anndataset import AnnDataset
from mtx_collection.mtx_collection import MtxCollection
from sc_logging import logger, _set_level


def main():

    parser = argparse.ArgumentParser("pai-data-curation")
    sub_parser = parser.add_subparsers(dest="command")

    # Build parser
    build = sub_parser.add_parser("build")
    build.add_argument(
        "--organism",
        required=True,
        type=str,
        choices=["human", "mouse", "multispecies"],
        help="organism to build the atlas for",
    )
    build.add_argument(
        "--raw_storage_directory",
        required=True,
        type=str,
        help="storage directory pointing to raw data, can be a S3 or local directory. Must abide by Phenomic raw data structure",
    )
    build.add_argument(
        "--atlas_name",
        required=True,
        type=str,
        help="Name of the atlas",
    )
    build.add_argument(
        "--atlas_storage_directory",
        required=True,
        type=str,
        help="The storage directory for the atlas. Must be local for now",
    )
    build.add_argument(
        "--log_filename",
        required=False,
        type=str,
        default="build.log",
        help="The storage directory for the atlas. Must be local for now",
    )

    # Update parser
    update = sub_parser.add_parser("update", parents=[build], add_help=False)
    update.add_argument(
        "--include",
        required=False,
        nargs="+",
        help="List of studies to include in the update process. Example: --include study1 study2",
    )

    args = parser.parse_args()

    # Set logging
    # Create a logs folder where you are
    log_dir = Path(args.atlas_storage_directory) / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    _set_level(level=20, add_file_handler=True, log_dir=log_dir, log_file=args.log_filename)

    if args.command == "build":
        logger.info(
            f"Building atlas for organism={args.organism} using raw data from {args.raw_storage_directory}. \
            Soma Experiment called {args.atlas_name} will be stored {args.atlas_storage_directory}. \
            Logs will be stored in {log_dir / args.log_filename}."
        )
    elif args.command == "update":
        logger.info(
            f"Updating atlas for organism={args.organism} using raw data from {args.raw_storage_directory}. \
            Soma Experiment called {args.atlas_name} will be updated in {args.atlas_storage_directory}. \
            Including only studies: {args.include}. \
            Logs will be stored in {log_dir / args.log_filename}."
        )
    # Get a schema, we can write multiple in the future, eg: human, mouse, cross-species, blah
    db_schema = get_schema(organism=args.organism)

    # Connect to a Phenomic Mtx collection, run help(MtxCollection) to see how this directory should be structured
    mtx_collection = MtxCollection(
        storage_directory=args.raw_storage_directory,
        db_schema=db_schema,
        include=args.include if args.command == "update" else None,
    )

    # Create an "Atlas" object, and then create the atlas via the create command
    am = AtlasManager(atlas_name=args.atlas_name, storage_directory=args.atlas_storage_directory, globals_=db_schema)
    if args.command == "build" and am.exists():
        logger.error(f"{args.atlas_name} exists already. Aborting...")
        sys.exit(1)
    elif args.command == "update" and not am.exists():
        logger.error(f"{args.atlas_name} does not exist. Aborting...")
        sys.exit(1)

    if args.command == "build":
        am.create()

    errors = {}

    for study in mtx_collection.list_studies():
        for sample in mtx_collection.list_samples(study_name=study):
            logger.info(f"Ingesting study: {study}, sample: {sample}...")
            try:
                # Get the AnnData
                logger.info("Get Anndata...")
                adata = mtx_collection.get_anndata(study_name=study, sample_name=sample)

                # Create a `Phenomic` AnnDataset with some checking etc
                logger.info("Converting to Phenomic Anndataset...")
                anndataset = AnnDataset(artifact=adata, database_schema=db_schema)

                # Standardize runs some functions eg: nnz, feature_presence, etc.
                logger.info("Standardize Anndataset...")
                anndataset.standardize()

                # Append to atlas
                logger.info("Append to Atlas...")
                am.append_anndatasets([anndataset])
                logger.info("Successfully ingested study: {study}, sample: {sample}")

            except Exception as e:
                logger.error(f"Error processing study: {study}, sample: {sample}. Detailed error: {e}")
                errors[study] = {}
                errors[study][sample] = e
                continue
    if errors:
        output_error_path = log_dir / "errors.pkl"
        with open(output_error_path, "wb") as f:
            pickle.dump(errors, f)
        logger.error(
            f"There were some errors in your pipeline. A detailed breakdown by sample has been written to {output_error_path}. Open using \
            import pickle \
            with open({output_error_path}, 'rb') as f: \
                errors = pickle.load(f) \
            print(errors)"
        )
    else:
        logger.info(
            f"Success! Your pipeline completed without any errors. Please check out your logs in {log_dir / args.log_filename}"
        )


if __name__ == "__main__":
    main()
