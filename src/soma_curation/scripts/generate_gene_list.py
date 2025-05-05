import argparse
from pathlib import Path
from cloudpathlib import AnyPath

from ..sc_logging import logger, _set_level
from ..collection import MtxCollection
from ..schema import DatabaseSchema, load_schema


def main():
    parser = argparse.ArgumentParser(description="Generate a list of all unique genes from a collection.")
    parser.add_argument(
        "--raw-storage-dir",
        type=str,
        required=True,
        help="Directory containing the raw data collection",
    )
    parser.add_argument(
        "--db-schema-fp",
        type=str,
        required=True,
        help="Filepath for database schema",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        required=True,
        help="Directory to store the output gene list",
    )
    parser.add_argument(
        "--include-studies",
        type=str,
        nargs="+",
        default=None,
        help="Optional list of studies to include. If not provided, all studies will be processed.",
    )
    parser.add_argument(
        "--log-dir",
        type=str,
        default="./logs",
        help="Directory to store logs",
    )
    args = parser.parse_args()

    # Set up logging
    log_dir = Path(args.log_dir)
    log_dir.mkdir(parents=True, exist_ok=True)
    _set_level(level=10, add_file_handler=True, log_dir=log_dir, log_file="generate_gene_list.log")

    # Load schema
    db_schema = load_schema(args.db_schema_fp)

    # Create collection
    collection = MtxCollection(
        storage_directory=args.raw_storage_dir,
        db_schema=db_schema,
        include=args.include_studies
    )

    # Get all genes
    logger.info("Starting to collect all unique genes...")
    all_genes = collection.get_all_genes()
    logger.info(f"Found {len(all_genes)} unique genes")

    # Write to output file
    output_path = AnyPath(args.output_dir) / "all_genes.txt.gz"
    logger.info(f"Writing gene list to {output_path}")
    
    import pandas as pd
    df = pd.DataFrame(all_genes)
    df.to_csv(output_path, sep="\t", header=False, index=False, compression="gzip")
    
    logger.info("Gene list generation complete")


if __name__ == "__main__":
    main() 