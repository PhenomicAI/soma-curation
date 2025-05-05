import argparse
from pathlib import Path
from cloudpathlib import AnyPath
from typing import List, Set
import pandas as pd

from ..sc_logging import logger, _set_level, init_worker_logging
from ..collection import MtxCollection
from ..executor.executors import MultiprocessingExecutor


def get_genes_from_sample(collection: MtxCollection, study: str, sample: str) -> Set[str]:
    """Get unique genes from a single sample.
    
    Args:
        collection: MtxCollection instance
        study: Study name
        sample: Sample name
        
    Returns:
        Set of unique genes found in the sample
    """
    try:
        _, _, features_df = collection.read_mtx(
            collection.storage_directory / study / "mtx" / sample,
            files=["features.tsv.gz"]
        )
        output_genes = set(features_df["gene"].tolist())
        logger.info(f"Found {len(output_genes)} genes in {study}/{sample}")
        return output_genes
    except Exception as e:
        logger.warning(f"Failed to read features for {study}/{sample}: {e}")
        return set()


def main():
    parser = argparse.ArgumentParser(description="Generate a list of all unique genes from a collection.")
    parser.add_argument(
        "--raw-storage-dir",
        type=str,
        required=True,
        help="Directory containing the raw data collection",
    )
    parser.add_argument(
        "--output-file",
        type=str,
        required=True,
        help="File to store the output gene list",
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
    parser.add_argument(
        "--processes",
        type=int,
        default=4,
        help="Number of processes to use for parallel processing",
    )
    args = parser.parse_args()

    # Set up logging
    log_dir = Path(args.log_dir)
    log_dir.mkdir(parents=True, exist_ok=True)
    _set_level(level=10, add_file_handler=True, log_dir=log_dir, log_file="generate_gene_list.log")

    # Create collection
    collection = MtxCollection(
        storage_directory=args.raw_storage_dir,
        include=args.include_studies
    )

    # Prepare tasks for parallel processing
    tasks = []
    for study in collection.list_studies():
        for sample in collection.list_samples(study):
            tasks.append((collection, study, sample))

    # Process samples in parallel
    logger.info(f"Starting parallel processing of {len(tasks)} samples using {args.processes} processes...")
    executor = MultiprocessingExecutor(
        processes=args.processes, 
        init_worker_logging=init_worker_logging, 
        init_args=(10, AnyPath(log_dir).as_posix(), "generate_gene_list.log")
    )
    result = executor.run(tasks, get_genes_from_sample)

    if not result.all_successful:
        logger.warning(f"Failed to process {result.num_failures} samples")

    # Combine results
    all_genes = set()
    for gene_set in result.successes:
        all_genes.update(gene_set)

    unique_genes = sorted(list(all_genes))
    logger.info(f"Found {len(unique_genes)} unique genes in collection")

    # Write to output file
    output_path = AnyPath(args.output_file)
    logger.info(f"Writing gene list to {output_path}")
    
    df = pd.DataFrame(unique_genes, columns=["gene"])
    df.to_csv(str(output_path), sep="\t", header=False, index=False, compression="gzip")
    
    logger.info("Gene list generation complete")


if __name__ == "__main__":
    main() 