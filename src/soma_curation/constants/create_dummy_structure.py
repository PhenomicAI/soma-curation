import uuid
import gzip
from pathlib import Path
from typing import List, Union

import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
from fast_matrix_market import mmwrite


def create_sample_metadata(file_path: Path, study_name: str, sample_names: List[str]):
    """
    Create a sample metadata TSV (gzipped) file for a given study.
    """
    data = {
        "sample_name": sample_names,
        "study_name": [study_name] * len(sample_names),
        "batch_name": ["Batch1"] * len(sample_names),
    }
    df = pd.DataFrame(data)
    df.to_csv(file_path, sep="\t", index=False, compression="gzip")


def create_cell_metadata_df(barcodes: List[str]) -> pd.DataFrame:
    """
    Create a small cell metadata DataFrame using the provided barcodes.
    """
    data = {
        "barcode": barcodes,
        "standard_true_celltype": ["T-cell", "B-cell", "Monocyte"],
        "authors_celltype": ["T-cell", "B-cell", "Monocyte"],
    }
    return pd.DataFrame(data)


def create_mtx_files(sample_path: Path, barcodes: List[str]):
    """
    Create a dummy matrix file (in Matrix Market format, gzipped), along with
    corresponding barcodes and features files.
    """
    # Create a 3x3 sparse matrix as a dummy matrix.
    matrix = coo_matrix(np.array([[1, 0, 0], [0, 1, 1], [1, 1, 0]]))
    matrix_file = sample_path / "matrix.mtx.gz"
    with gzip.open(matrix_file, "wb") as f:
        mmwrite(f, matrix)

    # Write barcodes.tsv.gz
    barcodes_df = pd.DataFrame({"barcode": barcodes})
    barcodes_df.to_csv(
        sample_path / "barcodes.tsv.gz",
        sep="\t",
        index=False,
        header=False,
        compression="gzip",
    )

    # Write features.tsv.gz (dummy gene features)
    features = pd.DataFrame({"feature": ["A", "B", "C"], "gene": ["D", "A", "G"]})
    features.to_csv(
        sample_path / "features.tsv.gz",
        sep="\t",
        index=False,
        header=False,
        compression="gzip",
    )


def create_test_data_structure(base_path: Path):
    """
    Create a test data structure mimicking a multi-study/single-cell dataset.

    The layout will be as follows for each study:

        study_X/
            sample_metadata/
                study_X.tsv.gz
            cell_metadata/
                study_X.tsv.gz
            mtx/
                sample_Y/
                    matrix.mtx.gz
                    barcodes.tsv.gz
                    features.tsv.gz
    """
    # Define a dummy structure: two studies, with each study containing one or more samples.
    raw_data = {
        "study_0": {
            "sample_0": [uuid.uuid4().hex for _ in range(3)],
            "sample_1": [uuid.uuid4().hex for _ in range(3)],
        },
        "study_1": {"sample_2": [uuid.uuid4().hex for _ in range(3)]},
    }

    for study, sample_dict in raw_data.items():
        study_path = base_path / study
        study_path.mkdir(parents=True, exist_ok=True)

        # Create sample_metadata folder and file.
        sample_metadata_path = study_path / "sample_metadata"
        sample_metadata_path.mkdir(parents=True, exist_ok=True)
        sample_metadata_file = sample_metadata_path / f"{study}.tsv.gz"
        create_sample_metadata(sample_metadata_file, study, list(sample_dict.keys()))

        # For collecting cell metadata from each sample.
        cell_metadata_dfs = []

        # Process each sample in the study.
        for sample, barcodes in sample_dict.items():
            # Create a cell metadata DataFrame for this sample.
            cell_metadata_dfs.append(create_cell_metadata_df(barcodes))

            # Create the MTX directory structure for the sample.
            sample_mtx_path = study_path / "mtx" / sample
            sample_mtx_path.mkdir(parents=True, exist_ok=True)
            create_mtx_files(sample_mtx_path, barcodes)

        # Concatenate the cell metadata for the study and save.
        cell_metadata_path = study_path / "cell_metadata"
        cell_metadata_path.mkdir(parents=True, exist_ok=True)
        cell_metadata_file = cell_metadata_path / f"{study}.tsv.gz"
        pd.concat(cell_metadata_dfs, axis=0).to_csv(cell_metadata_file, sep="\t", index=False, compression="gzip")


def create_dummy_structure(base_path: Union[str, Path]):
    """
    Create a dummy data structure on disk for testing purposes.

    This function wraps the creation of a test data structure (with sample metadata,
    cell metadata, and MTX files) under the specified base directory.

    Parameters
    ----------
    base_path : Path
        The base directory where the dummy structure will be created.
    """
    if isinstance(base_path, str):
        base_path = Path(base_path)
    create_test_data_structure(base_path)
