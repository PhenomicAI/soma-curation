import pytest
import anndata as ad
import pandas as pd
import numpy as np
import scipy.sparse as sp

from pydantic import ValidationError

from src.soma_curation.dataset.anndataset import AnnDataset
from src.soma_curation.constants.constants import dummy_anndata
from src.soma_curation.schema import load_schema


@pytest.fixture
def database_schema():
    # Load the default DatabaseSchema
    return load_schema()


@pytest.fixture
def faulty_anndata():
    # Create a faulty AnnData with missing required columns (study_name)
    faulty_obs = pd.DataFrame(
        {
            "barcode": ["cell1", "cell2", "cell3"],
            "sample_name": ["sample1", "sample2", "sample3"],
        }
    )
    faulty_var = pd.DataFrame({"gene": ["A", "B", "C"], "soma_joinid": [1, 2, 3]})
    faulty_X = sp.csr_matrix(np.array([[1, 0, 3], [4, 0, 6], [7, 8, 0]]))
    return ad.AnnData(X=faulty_X, obs=faulty_obs, var=faulty_var)


@pytest.fixture
def correct_anndata():
    # Create a valid AnnData
    return dummy_anndata()


@pytest.fixture
def non_normalized_anndata():
    # Create a faulty AnnData with normalized data
    correct_obs = pd.DataFrame(
        {
            "barcode": ["cell1", "cell2", "cell3"],
            "sample_name": ["sample1", "sample2", "sample3"],
            "study_name": ["study1", "study2", "study3"],
            "disease_name": ["disease1", "disease2", "disease3"],
        }
    )
    correct_var = pd.DataFrame(
        {
            "gene": ["gene1", "gene2", "gene3"],
        }
    )
    # Values slightly off to simulate non-normalized data
    correct_X = sp.csr_matrix(np.array([[1.1, 0.2, 3], [4, 0.5, 6], [7, 8, 0]]))
    return ad.AnnData(X=correct_X, obs=correct_obs, var=correct_var)


def test_faulty_anndata(database_schema, faulty_anndata):
    # Expect a ValidationError when required columns are missing.
    with pytest.raises(ValidationError):
        AnnDataset(artifact=faulty_anndata, db_schema=database_schema)


def test_non_normalized_anndata(database_schema, non_normalized_anndata):
    # Expect a ValidationError if the data is not normalized properly.
    with pytest.raises(ValidationError):
        AnnDataset(artifact=non_normalized_anndata, db_schema=database_schema)


def test_correct_anndata(database_schema, correct_anndata):
    # Should create the dataset without raising errors.
    dataset = AnnDataset(artifact=correct_anndata, db_schema=database_schema)
    # Optionally verify that the dataset has the expected number of observations.
    assert dataset.artifact.shape[0] == 3


def test_standardize_function(database_schema, correct_anndata):
    dataset = AnnDataset(artifact=correct_anndata, db_schema=database_schema)
    dataset.standardize()

    # Check if standardization is marked as complete.
    assert dataset.standardized is True

    # Check if expected layers are present.
    assert "col_raw" in dataset.artifact.layers
    assert "col_norm" in dataset.artifact.layers

    # Check if computed columns were added.
    for col in ["nnz", "umi_counts", "pct_mito", "pct_ribo"]:
        assert col in dataset.artifact.obs.columns

    # Verify computed column values.
    np.testing.assert_array_equal(dataset.artifact.obs["nnz"].values, [2, 2, 1])
    np.testing.assert_allclose(dataset.artifact.obs["umi_counts"].values, [3, 4, 2], atol=1e-2)
    np.testing.assert_allclose(dataset.artifact.obs["pct_mito"].values, [0.0, 0.0, 0.0], atol=1e-2)
    np.testing.assert_allclose(dataset.artifact.obs["pct_ribo"].values, [0.0, 0.0, 0.0], atol=1e-2)
