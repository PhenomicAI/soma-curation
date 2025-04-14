import pytest
import anndata as ad
import numpy as np
import pandas as pd
from src.soma_curation.collection.h5ad_collection import H5adCollection


@pytest.fixture
def h5ad_storage_dir(tmp_path):
    """Create a test directory with h5ad files."""
    storage_dir = tmp_path / "h5ad_storage"
    storage_dir.mkdir()

    # Create a few test h5ad files
    for i in range(3):
        # Create a simple AnnData object
        adata = ad.AnnData(
            X=np.random.rand(10, 5),
            obs=pd.DataFrame({"cell_type": [f"cell_{j}" for j in range(10)]}),
            var=pd.DataFrame({"gene_name": [f"gene_{j}" for j in range(5)]}, index=[f"gene_{j}" for j in range(5)]),
        )
        # Save it to a file
        adata.write(storage_dir / f"test_file_{i}.h5ad")

    # Create a non-h5ad file
    with open(storage_dir / "not_an_h5ad.txt", "w") as f:
        f.write("This is not an h5ad file")

    return storage_dir


def test_list_h5ad_files(h5ad_storage_dir):
    """Test that H5adCollection correctly lists h5ad files."""
    collection = H5adCollection(storage_directory=h5ad_storage_dir)
    files = collection.list_h5ad_files()

    # Should find our 3 test files
    assert len(files) == 3
    assert set(files) == {"test_file_0.h5ad", "test_file_1.h5ad", "test_file_2.h5ad"}

    # Non-h5ad files should be excluded
    assert "not_an_h5ad.txt" not in files


def test_include_filter(h5ad_storage_dir):
    """Test that include parameter correctly filters files."""
    include_files = ["test_file_0.h5ad", "test_file_2.h5ad"]
    collection = H5adCollection(storage_directory=h5ad_storage_dir, include=include_files)
    files = collection.list_h5ad_files()

    assert len(files) == 2
    assert set(files) == set(include_files)
    assert "test_file_1.h5ad" not in files


def test_get_h5ad_path(h5ad_storage_dir):
    """Test that get_h5ad_path returns correct paths."""
    collection = H5adCollection(storage_directory=h5ad_storage_dir)

    for i in range(3):
        filename = f"test_file_{i}.h5ad"
        path = collection.get_h5ad_path(filename)
        assert path == h5ad_storage_dir / filename
