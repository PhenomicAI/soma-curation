import pytest
import tiledbsoma as soma
import sys

from src.soma_curation.constants.create_dummy_structure import create_dummy_mtx_structure, create_dummy_h5ad_structure
from src.soma_curation.scripts.compute_feature_presence import (
    determine_sample_df_to_process,
    main as feature_presence_main,
)
from src.soma_curation.scripts.multiprocessing_ingest import main as ingest_main
from src.soma_curation.schema import load_schema


@pytest.fixture(scope="session")
def mock_schema():
    """Create a mock schema for testing."""
    schema = load_schema()
    return schema


@pytest.fixture
def test_dirs(tmp_path):
    """Create test directories."""
    dirs = {
        "raw": tmp_path / "raw",
        "raw_2": tmp_path / "raw_2",
        "h5ad": tmp_path / "h5ads",
        "atlas": tmp_path / "atlas",
        "rm_pkl": tmp_path / "rm.pkl",
        "filenames_pkl": tmp_path / "filenames.pkl",
        "logs": tmp_path / "logs",
    }

    return dirs


@pytest.fixture
def dummy_mtx_structure(test_dirs):
    """Create a dummy MTX structure for testing."""
    create_dummy_mtx_structure(test_dirs["raw"])
    return test_dirs["raw"]


@pytest.fixture
def dummy_h5ad_structure(test_dirs):
    """Create a dummy H5AD structure for testing."""
    create_dummy_h5ad_structure(test_dirs["raw_2"])
    return test_dirs["raw_2"]


@pytest.fixture
def test_experiment_name(test_dirs):
    """Test the experiment name."""
    return "test_atlas"


@pytest.fixture
def test_experiment(test_dirs, dummy_mtx_structure, mock_schema, monkeypatch, test_experiment_name):
    """Create a test experiment using the multiprocessing_ingest script."""
    # Prepare test arguments for ingest
    test_args = [
        "src.soma_curation.scripts.multiprocessing_ingest",
        "--processes",
        "1",
        "--atlas-name",
        test_experiment_name,
        "--raw-storage-dir",
        str(test_dirs["raw"]),
        "--h5ad-storage-dir",
        str(test_dirs["h5ad"]),
        "--atlas-storage-dir",
        str(test_dirs["atlas"]),
        "--log-dir",
        str(test_dirs["logs"]),
        "--registration-mapping-pkl",
        str(test_dirs["rm_pkl"]),
        "--filenames-pkl",
        str(test_dirs["filenames_pkl"]),
    ]

    # Mock sys.argv and run ingest
    monkeypatch.setattr(sys, "argv", test_args)
    ingest_main()
    test_dirs["rm_pkl"].unlink()
    test_dirs["filenames_pkl"].unlink()

    assert not test_dirs["rm_pkl"].exists(), "RM pickle should not exist"
    assert not test_dirs["filenames_pkl"].exists(), "Filenames pickle should not exist"

    return str(test_dirs["atlas"] / "test_atlas")


def test_determine_samples_to_process(test_dirs, mock_schema, test_experiment):
    """Test the determine_samples_to_process function with a new sample."""

    # Run determine_samples_to_process
    samples_df, resize_length = determine_sample_df_to_process(test_experiment, mock_schema)

    # Verify only the new sample is detected
    assert samples_df.shape[0] == 3, "Should detect all samples"
    assert resize_length == 3, "Should resize the presence matrix to the number of samples"


def test_feature_presence_script(
    test_dirs, test_experiment, mock_schema, monkeypatch, test_experiment_name, dummy_h5ad_structure
):
    """Test the full feature presence computation script."""
    # Prepare test arguments
    test_args = [
        "src.soma_curation.scripts.compute_feature_presence",
        "--n-processes",
        "1",
        "--exp-uri",
        test_experiment,
        "--raw-storage-dir",
        str(test_dirs["raw"]),
        "--raw-collection-type",
        "mtx",
    ]

    # Mock sys.argv
    monkeypatch.setattr(sys, "argv", test_args)

    # Run the feature presence computation
    feature_presence_main()

    # Verify results
    with soma.Experiment.open(test_experiment) as exp:
        presence_matrix = exp.ms["RNA"]["feature_presence_matrix"]

        total_samples = exp.obs.read(column_names=["sample_name"]).concat().to_pandas()["sample_name"].nunique()
        assert presence_matrix.shape[0] == total_samples
        assert presence_matrix.shape[1] == len(mock_schema.SORTED_CORE_GENES)
        assert (
            presence_matrix.read().coos().concat().to_scipy().todense()
            == [
                [1, 1, 0],
                [1, 1, 0],
                [1, 1, 0],
            ]
        ).all(), "Feature presence matrix is not correct"
    print("test dirs raw 2", test_dirs["raw_2"])
    test_args = [
        "src.soma_curation.scripts.multiprocessing_ingest",
        "--processes",
        "1",
        "--atlas-name",
        test_experiment_name,
        "--raw-storage-dir",
        str(test_dirs["raw_2"]),
        "--h5ad-storage-dir",
        str(test_dirs["h5ad"]),
        "--atlas-storage-dir",
        str(test_dirs["atlas"]),
        "--raw-collection-type",
        "h5ad",
        "--registration-mapping-pkl",
        str(test_dirs["rm_pkl"]),
        "--filenames-pkl",
        str(test_dirs["filenames_pkl"]),
    ]

    # Mock sys.argv
    monkeypatch.setattr(sys, "argv", test_args)

    ingest_main()

    with soma.Experiment.open(test_experiment) as exp:
        print("ARRAY", exp.ms["RNA"]["feature_presence_matrix"].shape)

    # Run the feature presence computation
    samples_df, resize_length = determine_sample_df_to_process(test_experiment, mock_schema)

    assert len(samples_df) == 1, "Should detect only the new sample"
    assert resize_length == 4, "Should resize the presence matrix to the number of samples"
