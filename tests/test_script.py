import sys
import pickle
import multiprocessing
import pytest
import tiledbsoma.io
from pathlib import Path

from src.soma_curation.atlas.crud import AtlasManager
from src.soma_curation.schema import load_schema
from src.soma_curation.constants.create_dummy_structure import create_dummy_mtx_structure
from src.soma_curation.scripts.multiprocessing_ingest import main as ingest_main
from src.soma_curation.collection import MtxCollection


@pytest.fixture(scope="session")
def schema():
    """Create a mock schema for testing."""
    schema = load_schema()
    return schema


@pytest.fixture(scope="session")
def test_dirs(tmp_path_factory):
    """Create test directories."""
    base = tmp_path_factory.mktemp("test_pipeline")
    dirs = {
        "raw": base / "raw",
        "h5ad": base / "h5ads",
        "atlas": base / "atlas",
        "rm_pkl": base / "rm.pkl",
        "filenames_pkl": base / "filenames.pkl",
        "logs": base / "logs",
    }

    # Create directories
    for dir_path in dirs.values():
        if not str(dir_path).endswith(".pkl"):
            dir_path.mkdir(parents=True, exist_ok=True)

    return dirs


@pytest.fixture(scope="session")
def dummy_mtx_structure(test_dirs):
    """Create a dummy MTX structure for testing."""
    create_dummy_mtx_structure(test_dirs["raw"])
    return test_dirs["raw"]


@pytest.fixture(autouse=True)
def set_spawn():
    """Ensure multiprocessing uses spawn method."""
    if multiprocessing.get_start_method(True) != "spawn":
        multiprocessing.set_start_method("spawn", True)


def verify_h5ad_files(h5ad_dir: Path, collection: MtxCollection):
    """Verify that H5AD files were created correctly."""
    h5ad_files = list(h5ad_dir.glob("*.h5ad"))
    assert len(h5ad_files) > 0, "No H5AD files were created"

    # Check that we have one H5AD file for each sample
    expected_count = sum(len(collection.list_samples(study)) for study in collection.list_studies())
    assert len(h5ad_files) == expected_count, f"Expected {expected_count} H5AD files, got {len(h5ad_files)}"


def test_pipeline_run(test_dirs, dummy_mtx_structure, schema, monkeypatch):
    """Test the full pipeline execution."""
    # Prepare test arguments
    test_args = [
        "src.soma_curation.scripts.multiprocessing_ingest",
        "--processes",
        "1",
        "--atlas-name",
        "test_atlas",
        "--raw-storage-dir",
        str(test_dirs["raw"]),
        "--h5ad-storage-dir",
        str(test_dirs["h5ad"]),
        "--atlas-storage-dir",
        str(test_dirs["atlas"]),
        "--registration-mapping-pkl",
        str(test_dirs["rm_pkl"]),
        "--filenames-pkl",
        str(test_dirs["filenames_pkl"]),
        "--log-dir",
        str(test_dirs["logs"]),
    ]

    # Mock sys.argv
    monkeypatch.setattr(sys, "argv", test_args)

    # Run the pipeline
    ingest_main()

    # Verify the atlas was created
    am = AtlasManager(atlas_name="test_atlas", storage_directory=str(test_dirs["atlas"]), db_schema=schema)
    assert am.exists(), "Atlas should have been created"

    # Verify H5AD files were created
    collection = MtxCollection(storage_directory=test_dirs["raw"], db_schema=schema)
    verify_h5ad_files(test_dirs["h5ad"], collection)

    # Verify registration mapping file was created
    assert test_dirs["rm_pkl"].is_file(), "Registration mapping file should exist"
    assert test_dirs["filenames_pkl"].is_file(), "Filenames pickle file should exist"

    # Verify log directory was created and contains logs
    assert test_dirs["logs"].is_dir(), "Log directory should exist"
    assert any(test_dirs["logs"].iterdir()), "Log directory should not be empty"

    # Load and verify the registration mapping
    with test_dirs["rm_pkl"].open("rb") as f:
        rm = pickle.load(f)
    assert isinstance(rm, tiledbsoma.io.ExperimentAmbientLabelMapping), "Registration mapping should be a dictionary"
    assert rm.get_obs_shape() == 9, "Registration mapping should have 9 observations"
    assert rm.get_var_shapes() == {"RNA": 3}, "Registration mapping should have 3 features"

    # Load and verify the filenames
    with test_dirs["filenames_pkl"].open("rb") as f:
        filenames = pickle.load(f)
    assert isinstance(filenames, list), "Filenames should be a list"
    assert len(filenames) > 0, "Filenames list should not be empty"
