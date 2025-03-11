import sys
import pickle
import multiprocessing
from pathlib import Path
import pytest

from src.soma_curation.atlas.crud import AtlasManager
from src.soma_curation.schema import load_schema


# Dummy result class to simulate successful parallel tasks.
class DummyResult:
    def __init__(self):
        self.successes = ["dummy.h5ad"]
        self.num_successes = 1
        self.num_failures = 0


# Patch the MtxCollection so that list_studies and list_samples return dummy values.
@pytest.fixture
def dummy_mtx_collection(monkeypatch):
    from src.soma_curation.mtx_collection.mtx_collection import MtxCollection

    monkeypatch.setattr(MtxCollection, "list_studies", lambda self: ["study1"])
    monkeypatch.setattr(MtxCollection, "list_samples", lambda self, study_name: ["sample1"])


# Patch the MultiprocessingExecutor.run method to always return a DummyResult.
@pytest.fixture
def patch_executor(monkeypatch):
    from src.soma_curation.executor.executors import MultiprocessingExecutor

    monkeypatch.setattr(MultiprocessingExecutor, "run", lambda self, tasks, func: DummyResult())


# Patch ingestion functions to bypass real processing.
@pytest.fixture
def patch_ingest(monkeypatch):
    # Patch create_registration_mapping to return a dummy mapping.
    monkeypatch.setattr(
        "src.soma_curation.ingest.ingestion_funcs.create_registration_mapping",
        lambda experiment_uri, filenames: {"dummy": "value"},
    )
    # Patch resize_experiment to do nothing.
    monkeypatch.setattr(
        "src.soma_curation.ingest.ingestion_funcs.resize_experiment",
        lambda experiment_uri, registration_mapping: None,
    )
    # The convert_and_std_mtx_to_h5ad and ingest_h5ad_soma functions are invoked via the executor,
    # so we don't need to patch them separately here.


# Patch any function that creates a dummy data structure if needed.
@pytest.fixture
def patch_dummy_structure(monkeypatch):
    # For example, if PipelineConfig or another module calls a function to create dummy raw data structure,
    # you can override it here. If not needed, this fixture can be empty.
    monkeypatch.setattr("src.soma_curation.constants.create_dummy_structure", lambda raw_storage_dir: None)


# Ensure that the multiprocessing start method is set to "spawn".
@pytest.fixture(autouse=True)
def set_spawn():
    if multiprocessing.get_start_method(True) != "spawn":
        multiprocessing.set_start_method("spawn", True)


def test_pipeline_run(tmp_path, monkeypatch, dummy_mtx_collection, patch_executor, patch_ingest, patch_dummy_structure):
    # Create temporary directories for raw data, H5AD files, and atlas storage.
    raw_storage_dir = tmp_path / "raw_storage"
    raw_storage_dir.mkdir()
    h5ad_storage_dir = tmp_path / "h5ads"
    h5ad_storage_dir.mkdir()
    atlas_storage_dir = tmp_path / "atlas_storage"
    atlas_storage_dir.mkdir()

    # Build command-line arguments for the pipeline.
    args = [
        "prog",
        "--processes",
        "1",
        "--atlas-name",
        "test_atlas",
        "--raw-storage-dir",
        str(raw_storage_dir),
        "--h5ad-storage-dir",
        str(h5ad_storage_dir),
        "--atlas-storage-dir",
        str(atlas_storage_dir),
    ]
    monkeypatch.setattr(sys, "argv", args)

    # Import and run the pipeline's main() function.
    from src.soma_curation.scripts.multiprocessing_ingest import main

    main()

    # Verify that the atlas was created.
    am = AtlasManager(atlas_name="test_atlas", storage_directory=str(atlas_storage_dir), db_schema=load_schema())
    assert am.exists(), "Atlas should have been created."

    # Verify that the registration mapping file was created.
    rm_file = Path("rm.pkl")
    assert rm_file.exists(), "Registration mapping file rm.pkl should exist."

    with rm_file.open("rb") as f:
        rm = pickle.load(f)
    assert rm == {"dummy": "value"}, "Registration mapping contents do not match expected dummy value."

    # Clean up the registration mapping file.
    rm_file.unlink()
