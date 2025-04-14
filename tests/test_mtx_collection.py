import pytest
from src.soma_curation.collection import MtxCollection
from src.soma_curation.schema import load_schema
from src.soma_curation.constants.create_dummy_structure import create_dummy_mtx_structure


@pytest.fixture
def db_schema():
    return load_schema()


@pytest.fixture
def valid_storage_dir(tmp_path):
    """Create a valid storage directory with no duplicate samples."""
    storage_dir = tmp_path / "valid_storage"
    storage_dir.mkdir()

    # Create a structure with unique sample names across studies
    create_dummy_mtx_structure(storage_dir)
    return storage_dir


@pytest.fixture
def duplicate_sample_storage_dir(tmp_path):
    """Create a storage directory with duplicate sample names across studies."""
    storage_dir = tmp_path / "duplicate_storage"
    storage_dir.mkdir()

    # Create study1 with sample1 and sample2
    study1 = storage_dir / "study1"
    mtx1 = study1 / "mtx" / "sample1"
    mtx1.mkdir(parents=True)
    mtx2 = study1 / "mtx" / "sample2"
    mtx2.mkdir(parents=True)

    # Create study2 with sample1 and sample3
    study2 = storage_dir / "study2"
    mtx3 = study2 / "mtx" / "sample1"  # Duplicate with study1
    mtx3.mkdir(parents=True)
    mtx4 = study2 / "mtx" / "sample3"
    mtx4.mkdir(parents=True)

    # Create study3 with sample2 and sample3
    study3 = storage_dir / "study3"
    mtx5 = study3 / "mtx" / "sample2"  # Duplicate with study1
    mtx5.mkdir(parents=True)
    mtx6 = study3 / "mtx" / "sample3"  # Duplicate with study2
    mtx6.mkdir(parents=True)

    return storage_dir


def test_valid_mtx_collection(valid_storage_dir, db_schema):
    """Test that MtxCollection works with valid storage directory (no duplicate samples)."""
    collection = MtxCollection(storage_directory=valid_storage_dir, db_schema=db_schema)
    assert collection is not None
    assert len(collection.list_studies()) > 0

    # Verify that all samples are unique across studies
    sample_to_studies = {}
    for study in collection.list_studies():
        for sample in collection.list_samples(study):
            if sample in sample_to_studies:
                sample_to_studies[sample].append(study)
            else:
                sample_to_studies[sample] = [study]

    # Check that no sample appears in more than one study
    for sample, studies in sample_to_studies.items():
        assert len(studies) == 1, f"Sample {sample} appears in multiple studies: {studies}"


def test_duplicate_samples_raises_error(duplicate_sample_storage_dir, db_schema):
    """Test that MtxCollection raises an error when duplicate sample names are found."""
    with pytest.raises(ValueError) as exc_info:
        MtxCollection(storage_directory=duplicate_sample_storage_dir, db_schema=db_schema)

    error_message = str(exc_info.value)

    # Check that the error message contains the expected content
    assert "Found duplicate sample names across studies" in error_message

    # Check that all duplicate samples are mentioned
    assert "sample1" in error_message
    assert "sample2" in error_message
    assert "sample3" in error_message

    # Check that the error message shows which studies each duplicate appears in
    assert "study1" in error_message
    assert "study2" in error_message
    assert "study3" in error_message
