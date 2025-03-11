import pytest
import tiledbsoma as soma
from src.soma_curation.atlas.crud import AtlasManager
from src.soma_curation.schema import load_schema


@pytest.fixture
def global_schema():
    return load_schema()


@pytest.fixture
def atlas_manager(global_schema, tmp_path):
    return AtlasManager(atlas_name="test_atlas", db_schema=global_schema, storage_directory=tmp_path)


def test_create(atlas_manager: AtlasManager):
    atlas_manager.create()

    # Check if the experiment was created
    exp_path = atlas_manager.experiment_path
    assert exp_path.exists()

    # Check metadata from the experiment
    with soma.Experiment.open(exp_path.as_posix()) as exp:
        assert exp.metadata["created_on"] is not None
        assert exp.metadata["pai_schema_version"] == atlas_manager.db_schema.PAI_SCHEMA_VERSION
        assert exp.metadata["atlas_name"] == atlas_manager.atlas_name

    # Check that the obs schema is created
    obs_uri = exp_path / "obs"
    with soma.DataFrame.open(obs_uri.as_posix()) as arr:
        assert arr.schema is not None

    # Check that the var schema is created
    var_uri = exp_path / "ms" / "RNA" / "var"
    with soma.DataFrame.open(var_uri.as_posix()) as arr:
        assert arr.schema is not None

    # Check that the row_raw schema is created
    row_raw_uri = exp_path / "ms" / "RNA" / "X" / "row_raw"
    with soma.SparseNDArray.open(row_raw_uri.as_posix()) as arr:
        assert arr.schema is not None


def test_delete_experiment(atlas_manager: AtlasManager):
    atlas_manager.create()

    # Ensure experiment was created
    exp_path = atlas_manager.experiment_path
    assert exp_path.exists()

    # Delete experiment and verify deletion
    atlas_manager.delete()
    assert not exp_path.exists()
