# Quick Start Guide

This guide will help you get started with SOMA Curation quickly. We'll walk through the basic workflow of setting up a schema, organizing data, and creating a collection.

## 1. Define Your Schema

The simplest way to get started is to use the default schema:

```python
from soma_curation.schema import load_schema

# loads the default schema
schema = load_schema()
```

For custom schemas, you can create your own:

```python
from soma_curation.schema import DatabaseSchema

custom_schema = DatabaseSchema(
    PAI_SCHEMA_VERSION="1.0",
    PAI_OBS_SAMPLE_COLUMNS=[
        ("sample_name", "string"),
        ("study_name", "string"),
        # Add more columns as needed
    ],
    PAI_OBS_CELL_COLUMNS=[
        ("barcode", "string"),
    ]
)
```

## 2. Organize Your Data

You can quickly create a test data structure using the provided utility:

```python
from soma_curation.utils import test_dummy_structure
test_dummy_structure()
```

This will create a directory structure like:

```
raw_data/
├── study_1/
│   ├── mtx/
│   │   ├── sample_1/
│   │   │   ├── matrix.mtx.gz
│   │   │   ├── barcodes.tsv.gz
│   │   │   └── features.tsv.gz
│   │   └── sample_2/
│   ├── cell_metadata/
│   │   └── study_1.tsv.gz
│   └── sample_metadata/
│       └── study_1.tsv.gz
└── study_2/
    └── ...
```

## 3. Create a Collection and Atlas

```python
from soma_curation.mtx_collection import MtxCollection
from soma_curation.atlas.crud import AtlasManager

# Create a collection to manage raw data
collection = MtxCollection(
    storage_directory="path/to/raw_data",
    db_schema=schema
)

# Create an atlas for processed data
am = AtlasManager(
    atlas_name="my_atlas",
    db_schema=schema,
    storage_directory="path/to/atlas"
)
am.create()
```

## 4. Process Your Data

### Using Python API

```python
from soma_curation.dataset.anndataset import AnnDataset

# Create and standardize a dataset
dataset = AnnDataset(
    atlas_manager=am,
    collection=collection
)
dataset.standardize()
```

## 5. Ingest Data into Atlas

After standardizing your data, you can access the AnnData object and ingest it into your atlas with a single function call from `tiledbsoma.io`:

```python
import tiledbsoma

# Access the AnnData object from the dataset and ingest it into the atlas
tiledbsoma.io.from_anndata(
    experiment_uri=str(am.experiment_path),
    measurement_name="RNA",
    anndata=dataset.artifact
)

print("Data successfully ingested into atlas!")
```

### Using Command Line Interface

For large-scale data processing, use the multiprocessing ingestion script:

```bash
python -m soma_curation.scripts.multiprocessing_ingest \
    --processes 4 \
    --atlas-name my_atlas \
    --raw-storage-dir path/to/raw_data \
    --h5ad-storage-dir path/to/h5ads \
    --atlas-storage-dir path/to/atlas \
    --log-dir path/to/logs
```

Key command line options:

- `--processes`: Number of worker processes (default: 4)
- `--atlas-name`: Name of the atlas
- `--raw-storage-dir`: Directory containing raw data
- `--h5ad-storage-dir`: Directory for H5AD file storage
- `--atlas-storage-dir`: Directory for atlas storage
- `--log-dir`: Directory for logs
- `--db-schema-fp`: Optional path to custom schema file

## Next Steps

- Learn more about [Schema Definition](user_guide/schema.md)
- Explore [Data Organization](user_guide/data_organization.md)
- Check out [Advanced Usage](examples/advanced.md) for more complex scenarios
