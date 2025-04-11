# SOMA Curation

[![PyPI version](https://badge.fury.io/py/soma-curation.svg)](https://badge.fury.io/py/soma-curation)
[![TestPyPI version](https://img.shields.io/badge/TestPyPI-pre--releases-blue)](https://test.pypi.org/project/soma-curation/)
[![Documentation Status](https://img.shields.io/github/actions/workflow/status/PhenomicAI/soma-curation/docs.yml?label=docs)](https://github.com/PhenomicAI/soma-curation/actions/workflows/docs.yml)
[![Tests](https://img.shields.io/github/actions/workflow/status/PhenomicAI/soma-curation/test.yml?label=tests)](https://github.com/PhenomicAI/soma-curation/actions/workflows/test.yml)
[![Build](https://img.shields.io/github/actions/workflow/status/PhenomicAI/soma-curation/build.yml?label=build)](https://github.com/PhenomicAI/soma-curation/actions/workflows/build.yml)

## Overview

`soma-curation` is a light-weight Python package used at Phenomic to streamline the curation and management of single-cell RNA sequencing (scRNA-seq) atlases using TileDB-SOMA. It's still in its early stages, but the hope is to allow bioinformaticians and ML practitioners to organize their SOMA atlases and access their raw data a bit better. There are assumptions of raw storage organization baked into the package that mimic practices at Phenomic.

## Installation

```bash
# Create and activate a virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install the package
pip install soma-curation
```

## Quick Start

1. Define your schema:

```python
from soma_curation.schema import load_schema

# loads the default schema
schema = load_schema()
```

2. Organize your raw data:

```python
# You can simulate this structure with the following commands
from soma_curation.utils import test_dummy_structure
test_dummy_structure()
```

It should give you a structure like this:

```
raw_data/
├── study_1/
│   ├── mtx/
│   │   ├── sample_1/
│   │   └── sample_2/
│   ├── cell_metadata/
│   │   ├── study_1.tsv.gz
│   └── sample_metadata/
│   │   ├── study_1.tsv.gz
└── study_2/
    └── ...
```

3. Create and use your collection:

```python
from soma_curation.collection import MtxCollection

# For MTX files
collection = MtxCollection(
    storage_directory="path/to/raw_data",
    db_schema=schema
)

# Access AnnDatas from MTX files
adata = collection.get_anndata(study_name="study_1", sample_name="sample_1")

# For H5AD files
h5ad_collection = H5adCollection(
    storage_directory="path/to/h5ad_files"
)

# List all H5AD files
h5ad_files = h5ad_collection.list_h5ad_files()

# Access AnnData directly from an H5AD file
adata = h5ad_collection.get_anndata(filename="file1.h5ad")
```

4. Create a TileDB-SOMA Experiment

```python
from soma_curation.atlas.crud import AtlasManager

# Create an atlas
am = AtlasManager(atlas_name="...", db_schema=db_schema, storage_directory="...")
am.create()

# Delete an atlas
# am.delete()
```

5. Create a Dataset according to your schema and standardize it

```python
from soma_curation.dataset.anndataset import AnnDataset

# Create a Phenomic Dataset
# Original anndata is stored under the `.artifact` attribute
dataset = AnnDataset(
    atlas_manager=am,
    collection=collection
)
dataset.standardize()
```

6. Ingest your data into the TileDB-SOMA Experiment using traditional TileDB-SOMA syntax documented [here](https://documentation.cloud.tiledb.com/academy/structure/life-sciences/single-cell/tutorials/data-ingestion/)

## Documentation

For detailed documentation, including API reference and usage examples, visit our [documentation site](https://phenomicai.github.io/soma-curation/).

## Contributing

We welcome contributions! Please see our [Contributing Guide](CONTRIBUTING.md) for details.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

This work was inspired by the TileDB-SOMA and CellxGene Census teams. We extend our gratitude to the TileDB team for their valuable feedback and support.

## Setup

Below are setup instructions. If you're working in VSCode, we highly recommend installing the Python extension.

### Cloning the Repository

Clone this repository to your local machine:

```bash
git clone https://github.com/PhenomicAI/soma-curation.git
cd soma-curation/
```

### Developer Setup

You only need to create a virtual environment once.

Create and activate a virtual environment:

```bash
virtualenv venv
source venv/bin/activate
pip install ".[dev]"
```
