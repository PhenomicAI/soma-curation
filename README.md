# SOMA Curation

## Overview

`soma-curation` is a package designed for bioinformaticians to facilitate their data curation workflows for single cell RNA sequencing data. This README provides instructions on setting up the environment and using Jupytext notebooks to test and develop functions within the package.

The general idea is:

1. Define a schema: Outline the structure of your data up front.
2. Structure your raw data: Organize your data in a consistent, predictable way.
3. Ingest and analyze: Once your data follows the schema, the rest of your pipeline “just works.”

## Inspiration & Acknowledgments

This work was heavily inspired by TileDB-SOMA and the cellxgene census. We extend our gratitude to the TileDB team for their valuable feedback. We hope this package will serve as a helpful community addition by promoting a schema-first approach for TileDB ingestion.

## Setup

Below are setup instructions. If you’re working in VSCode, we highly recommend installing the Python extension.

### Cloning the Repository

Clone this repository to your local machine:

```bash
git clone https://github.com/PhenomicAI/soma-curation.git
cd data_curation/
```

### Creating a Virtual Environment

You only need to create a virtual environment once.

Create and activate a virtual environment:

```bash
virtualenv venv
source venv/bin/activate
pip install .
```
