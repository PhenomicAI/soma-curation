# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.16.4
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# # Create an Atlas
#
# This notebook demonstrates how to create a Phenomic SOMA atlas with a predefined schema, create a `Dataset` that adheres to that schema, and append information to it.

# +
import data_curation.logging
from data_curation.soma.schema import get_schema
from data_curation.soma.atlas.crud import AtlasManager
from data_curation.soma.dataset.anndataset import AnnDataset

# Set the log level to info
data_curation.logging.info()
# -

# ## Pull a Schema

db_schema = get_schema(organism="human", data_type="scrna")

print(db_schema)

# +
print("OBS columns:")
print(list(db_schema.PAI_OBS_TERM_COLUMNS.keys()))

print("\nVAR columns:")
print(list(db_schema.PAI_VAR_COLUMNS))
# -
print("Validation Schema:")
print(db_schema.VALIDATION_SCHEMA)


# ## Create an Experiment/Atlas 

# !mkdir ~/data_curation/soma

am = AtlasManager(atlas_name="human_v6.0.5", globals_=db_schema, storage_directory="/Users/dylanmendonca/data_curation/")
am.create()
# am.delete() # to delete

# ## Find an AnnDataset, Standardize, Append

# +
# Assuming that you fetch a cool dataset from somewhere
import scanpy as sc

anndata = sc.datasets.pbmc3k()
# -

# Create a Phenomic dataset using our schema
# You should get an error in this case, since the schema does not allow this
# In this case, we're missing sample_name, disease_name, study_name, gene columns
dataset = AnnDataset(artifact=anndata, database_schema=db_schema)

# Fixing issues
anndata.obs["sample_name"] = "sample"
anndata.obs["study_name"] = "study"
anndata.obs["disease_name"] = "disease"
anndata.var["gene"] = anndata.var.index
anndata.obs["barcode"] = anndata.obs.index

# You met the basic criteria. Great success!
dataset = AnnDataset(artifact=anndata, database_schema=db_schema)

# Standardize the columns
# This will run a general standardization pipeline, which involves:
# (a) Fixing the gene names
# (b) Embedding and predicting cell types
# (c) Normalizing the arrays
dataset.standardize()

am.append_anndatasets([dataset])
