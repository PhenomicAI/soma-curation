# src/soma_curation/constants/constants.py
import importlib.resources
import pandas as pd
import anndata as ad
import numpy as np
import scipy.sparse as sp


def dummy_anndata() -> ad.AnnData:
    """
    Create a minimal dummy AnnData object using the dummy_core_geneset.tsv.gz file.

    Returns:
    -------
    ad.AnnData
        An AnnData object with minimal .obs, .var, and .X to satisfy the schema.
    """

    dummy_path = importlib.resources.files("soma_curation.constants").joinpath("dummy_core_geneset.tsv.gz")
    genes = pd.read_csv(dummy_path, sep="\t", header=None)[0].tolist()

    obs = pd.DataFrame(
        {
            "barcode": ["cell1", "cell2", "cell3"],
            "sample_name": ["sample1000", "sample1000", "sample1000"],
            "study_name": ["study1000", "study1000", "study1000"],
        }
    )
    var = pd.DataFrame({"gene": genes}, index=genes)
    data = np.array(
        [
            [1, 0, 2],
            [0, 3, 1],
            [2, 0, 0],
        ]
    )
    X = sp.csr_matrix(data)

    adata = ad.AnnData(X=X, obs=obs, var=var)

    return adata
