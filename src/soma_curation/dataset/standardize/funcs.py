import scipy.sparse as sp
import anndata as ad
import numpy.typing as npt
import numpy as np
import pandas as pd
import sys
import re

from ...sc_logging import logger


def normalize_raw_array(X: sp.spmatrix, coeff: float = 10000) -> sp.spmatrix:
    """
    Normalize a sparse matrix of raw counts.

    Parameters:
    - X: sp.spmatrix
        Sparse matrix of raw counts to be normalized.
    - coeff: float, default=10000
        Coefficient used to scale the normalized counts.

    Returns:
    - sp.spmatrix
        Normalized and log-transformed sparse matrix.
    """
    logger.info("Normalizing raw array...")
    X = X.multiply(1 / X.sum(axis=1)).multiply(coeff).log1p().tocsr()
    return X


def compute_nnz(adata: ad.AnnData) -> npt.NDArray[np.int32]:
    """
    Compute the number of non-zero elements per cell.

    Parameters:
    - adata: ad.AnnData
        AnnData object containing the data matrix.

    Returns:
    - npt.NDArray[np.int32]
        Array containing the number of non-zero elements for each cell.
    """
    logger.info("Computing nnz...")
    return adata.X.getnnz(axis=1)


def compute_umi_counts(adata: ad.AnnData) -> npt.NDArray[np.float64]:
    """
    Compute the UMI counts per cell.

    Parameters:
    - adata: ad.AnnData
        AnnData object containing the data matrix.

    Returns:
    - npt.NDArray[np.float64]
        Array containing the UMI counts for each cell.
    """
    logger.info("Computing umi counts...")
    return np.asarray(adata.X.sum(axis=1)).flatten()


def compute_pct_mito(adata: ad.AnnData) -> npt.NDArray[np.float64]:
    """
    Compute the percentage of mitochondrial gene counts per cell.

    Parameters:
    - adata: ad.AnnData
        AnnData object containing the data matrix.

    Returns:
    - npt.NDArray[np.float64]
        Array containing the percentage of mitochondrial gene counts for each cell.
    """
    logger.info("Computing pct mito...")
    mito_mask = adata.var["gene"].str.startswith("MT-")
    return np.asarray((adata.X[:, mito_mask].sum(axis=1) / adata.X.sum(axis=1)).flatten() * 100).flatten()


def compute_pct_ribo(adata: ad.AnnData) -> npt.NDArray[np.float64]:
    """
    Compute the percentage of ribosomal protein gene counts per cell.

    Parameters:
    - adata: ad.AnnData
        AnnData object containing the data matrix.

    Returns:
    - npt.NDArray[np.float64]
        Array containing the percentage of ribosomal protein gene counts for each cell.
    """
    logger.info("Computing pct ribo...")
    ribo_mask = adata.var["gene"].str.startswith(("RPS", "RPL"))
    return np.asarray((adata.X[:, ribo_mask].sum(axis=1) / adata.X.sum(axis=1)).flatten() * 100).flatten()


def get_main_gene_name_from_2_pd_cols(
    input_genes: pd.DataFrame,
    gene_aliases: pd.DataFrame,
) -> pd.DataFrame:
    """
    Entering an input Pandas dataframe (input_genes) with gene names and/or ENSEMBLE id's (ENS*),
    this function obtains the main gene name (GeneName_main) for each row in input_genes using
    a gene identifier lookup dataframe (gene_aliases).

    Each row in input_genes will be assigned a score (1 to 5) representing the score of
    gene id mapping as follows (1=lowest, 5=highest):
    (5) if the <query gene name> matches a <subject GeneName_main>, then the output is such GeneName_main
    (4) else, if the <query ENS code> matches a <subject ENS_main>, then the output is its associated GeneName_main
    (3) else, if the <query ENS code> matches a <subject ENS_alias>, then the output is its associated GeneName_main
    (2) else, if the <query gene name> matches a <subject GeneName_alias>, then the output is its associated GeneName_main
    (1) else, the output is the <query gene name>
    Notes:
    a) For each group of rows mapping to the same GeneName_main, the row with the highest score will be assigned a keep=True flag
    and the other rows will be assigned keep=False flags. This will avoid having redundant main gene names in downstream applications.
    b) Uppercase is used to compare queries vs. subjects

    Parameters
    ----------
    input_genes: (pd.DataFrame): Pandas dataframe with two columns: i) gene names or ENS codes, ii) gene names, like:
        ENSMUSG00000051951  Xkr4
        ENSMUSG00000089699  Gm1992
        ENSMUSG00000102331  Gm19938
        OR
        Xkr4     Xkr4
        Gm1992   Gm1992
        Gm19938  Gm19938

    gene_aliases: (pd.DataFrame): Pandas dataframe with first 4 columns being:
        GeneName_main  ENS_main         GeneName_aliases  ENS_aliases
        ALG9           ENSG00000086848  ALG9;DIBD1        ENSG00000086848;ENSG00000258529
        ANX8           ENSG00000265190  ANX8              ENSG00000165390;ENSG00000265190
        AQP1           ENSG00000240583  AQP1;CHIP28       ENSG00000240583
        Notes: 1) gene identifier aliases in columns 3 and 4 must be separated by ';'
               2) Human and mouse files can be obtained from: s3://pai-scrnaseq/sctx_gui/gene_level_metadata/

    Returns
    -------
        input_genes pd.Dataframe with added columns 'score', 'GeneName_main' and 'keep'
    """

    ## Load gene aliases
    print("Load gene aliases")
    if len(gene_aliases.columns) == 4:
        gene_aliases = gene_aliases.iloc[:, range(0, 4)].apply(lambda x: x.str.upper())
        gnamemain_dict = {}
        ensmain_to_gnamemain_dict = {}
        genealt_to_gnamemain_dict = {}
        ensalt_to_gnamemain_dict = {}
        for idx, fields in gene_aliases.iterrows():
            gnamemain = fields[0]
            ensmain = fields[1]
            gnamealts = fields[2].split(";")
            ensalts = fields[3].split(";")
            gnamemain_dict[gnamemain] = 1
            ensmain_to_gnamemain_dict[ensmain] = gnamemain
            for gnamealt in gnamealts:
                genealt_to_gnamemain_dict[gnamealt] = gnamemain
            for ensalt in ensalts:
                ensalt_to_gnamemain_dict[ensalt] = gnamemain
    else:
        kill_message = "\nERROR!!! unexpectd number of columns in gene_aliases"
        sys.exit(kill_message)

    ## Mapping query vs. subject gene identifiers
    print("Mapping query vs. subject gene identifiers")
    if len(input_genes.columns) == 2:
        input_genes = input_genes.iloc[:, range(0, 2)].apply(lambda x: x.str.upper())
        for idx, fields in input_genes.iterrows():
            ensorig = str(fields[0])
            gnameorig = str(fields[1])
            gnamenew = ""
            ### Remove trailing version to ENS codes
            ENSversions = re.compile(r"^ENSG\d+\.\d+$|^ENSMUSG\d+\.\d+$")
            if bool(ENSversions.match(ensorig)) == True:
                ensorig = re.sub(r"\.[0-9]+$", "", ensorig)
            if gnameorig in gnamemain_dict.keys():
                score = 5
                gnamenew = gnameorig
            elif ensorig.startswith("ENS") == True:
                if ensorig in ensmain_to_gnamemain_dict.keys():
                    score = 4
                    gnamenew = ensmain_to_gnamemain_dict[ensorig]
                elif ensorig in ensalt_to_gnamemain_dict.keys():
                    score = 3
                    gnamenew = ensalt_to_gnamemain_dict[ensorig]
            if gnamenew == "":
                if gnameorig in genealt_to_gnamemain_dict.keys():
                    score = 2
                    gnamenew = genealt_to_gnamemain_dict[gnameorig]
                else:
                    score = 1
                    gnamenew = gnameorig
            input_genes.at[idx, "score"] = score
            input_genes.at[idx, "GeneName_main"] = gnamenew
    else:
        kill_message = "\nERROR!!! unexpectd number of columns in gene_aliases"
        sys.exit(kill_message)

    ## Select a gene for each group of genes with the same GeneName_main
    input_genes["idx"] = input_genes.index.values
    input_genes_dic = input_genes.groupby("GeneName_main")["idx"].agg(list).to_dict()
    input_genes["keep"] = False
    for gnamenew in input_genes_dic:
        if len(input_genes_dic[gnamenew]) > 1:
            idxs = input_genes_dic[gnamenew]
            topHit_idx = input_genes.iloc[idxs][["score"]].sort_values(by="score", ascending=False).index[0]
            input_genes.at[topHit_idx, "keep"] = True
        else:
            idx = input_genes_dic[gnamenew][0]
            input_genes.at[idx, "keep"] = True

    input_genes = input_genes[[0, 1, "score", "GeneName_main", "keep"]]
    return input_genes
