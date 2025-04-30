from pydantic import BaseModel, ConfigDict, model_validator
import pandas as pd
import scipy.sparse as sp
import anndata as ad
import numpy as np

from typing import List, Tuple, Optional, Union, Generator
from fast_matrix_market import mmread
from cloudpathlib import AnyPath

from ..schema import DatabaseSchema
from ..sc_logging import logger
from ..types.path import ExpandedPath


class MtxCollection(BaseModel):
    """A mapping for local or S3 based raw data for Phenomic

    This class assumes the following directory structure for raw data:

    <storage_directory>
    |_ study_0 (required)
        |_ cell_metadata (not required to specify)
            |_ external_slyper_natmed_2020_32405060.tsv.gz
        |_ sample_metadata (not required to specify)
            |_ external_slyper_natmed_2020_32405060.tsv.gz
        |_ mtx (required)
            |_ GSM4186956
            |_ GSM4186995
    |_ study_1
        ...

    - It also places hard assumptions on the MTX file contents; refer to `read_mtx` method
    - Sample/cell columns come from the schema provided
    - Every sub file is gzipped (eg: `.tsz.gz` or `.mtx.gz`)
    - Tab separated for dataframe-esque files (i.e. `.tsz.gz`)

    Args:
        storage_directory (DirectoryPath): Path to the directory
        schema (DatabaseSchema): Schema to read from
        include (Optional[List[str]]): A list of studies to include. Useful when you want a few studies to be added on

    Returns:
        _type_: CollectionSpec
    """

    storage_directory: ExpandedPath
    db_schema: DatabaseSchema
    include: Optional[List[str]] = None

    model_config = ConfigDict(arbitrary_types_allowed=True)

    @model_validator(mode="after")
    def check_duplicate_samples(self) -> "MtxCollection":
        """Validator to ensure no sample names overlap among studies."""
        sample_to_studies = {}  # Track which studies each sample appears in

        for study in self.list_studies():
            for sample in self.list_samples(study):
                if sample in sample_to_studies:
                    sample_to_studies[sample].append(study)
                    logger.warning(f"Sample '{sample}' appears in studies: {sample_to_studies[sample]}")
                else:
                    sample_to_studies[sample] = [study]

        duplicate_samples = {sample: studies for sample, studies in sample_to_studies.items() if len(studies) > 1}
        if duplicate_samples:
            raise ValueError(
                f"Found duplicate sample names across studies: {duplicate_samples}. "
                "Each sample name must be unique across all studies."
            )

        return self

    def list_studies(self) -> List[str]:
        all_studies = [study_path.parts[-1] for study_path in self.clean_listdir(self.storage_directory)]
        if self.include:
            filtered_studies = [study for study in all_studies if study in self.include]
            logger.info(f"Filtered studies using include parameter: {filtered_studies}")
            return filtered_studies
        return all_studies

    def list_samples(self, study_name: str) -> List[str]:
        study_path = self.storage_directory / study_name / "mtx"
        return [sample_path.parts[-1] for sample_path in self.clean_listdir(study_path)]

    def get_sample_metadata(self, study_name: str) -> pd.DataFrame:
        logger.info(f"Reading sample metadata for {study_name} from {self.storage_directory}...")
        sample_metadata_path = self.storage_directory / study_name / "sample_metadata" / f"{study_name}.tsv.gz"

        return self.read_metadata_file(
            sample_metadata_path, reindex_columns=[x[0] for x in self.db_schema.PAI_OBS_SAMPLE_COLUMNS]
        )

    def get_cell_metadata(self, study_name: str) -> pd.DataFrame:
        logger.info(f"Reading cell metadata for {study_name} from {self.storage_directory}...")
        cell_metadata_path = self.storage_directory / study_name / "cell_metadata" / f"{study_name}.tsv.gz"

        return self.read_metadata_file(
            cell_metadata_path, reindex_columns=[x[0] for x in self.db_schema.PAI_OBS_CELL_COLUMNS]
        )

    def get_mtx(
        self, study_name: str, sample_name: str
    ) -> Tuple[Optional[sp.coo_matrix], Optional[pd.DataFrame], Optional[pd.DataFrame]]:
        logger.info(f"Reading mtx for study: {study_name}, sample: {sample_name} from {self.storage_directory}...")
        mtx_dir = self.storage_directory / study_name / "mtx" / sample_name

        return self.read_mtx(mtx_dir)

    def get_anndata(
        self, study_name: str, sample_name: str, add_cell_metadata: bool = True, add_sample_metadata: bool = True
    ) -> ad.AnnData:
        logger.info(
            f"Assembling AnnData for study: {study_name}, sample: {sample_name} from {self.storage_directory}..."
        )
        mtx, barcodes, features = self.get_mtx(study_name=study_name, sample_name=sample_name)
        if add_cell_metadata:
            cell_metadata = self.get_cell_metadata(study_name=study_name)
            barcodes = self.add_metadata_to_df(
                dataframe=barcodes,
                metadata_df=cell_metadata,
                join=["barcode"],
                columns_to_add=[x[0] for x in self.db_schema.PAI_OBS_CELL_COLUMNS],
            )
        if add_sample_metadata:
            sample_metadata = self.get_sample_metadata(study_name=study_name)
            barcodes = self.add_metadata_to_df(
                dataframe=barcodes,
                metadata_df=sample_metadata,
                join=["sample_name"],
                columns_to_add=[x[0] for x in self.db_schema.PAI_OBS_SAMPLE_COLUMNS],
            )

        barcodes.index = barcodes["barcode"].astype(str)
        barcodes.index.name = "index"
        anndata = ad.AnnData(X=mtx, obs=barcodes, var=features)

        return anndata

    @staticmethod
    def clean_listdir(
        path: AnyPath, ignore_patterns: List[str] = [".DS_Store", ".log"]
    ) -> Generator[AnyPath, None, None]:
        """
        A generator that yields directory contents, ignoring specific files or patterns.

        Args:
            path (AnyPath): The directory to list contents from.
            ignore_patterns (List[str]): List of patterns to ignore (e.g., ".DS_Store").

        Yields:
            AnyPath: Valid directory components.
        """
        for item in path.iterdir():
            if not any(pattern in item.parts[-1] for pattern in ignore_patterns):
                yield item

    def read_mtx(
        self, root_fp: AnyPath, files: Optional[List[str]] = None
    ) -> Tuple[Optional[sp.csr_matrix], Optional[pd.DataFrame], Optional[pd.DataFrame]]:
        matrix, barcodes_df, features_df = None, None, None
        if files is None:
            files = ["matrix.mtx.gz", "barcodes.tsv.gz", "features.tsv.gz"]
        file_paths = [root_fp / fp for fp in files]
        for fp in file_paths:
            if "matrix.mtx.gz" in str(fp):
                matrix = self.mmread(fp).tocsr()
            elif "barcodes.tsv.gz" in str(fp):
                barcodes_df = self.read_csv(fp, sep="\t", names=["barcode"])
                sample_name = fp.parent.parts[-1]
                barcodes_df["sample_name"] = sample_name
            elif "features.tsv.gz" in str(fp):
                features_df = self.read_csv(fp, sep="\t", usecols=[0, 1], index_col=0, names=["index", "gene"])

        return matrix, barcodes_df, features_df

    def read_metadata_file(self, fp: AnyPath, reindex_columns: List[str]) -> pd.DataFrame:
        try:
            metadata_df = self.read_csv(fp, sep="\t")
            metadata_df = metadata_df.reindex(columns=reindex_columns, fill_value="Unknown")
        except Exception as e:
            logger.info(f"File not found: {fp}")
            metadata_df = pd.DataFrame(columns=reindex_columns)

        return metadata_df

    @staticmethod
    def read_csv(filepath: AnyPath, **kwargs) -> pd.DataFrame:
        # Long form of reading versus supplying path because that's how CloudPath reads are supported
        # Pandas does not support CloudPaths atm
        # This helps eliminate the if-else for s3 versus normal posix paths
        # https://s3pathlib.readthedocs.io/en/latest/03-S3-Write-API.html#Pandas
        if isinstance(filepath, AnyPath):
            df = pd.read_csv(filepath, **kwargs)
        else:
            raise ValueError("Unsupported filepath type. Filepath needs to be cloudpathlib AnyPath.")
        return df

    @staticmethod
    def mmread(filepath: AnyPath) -> sp.csr_matrix:
        if not isinstance(filepath, AnyPath):
            raise ValueError("Unsupported filepath type. Filepath needs to be cloudpathlib AnyPath.")
        matrix = mmread(filepath).T.tocsr()

        return matrix

    @staticmethod
    def add_metadata_to_df(
        dataframe: pd.DataFrame,
        metadata_df: pd.DataFrame,
        join: List[str],
        columns_to_add: List[str],
    ) -> pd.DataFrame:
        """Adds cell/sample-level metadata to DataFrame object. Size and order of the original dataframe is maintained.
        N/A values are left as is.

        Args:
            dataframe (pd.DataFrame): DataFrame object with obs attribute to add cell metadata to
            metadata_df (pd.DataFrame): Dataframe containing columns as barcodes, and associated columns of cell metadata
            join (List[str]): Choose what key to join the dataframe to the observation column in the adata object. Eg: `barcode` or `sample_name`
            columns_to_add (List[str]): Columns to add in metadata dataframe. These columns must be present in `cell_metadata_df`. Do not include join columns here!

        Raises:
            ValueError: Mismatch in shape of original adata obs and merged obs if merge is incorrect

        Returns:
            pd.DataFrame: DataFrame object with cell metadata added from dataframe
        """
        original_columns = list(dataframe.columns)

        for col in join:
            if not (col in metadata_df.columns and col in dataframe.columns):
                raise ValueError(f"{col} not in both dataframes...")
            if col in columns_to_add:
                logger.warning(f"{col} specified in `join` and `columns_to_add`, this should be handled gracefully...")

        # Deals with dropping columns that are in both dataframes but not mentioned in `join`
        for col in columns_to_add:
            if col in original_columns and col not in join:
                dataframe.drop(columns=col, inplace=True, axis=1)

        merged_obs = dataframe.merge(metadata_df, how="left", on=join).reindex(
            list(set(original_columns + columns_to_add)), axis=1
        )

        if merged_obs.shape[0] != dataframe.shape[0]:
            raise ValueError("Mismatch in shape of original dataframe and merged dataframe")

        return merged_obs

    def presence_matrix(self, study_name: str, sample_name: str, global_var_list: List[str]) -> sp.coo_matrix:
        """Given a list of global features, return a dataframe with the presence of each feature in the study/sample

        Args:
            study_name (str): Name of the study
            sample_name (str): Name of the sample
            global_var_list (List[str]): List of global features. The presence matrix is returned in sorted order of this list.

        Returns:
            sp.coo_matrix: Presence matrix with the presence of each feature in the study-sample
        """
        root_fp = self.storage_directory / study_name / "mtx" / sample_name
        _, _, features_df = self.read_mtx(root_fp, files=["features.tsv.gz"])

        presence_matrix = np.zeros((1, len(global_var_list)))
        feature_set = set(features_df["gene"])
        for i, feature in enumerate(global_var_list):
            if feature in feature_set:
                presence_matrix[:, i] = 1

        return sp.coo_matrix(presence_matrix)
