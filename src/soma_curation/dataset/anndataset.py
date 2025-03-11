import anndata as ad
import anndata.io
import numpy as np
import scipy.sparse as sp
import pyarrow as pa
import pandas as pd

from typing import List, Tuple, Union
from typing_extensions import Self
from pydantic import BaseModel, ConfigDict, model_validator, Field, computed_field
from pathlib import Path
from cloudpathlib import CloudPath

from ..schema import DatabaseSchema
from ..sc_logging import logger
from .standardize.funcs import normalize_raw_array


class AnnDataset(BaseModel):
    """
    AnnDataset class is designed to validate and standardize AnnData objects according to a given database schema.

    Attributes:
    - artifact: ad.AnnData
        The AnnData object containing the data matrix.
    - db_schema: DatabaseSchema
        The schema defining the structure and validation rules for the data.

    Methods:
    - validate: Validates the AnnData object against the provided schema.
    - standardize: Standardizes the AnnData object by normalizing and computing required columns.
    """

    artifact: ad.AnnData = Field(repr=False)
    db_schema: DatabaseSchema = Field(repr=False)
    standardized: bool = Field(default=False)

    model_config = ConfigDict(arbitrary_types_allowed=True)

    @computed_field
    @property
    def shape(self) -> Tuple[int, int]:
        return self.artifact.shape

    @model_validator(mode="after")
    def validate(self) -> Self:
        """
        Validate the AnnData object.

        Returns:
        - Self
            The validated AnnDataset object.
        """
        errors = []
        self._validate_obs(errors=errors)
        self._validate_var(errors=errors)
        self._validate_X(errors=errors)

        if errors:
            raise ValueError("\n".join(errors))

        return self

    def _validate_obs(self, errors: List[str]):
        """
        Validate the .obs attribute of the AnnData object.

        Parameters:
        - errors: List[str]
            List to collect error messages.
        """
        if not hasattr(self.artifact, "obs"):
            errors.append("AnnData missing .obs attribute")
            return

        assert (
            "barcode" in self.db_schema.VALIDATION_SCHEMA.REQUIRED_OBS_COLUMNS
        ), "barcode needs to be present in required obs columns"
        for col in self.db_schema.VALIDATION_SCHEMA.REQUIRED_OBS_COLUMNS:
            if col not in self.artifact.obs:
                errors.append(f"Missing col in obs: `{col}`")
            else:
                if self.artifact.obs[col].isna().any():
                    errors.append(f"Column `{col}` in obs contains NaN values")
                if (self.artifact.obs[col] == "").any():
                    errors.append(f"Column `{col}` in obs contains empty strings")
                if self.artifact.obs[col].isnull().any():
                    errors.append(f"Column `{col}` in obs contains null values")

    def _validate_var(self, errors: List[str]):
        """
        Validate the .var attribute of the AnnData object.

        Parameters:
        - errors: List[str]
            List to collect error messages.
        """
        assert (
            "gene" in self.db_schema.VALIDATION_SCHEMA.REQUIRED_VAR_COLUMNS
        ), "Gene needs to be present in required var columns"
        if not hasattr(self.artifact, "var"):
            errors.append("AnnData missing .var attribute")
            return

        for col in self.db_schema.VALIDATION_SCHEMA.REQUIRED_VAR_COLUMNS:
            if col not in self.artifact.var:
                errors.append(f"Missing col in var: `{col}`")
            else:
                if self.artifact.var[col].isna().any():
                    errors.append(f"Column `{col}` in var contains NaN values")
                if (self.artifact.var[col] == "").any():
                    errors.append(f"Column `{col}` in var contains empty strings")
                if self.artifact.var[col].isnull().any():
                    errors.append(f"Column `{col}` in var contains null values")

        # Checking intersection of genes is above threshold
        if "gene" not in self.artifact.var:
            errors.append("Gene not in .var attribute")
            return
        genes = set(self.artifact.var["gene"])
        intersection = len(self.db_schema.CORE_GENES.intersection(genes)) / len(self.db_schema.CORE_GENES)
        # Sometimes you might want to explcitly set this to 0 for testing purposes
        if self.db_schema.VALIDATION_SCHEMA.GENE_INTERSECTION_THRESHOLD_FRAC > 0:
            if intersection < self.db_schema.VALIDATION_SCHEMA.GENE_INTERSECTION_THRESHOLD_FRAC:
                errors.append(
                    f"Gene intersection >= {self.db_schema.VALIDATION_SCHEMA.GENE_INTERSECTION_THRESHOLD_FRAC} required"
                )
        else:
            pass

    def _validate_X(self, errors: List[str]):
        """
        Validate the .X attribute of the AnnData object.

        Parameters:
        - errors: List[str]
            List to collect error messages.
        """
        x = None
        if hasattr(self.artifact, "raw") and hasattr(
            self.artifact.raw, "X"
        ):  # CellxGene: raw: adata.raw.X; normalized: adata.X
            logger.info("Found raw.X in AnnData")
            x = self.artifact.raw.X
            logger.info("Setting .X and deleting raw.X")
            self.artifact.X = self.artifact.raw.X
            del self.artifact.raw.X

        elif hasattr(self.artifact, "X"):  # Phenomic: raw: adata.X; normalized: adata.layers["X_norm"]
            logger.info("Found .X in AnnData")
            x = self.artifact.X

        if x is None:
            errors.append("Raw counts are needed within raw.X or .X of AnnData object")
            return

        if not sp.issparse(x):
            errors.append("X matrix is not sparse")
            return

        x_sum = (x.data - x.data.astype(np.int32)).sum()
        if x_sum > 0.0:  # if not integer x_sum >> 1.0
            errors.append("AnnData counts are not integer values")

    def standardize(self) -> bool:
        """
        Standardize the AnnData object by normalizing and computing required columns.

        Returns:
        - bool
            True if standardization is successful.
        """
        logger.info(
            "Running standardization pipeline which will include embedding, cell labeling, computed columns, etc."
        )

        self._standardize_obs()
        self._standardize_var()
        self._standardize_X()
        self._standardize_obsm()
        self.standardized = True

        return True

    def _standardize_obs(self):
        """
        Standardize the .obs attribute of the AnnData object.
        """
        non_index_columns = [x[0] for x in self.db_schema.PAI_OBS_CELL_COLUMNS] + [
            y[0] for y in self.db_schema.PAI_OBS_SAMPLE_COLUMNS
        ]

        self.artifact.obs = self.artifact.obs.reindex(non_index_columns, axis=1)
        for col, dtype_ in self.db_schema.PAI_OBS_CELL_COLUMNS + self.db_schema.PAI_OBS_SAMPLE_COLUMNS:
            try:
                dtype = dtype_.to_pandas_dtype()
                self.artifact.obs[col] = self.artifact.obs[col].fillna(None).astype(dtype)

            # Strings and categoricals do not have explicit pandas dtype conversions from Arrow so we have to check the exception
            except Exception as e:
                dtype = "str"
                series_ = self.artifact.obs[col].fillna("Unknown").astype(dtype)

                # Categoricals need explicit conversion
                if isinstance(dtype_, pa.DictionaryType):
                    series_ = series_.astype("category")

                self.artifact.obs[col] = series_

        for col, _ in self.db_schema.PAI_OBS_COMPUTED_COLUMNS:
            try:
                func = self.db_schema.COMPUTED_COLUMN_FUNCTIONS[col]
            except Exception as e:
                logger.warning(f"Computation function for {col} not found, skipping, {e}...")
                continue
            logger.info(f"Applying function {func.__name__} for {col}...")
            self.artifact.obs[col] = func(self.artifact)

    def _standardize_var(self):
        """
        Standardize the .var attribute of the AnnData object.
        """
        self.artifact.var = self.artifact.var.reindex([x[0] for x in self.db_schema.PAI_VAR_COLUMNS], axis=1)
        for col, dtype_ in self.db_schema.PAI_VAR_COLUMNS:
            try:
                dtype = dtype_.to_pandas_dtype()
                self.artifact.var[col] = self.artifact.var[col].fillna(None)
            except Exception as e:
                logger.warning(
                    f"Cannot standardize col {col}: {e} since it does not have associated pandas dtype. Casting to string."
                )
                dtype = "str"
                self.artifact.var[col] = self.artifact.var[col].fillna("").astype(dtype)
        self.artifact = self.artifact[:, self.artifact.var["gene"].isin(self.db_schema.SORTED_CORE_GENES)]

    def _standardize_X(self):
        """
        Standardize the .X attribute of the AnnData object.
        """
        del self.artifact.layers
        for layer_name, _ in self.db_schema.PAI_X_LAYERS:
            if layer_name == "row_raw":
                continue
            if layer_name.endswith("_raw"):
                self.artifact.layers[layer_name] = self.artifact.X.tocsr()
            elif layer_name.endswith("_norm"):
                self.artifact.layers[layer_name] = normalize_raw_array(self.artifact.X)

    def _standardize_obsm(self):
        """
        Standardize the .obsm attribute of the AnnData object.
        """
        pass

    def _generate_feature_presence(self) -> sp.coo_matrix:
        """
        Generate the feature presence of the AnnData object.
        """
        temp_var = self.artifact.var.copy()
        temp_var["presence"] = True
        merged = pd.merge(self.db_schema.VAR_DF, temp_var, on="gene", how="left")
        index = merged.index[~merged["presence"].isna()]

        presence_matrix = sp.lil_matrix(
            (self.artifact.obs["sample_name"].nunique(), len(merged)),
            dtype=self.db_schema.PAI_PRESENCE_LAYER.to_pandas_dtype(),
        )

        presence_matrix[:, index] = 1

        return presence_matrix.tocoo()

    def write(self, output_filepath: Union[str, Path, CloudPath]) -> Union[Path, CloudPath]:
        """Write the AnnDataset to H5AD format at a specific filepath

        Args:
            output_filepath (Union[str, Path, CloudPath]): Needs to end with `.h5ad`

        Returns:
            Union[Path, CloudPath]: Location of written H5AD object
        """

        # TODO: method to deal with string path logic
        if isinstance(output_filepath, str):
            if output_filepath.startswith("s3://"):
                output_filepath = CloudPath(output_filepath)
            else:
                output_filepath = Path(output_filepath)

        if output_filepath.suffix != ".h5ad":
            raise ValueError("Please ensure the output filepath is a H5AD path i.e. ends with the `.h5ad` extension.")

        logger.info(f"Saving AnnData as H5AD to {output_filepath}...")
        output_filepath.parent.mkdir(parents=True, exist_ok=True)

        if output_filepath.exists():
            logger.info(f"H5AD exists at {output_filepath.as_posix()}, overwriting...")

        if isinstance(output_filepath, Path):
            # Strings to categoricals needs to be false so we ensure data is saved in the same way that it started as
            anndata.io.write_h5ad(filepath=output_filepath, adata=self.artifact, convert_strings_to_categoricals=False)
        elif isinstance(output_filepath, CloudPath):
            temp_path = Path(f"/tmp/{output_filepath.basename}")
            # Strings to categoricals needs to be false so we ensure data is saved in the same way that it started as
            anndata.io.write_h5ad(filepath=temp_path, adata=self.artifact, convert_strings_to_categoricals=False)
            output_filepath.upload_file(temp_path, overwrite=True)

        return output_filepath
