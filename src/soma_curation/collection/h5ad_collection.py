from pydantic import BaseModel, ConfigDict
import pandas as pd
import anndata as ad

from typing import List, Generator, Optional
from cloudpathlib import AnyPath

from ..sc_logging import logger
from ..types.path import ExpandedPath


class H5adCollection(BaseModel):
    """A collection manager for H5AD files in a directory

    This class manages a directory containing H5AD files and provides methods to
    list and access them. It's designed for pipelines that need to process
    multiple H5AD files in a directory.

    Args:
        storage_directory (ExpandedPath): Path to the directory containing H5AD files
        include (Optional[List[str]]): A list of filenames to include. Useful when you want
                                      to process only specific files

    Returns:
        H5adCollection: Collection manager for H5AD files
    """

    storage_directory: ExpandedPath
    include: Optional[List[str]] = None

    model_config = ConfigDict(arbitrary_types_allowed=True)

    def list_h5ad_files(self) -> List[str]:
        """List all H5AD files in the storage directory.

        Returns:
            List[str]: List of H5AD filenames in the storage directory
        """
        all_files = [
            file_path.parts[-1]
            for file_path in self.clean_listdir(self.storage_directory)
            if file_path.parts[-1].endswith(".h5ad")
        ]

        if self.include:
            filtered_files = [file for file in all_files if file in self.include]
            logger.info(f"Filtered files using include parameter: {filtered_files}")
            return filtered_files
        return all_files

    def get_h5ad_path(self, filename: str) -> AnyPath:
        """Get the full path to an H5AD file.

        Args:
            filename (str): The name of the H5AD file

        Returns:
            AnyPath: Full path to the H5AD file
        """
        return self.storage_directory / filename

    def get_anndata(self, filename: str) -> ad.AnnData:
        """Read an H5AD file into an AnnData object.

        Args:
            filename (str): The name of the H5AD file to read

        Returns:
            ad.AnnData: AnnData object containing the data from the H5AD file
        """
        logger.info(f"Reading H5AD file: {filename} from {self.storage_directory}...")
        file_path = self.get_h5ad_path(filename)

        try:
            adata = ad.read_h5ad(file_path)
            return adata
        except Exception as e:
            logger.error(f"Error reading H5AD file {filename}: {e}")
            raise

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
