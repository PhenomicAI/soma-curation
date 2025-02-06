import tiledbsoma.io

from typing import List, Literal

from ..sc_logging import logger
from ..schema import SOMA_TileDB_Context


def _create_registration_mapping(
    experiment_uri: str,
    filenames: List[str],
    measurement_name: Literal["RNA"] = "RNA",
    obs_field_name: str = "barcode",
    var_field_name: str = "gene",
    context: tiledbsoma.SOMATileDBContext = SOMA_TileDB_Context(),
) -> tiledbsoma.io.ExperimentAmbientLabelMapping:

    rm = tiledbsoma.io.register_h5ads(
        experiment_uri=experiment_uri,
        h5ad_file_names=filenames,
        measurement_name=measurement_name,
        obs_field_name=obs_field_name,
        var_field_name=var_field_name,
        context=context,
    )

    return rm
