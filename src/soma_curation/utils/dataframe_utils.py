import pandas as pd

from typing import List, Optional


def add_metadata_to_df(
    df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    join: List[str],
    columns_to_add: List[str],
    fillna: Optional[str] = None,
) -> pd.DataFrame:
    """Adds cell/sample-level metadata to DataFrame obs object. All barcodes from dataframe kept - if there are fewer barcodes in the cell_metadata_df,
    those columns are filled with Unknown

    Args:
        df (pd.DataFrame): Dataframe object with obs attribute to add cell metadata to
        metadata_df (pd.DataFrame): Dataframe containing columns as barcodes, and associated columns of cell metadata
        join (List[str]): Choose what key to join the dataframe to the observation column in the adata object. Eg: `barcode` or `sample_name`
        columns_to_add (List[str]): Columns to add in AnnData observation column. These columns must be present in `cell_metadata_df`. Do not include barcode/sample_name here!
        fillna (Optional[str]): Specify if you'd like to fillna values in the merged observation dataframe with a string/value

    Raises:
        ValueError: Mismatch in shape of original df obs and merged obs if merge is incorrect
        ValueError: The merge did not preserve order of barcodes if the left join messes up the order of barcodes

    Returns:
        pd.DataFrame: DataFrame object with metadata added from dataframe
    """
    original_df_columns = list(df.columns)

    for col in columns_to_add:
        if col in df.columns and col not in join:
            df.drop(columns=col, inplace=True, axis=1)

    merged_obs = pd.merge(df, metadata_df, how="left", on=join).reindex(
        list(set(original_df_columns + columns_to_add)), axis=1
    )

    assert merged_obs.shape[0] == df.shape[0], "Mismatch in shape of original adata obs and merged obs"
    assert (merged_obs["barcode"].values == df["barcode"].values).all(), "The merge did not preserve order of barcodes"

    if fillna is not None:
        merged_obs.fillna(fillna, inplace=True)

    return merged_obs
