# TODO: Investigate dense arrays for presence matrix # atlas/crud.py
import argparse
import tiledbsoma as soma

from ..collection import MtxCollection, H5adCollection
from ..schema import load_schema
from ..executor.executors import MultiprocessingExecutor


def main():
    parser = argparse.ArgumentParser(description="Run feature presence computation.")
    parser.add_argument("--exp-uri", type=str, default="test", help="URI of the atlas.")
    parser.add_argument("--raw-storage-dir", type=str, default="human", help="Directory to raw storage.")
    parser.add_argument("--db-schema-fp", type=str, default=None, help="Filepath for database schema.")
    parser.add_argument(
        "--raw-collection-type",
        type=str,
        choices=["mtx", "h5ad"],
        default="mtx",
        help="Type of collection to process (mtx or h5ad)",
    )
    args = parser.parse_args()

    schema = load_schema(args.db_schema_fp)

    if args.raw_collection_type == "mtx":
        collection = MtxCollection(args.raw_storage_dir)
    else:
        collection = H5adCollection(args.raw_storage_dir)

    global_var_list = schema.SORTED_CORE_GENES

    # Determine samples needed
    exp = soma.Experiment.open(args.exp_uri)
    sample_names_df = (
        exp.obs.read(column_names=["sample_name"]).concat().to_pandas().drop_duplicates().reset_index(drop=True)
    )
    sample_names_df["sample_idx"] = sample_names_df["sample_name"].cat.codes
    last_sample_idx = exp.ms["RNA"][schema.PAI_PRESENCE_MATRIX_NAME].non_empty_domain()[0][1]
    exp.close()

    samples_to_generate_df = sample_names_df[~sample_names_df["sample_idx"].isin(range(last_sample_idx))]
    samples_to_generate_df.to_csv("samples_to_generate_df.tsv.gz", sep="\t", index=False)
    # presence_matrix = collection.presence_matrix(global_var_list)


if __name__ == "__main__":
    main()
