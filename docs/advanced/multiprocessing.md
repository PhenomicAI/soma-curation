# Multiprocessing Ingestion

SOMA Curation supports parallel processing for efficient data ingestion. This is especially useful when working with large datasets that contain many samples.

## Command-Line Interface

The simplest way to use multiprocessing is through the command-line interface:

```bash
python -m soma_curation.scripts.multiprocessing_ingest \
    --processes 4 \
    --atlas-name my_atlas \
    --raw-storage-dir path/to/raw_data \
    --h5ad-storage-dir path/to/h5ads \
    --atlas-storage-dir path/to/atlas \
    --log-dir path/to/logs
```

### Key Parameters

- `--processes`: Number of worker processes to use (default: 4)
- `--atlas-name`: Name of the atlas to create or update
- `--raw-storage-dir`: Directory containing raw data
- `--h5ad-storage-dir`: Directory for temporary H5AD file storage
- `--atlas-storage-dir`: Directory for atlas storage
- `--log-dir`: Directory for logs
- `--db-schema-fp`: Optional path to custom schema file
- `--filenames-pkl`: Path to store/load successful H5AD filenames
- `--registration-mapping-pkl`: Path to store/load registration mapping

## Pipeline Steps

The multiprocessing ingestion pipeline consists of four main steps:

1. **Convert to H5AD**: Convert each study-sample pair to H5AD format in parallel
2. **Create Registration Mapping**: Map each sample to its position in the atlas (serial step)
3. **Resize Experiment**: Allocate space in the atlas (serial step)
4. **Ingest H5AD**: Load H5AD files into the atlas in parallel

## Performance Considerations

- **Process Count**: Set `--processes` to match your CPU cores (or slightly less)
- **Memory Usage**: Each process requires memory for loading a sample
- **Storage Space**: Ensure sufficient space for both raw data and H5AD files
- **Checkpointing**: Use the pickle files to resume interrupted ingestion

## Error Handling

The pipeline logs failures but continues processing other samples. You can check:

- The log files for detailed error information
- The success/failure counts in the command output
- Failed samples by comparing input data with successfully processed files

## Example Workflow

1. Organize your data following the [Data Organization](../usage/data_organization.md) guidelines
2. Run the multiprocessing ingestion script
3. Check logs for any errors or warnings
4. Access your data using the AtlasManager

```python
from soma_curation.atlas.crud import AtlasManager
from soma_curation.schema import load_schema

# Load the completed atlas
am = AtlasManager(
    atlas_name="my_atlas",
    storage_directory="path/to/atlas",
    db_schema=load_schema()
)

# Access the atlas data
with am.open(mode="r") as exp:
    # Work with the experiment data
    print(exp.ms.keys())
```
