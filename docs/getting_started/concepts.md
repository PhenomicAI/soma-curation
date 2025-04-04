# Basic Concepts

This guide explains the fundamental concepts and components of SOMA Curation.

## Core Components

### Schema

The schema defines the structure and validation rules for your data. It includes:

- Column definitions for samples and cells
- Data types and constraints
- Version information
- Platform configurations

Example schema components:

```python
PAI_OBS_SAMPLE_COLUMNS = [
    ("sample_name", "string"),
    ("study_name", "string"),
    # ...
]

PAI_OBS_CELL_COLUMNS = [
    ("barcode", "string"),
    # ...
]
```

### Collections

Collections (`MtxCollection`) manage your raw data:

- Organize data by studies and samples
- Handle both local and cloud storage (S3)
- Validate data structure and integrity
- Provide access to metadata and matrices

### Atlas

The Atlas (`AtlasManager`) manages processed data:

- Stores data in TileDB-SOMA format
- Handles versioning and metadata
- Provides efficient data access
- Supports parallel processing

## Data Organization

### Directory Structure

```
raw_data/
├── study_1/
│   ├── mtx/              # Matrix files
│   ├── cell_metadata/    # Cell-level metadata
│   └── sample_metadata/  # Sample-level metadata
└── study_2/
    └── ...
```

### File Formats

- Matrix files: `.mtx.gz` (Matrix Market format)
- Metadata files: `.tsv.gz` (Tab-separated values)
- Features and barcodes: `.tsv.gz`

## Data Flow

1. **Raw Data**

   - Organized in studies and samples
   - Stored in specified directory structure
   - Validated against schema

2. **Collection**

   - Manages access to raw data
   - Validates data integrity
   - Provides data access methods

3. **Atlas**
   - Stores processed data
   - Manages versioning
   - Provides efficient querying

## Key Features

### Validation

- Schema compliance
- Data integrity checks
- Duplicate detection
- Format validation

### Flexibility

- Support for local and cloud storage
- Customizable schema
- Extensible architecture
- Parallel processing support

### Integration

- Compatible with TileDB-SOMA
- Works with common single-cell tools
- Supports standard file formats
- Integrates with existing pipelines

## Best Practices

1. **Schema Design**

   - Define clear column names and types
   - Include necessary metadata
   - Plan for future extensions

2. **Data Organization**

   - Follow the recommended structure
   - Use consistent naming
   - Maintain clear documentation

3. **Processing**
   - Validate data early
   - Use appropriate resources
   - Monitor progress
   - Keep logs

## Next Steps

- Learn about [Schema Definition](user_guide/schema.md)
- Explore [Data Organization](user_guide/data_organization.md)
- Check out [Advanced Usage](examples/advanced.md)
