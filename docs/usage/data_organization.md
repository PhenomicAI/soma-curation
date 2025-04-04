# Data Organization

SOMA Curation expects a specific directory structure for your raw data. Following this structure ensures that your data can be properly ingested and processed.

## Directory Structure

```
raw_data/
├── study_1/
│   ├── mtx/                 # Matrix files
│   │   ├── sample_1/        # Each sample has its own directory
│   │   │   ├── matrix.mtx.gz
│   │   │   ├── barcodes.tsv.gz
│   │   │   └── features.tsv.gz
│   │   └── sample_2/
│   │       └── ...
│   ├── cell_metadata/       # Cell-level metadata
│   │   └── study_1.tsv.gz   # Named after the study
│   └── sample_metadata/     # Sample-level metadata
│       └── study_1.tsv.gz   # Named after the study
└── study_2/
    └── ...
```

## File Formats

### Matrix Files

Matrix files use the Matrix Market format (`.mtx.gz`), which efficiently stores sparse matrices. Each sample directory must contain:

- `matrix.mtx.gz`: The gene expression matrix (genes × cells)
- `barcodes.tsv.gz`: Cell barcodes
- `features.tsv.gz`: Gene features

### Metadata Files

Metadata files use tab-separated values (`.tsv.gz`):

- `sample_metadata/<study_name>.tsv.gz`: Contains sample-level information
- `cell_metadata/<study_name>.tsv.gz`: Contains cell-level information

## Validation

SOMA Curation performs several validations when loading data:

- Checks for directory structure correctness
- Validates file formats and content
- Ensures metadata columns match the schema
- Verifies no duplicate sample names across studies

## Best Practices

1. **Consistent Naming**: Use consistent naming conventions for studies and samples
2. **Compression**: Always compress files (`.gz`) to save space
3. **Metadata Organization**: Keep metadata well-organized and consistent
4. **Schema First**: Define your schema before organizing data
5. **Sample Uniqueness**: Ensure sample names are unique across all studies
