# Installation

## Prerequisites

Before installing SOMA Curation, ensure you have:

- Python 3.8 or higher
- pip (Python package installer)
- Virtual environment (recommended)

## Installation Methods

### Using pip (Stable Release)

The simplest way to install the stable release of SOMA Curation is using pip:

```bash
pip install soma-curation
```

### Using TestPyPI (Pre-releases)

To install pre-release versions for testing:

```bash
# Install from TestPyPI
pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple soma-curation
```

Pre-releases are published when a GitHub release is marked as a pre-release. These versions may contain new features that are still under development or testing.

### From Source

If you want to install from source or contribute to the development:

1. Clone the repository:

```bash
git clone https://github.com/PhenomicAI/soma-curation.git
cd soma-curation
```

2. Create and activate a virtual environment:

```bash
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

3. Install the package in development mode:

```bash
pip install -e .
```

## Dependencies

SOMA Curation has the following main dependencies:

- pandas
- numpy
- scipy
- anndata
- tiledbsoma
- cloudpathlib
- pydantic

These will be automatically installed when you install the package.

## Verifying Installation

To verify your installation, you can run:

```python
import soma_curation
print(soma_curation.__version__)
```

## Troubleshooting

If you encounter any issues during installation:

1. Ensure you're using a supported Python version
2. Try upgrading pip: `pip install --upgrade pip`
3. Check if all dependencies are properly installed
4. If using a virtual environment, ensure it's activated

For additional help, see the [Troubleshooting](examples/troubleshooting.md) guide.
