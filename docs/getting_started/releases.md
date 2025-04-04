# Releases and Versioning

SOMA Curation follows [Semantic Versioning](https://semver.org/) with a "v" prefix.

## Release Types

- **Stable Releases**: Format `v1.0.0`
- **Release Candidates**: Format `v1.0.0-rc`, `v1.0.0-rc1`, etc.

## Installation

### Stable Releases (PyPI)

```bash
pip install soma-curation
```

### Release Candidates (TestPyPI)

```bash
pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple soma-curation
```

## Documentation Versions

- **latest**: Most recent stable release
- **next**: Most recent release candidate
- **dev**: Latest changes from the main branch
- **vX.Y**: Major.minor version (e.g., `v1.0`)
- **vX.Y.Z**: Specific version (e.g., `v1.0.0`)
