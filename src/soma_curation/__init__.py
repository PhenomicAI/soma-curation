from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("soma-curation")
except PackageNotFoundError:
    # package is not installed
    pass
