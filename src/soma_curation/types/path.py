from typing import Union
from pydantic import AfterValidator
from typing_extensions import Annotated
from cloudpathlib import AnyPath


def expand_paths(value: Union[str, AnyPath]) -> AnyPath:
    """Convert string paths to AnyPath objects."""
    if isinstance(value, str):
        value = AnyPath(value)
    return value


ExpandedPath = Annotated[Union[str, AnyPath], AfterValidator(expand_paths)]
