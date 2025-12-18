"""Utility functions related to getting paths to various important places"""

import sys
import importlib
from pathlib import Path

# path to the root directory of FATES, based on the path of this file
_FATES_ROOT = Path(__file__).resolve().parents[3]


def add_cime_lib_to_path():
    """Adds the CIME python library to the python path, to allow importing
       modules from that library

    Returns:
        str: path to top-level cime directory
    """
    cime_path = path_to_cime()
    prepend_to_python_path(str(cime_path))

    cime_lib_path = (cime_path / "CIME" / "Tools").resolve()
    prepend_to_python_path(str(cime_lib_path))

    try:
        import CIME
    except ImportError as e:
        raise ImportError(
            f"Could not find CIME at {cime_lib_path}. "
            f"Ensure git_fleximod has been run. Error: {e}"
        ) from e


def path_to_cime() -> Path:
    """Returns the path to cime, if it can be found

    Raises:
        RuntimeError: can't find path to cime

    Returns:
        str: full path to cime
    """
    cime_path = (path_to_fates_root() / "../../cime").resolve()
    if cime_path.is_dir():
        return cime_path
    raise RuntimeError("Cannot find cime.")


def path_to_fates_root() -> Path:
    """Returns Returns the path to the root directory of FATES

    Returns:
        str: path to the root directory of FATES
    """
    return _FATES_ROOT


def prepend_to_python_path(path: str):
    """Adds the given path to python's sys.path if not already there
       Path is added near the beginning, so that it takes precedence over existing
       entries in path.

    Args:
        path (str): input path
    """
    if not path in sys.path:
        # insert at location 1 rather than 0, because 0 is special
        sys.path.insert(1, path)


def get_cime_module(module_path: str):
    """Safely retrieves a CIME module.
    Usage: utils = get_cime_module('CIME.utils')

    Args:
        module_path (str): module path string

    Raises:
        ImportError: Failed to load CIME module

    Returns:
        ModuleType: module
    """

    add_cime_lib_to_path()
    try:
        return importlib.import_module(module_path)
    except ImportError as e:
        raise ImportError(f"Failed to load {module_path} from CIME: {e}") from e
