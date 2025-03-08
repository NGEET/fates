"""Utility functions related to getting paths to various important places
"""

import os
import sys

# path to the root directory of FATES, based on the path of this file
# it's important that this NOT end with a trailing slash
_FATES_ROOT = os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir)
)


def add_cime_lib_to_path() -> str:
    """Adds the CIME python library to the python path, to allow importing
       modules from that library

    Returns:
        str: path to top-level cime directory
    """
    cime_path = path_to_cime()
    prepend_to_python_path(cime_path)

    cime_lib_path = os.path.join(cime_path, "CIME", "Tools")
    prepend_to_python_path(cime_lib_path)

    return cime_path


def path_to_cime() -> str:
    """Returns the path to cime, if it can be found

    Raises:
        RuntimeError: can't find path to cime

    Returns:
        str: full path to cime
    """
    cime_path = os.path.join(path_to_fates_root(), "../../cime")
    if os.path.isdir(cime_path):
        return cime_path
    raise RuntimeError("Cannot find cime.")


def path_to_fates_root():
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
