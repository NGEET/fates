"""Utility functions related to getting paths to various important places
"""

import os
import sys

# Path to the root directory of FATES, based on the path of this file
# Note: It's important that this NOT end with a trailing slash;
_FATES_ROOT = os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir)
)

def path_to_fates_root():
    """Returns the path to the root directory of FATES"""
    return _FATES_ROOT

def path_to_cime():
    """Returns the path to cime, if it can be found

    Raises a RuntimeError if it cannot be found

    """
    cime_path = os.path.join(path_to_fates_root(), "../../cime")
    if os.path.isdir(cime_path):
        return cime_path
    raise RuntimeError("Cannot find cime.")

def prepend_to_python_path(path):
    """Adds the given path to python's sys.path if it isn't already in the path

    The path is added near the beginning, so that it takes precedence over existing
    entries in the path
    """
    if not path in sys.path:
        # Insert at location 1 rather than 0, because 0 is special
        sys.path.insert(1, path)

def add_cime_lib_to_path():
    """Adds the CIME python library to the python path, to allow importing
    modules from that library

    Returns the path to the top-level cime directory

    For documentation on standalone_only: See documentation in
    path_to_cime
    """
    cime_path = path_to_cime()
    prepend_to_python_path(cime_path)
    cime_lib_path = os.path.join(cime_path, "CIME", "Tools")
    prepend_to_python_path(cime_lib_path)
    return cime_path
