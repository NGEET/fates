"""Checks python environment"""

import sys

# minimum python version
MIN_PYTHON = (3, 12)

def validate():
    """Checks python version and exits with error if not satisfied
    """
    if sys.version_info < MIN_PYTHON:
        sys.exit(
            f"[-] Failure: Python {MIN_PYTHON[0]}.{MIN_PYTHON[1]}+ is required.\n"
            f"[-] You are currently running {sys.version}."
        )
       
# run on import
validate()
