"""
Responsible for running tests within the FATES repository
"""

import os
import sys
import shutil
import logging
from pathlib import Path
from dataclasses import dataclass, field
from typing import List, Dict, Type, Optional, Any

# local imports
from framework.builder import build_tests, build_exists
from framework.utils.path import add_cime_lib_to_path

# initialize CIME path
add_cime_lib_to_path()

try:
    from CIME.utils import run_cmd
except ImportError as e:
    # fail fast if environment isn't set up
    raise ImportError(
        f"CIME dependencies missing. Ensure environment is configured. Error: {e}"
    ) from e

# setup logging
logging.basicConfig(
    level=logging.WARNING, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

