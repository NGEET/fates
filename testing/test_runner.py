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

from fortran_test_builder import build_tests, build_exists

# local imports
from path_utils import add_cime_lib_to_path

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
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

# constants
_FILE_DIR = Path(__file__).parent.resolve()
_DEFAULT_CONFIG_FILE = _FILE_DIR / "functional_tests.cfg"
_DEFAULT_CDL_PATH = (
    _FILE_DIR.parent / "parameter_files" / "fates_params_default.cdl"
).resolve()
_CMAKE_BASE_DIR = _FILE_DIR.parent
_TEST_SUB_DIR = "testing"

@dataclass
class RunConfig:
    """Holds all configuration for the test run."""
    parameter_file: Path
    config_file: Path
    build_dir: Path
    run_dir: Path
    make_j: int
    clean: bool = False
    skip_build: bool = False
    skip_run: bool = False
    save_figs: bool = False
    verbose_make: bool = False
    test_list: str = "all"
    
    def validate(self):
        """Self-validation logic."""
        if not self.parameter_file.exists():
            raise FileNotFoundError(f"Parameter file not found: {self.parameter_file}")
        
        valid_exts = {'.cdl', '.nc'}
        if self.parameter_file.suffix not in valid_exts:
            raise ValueError(f"Parameter file must be one of {valid_exts}")

        if not self.config_file.is_file():
            raise FileNotFoundError(f"Config file not found: {self.config_file}")

        if self.skip_build and self.verbose_make:
            raise ValueError("Cannot run verbose make if skipping build.")
        
class TestDiscovery:
    """Responsible for discovering and loading test classes."""
    
    @staticmethod
    def load_tests(config_file: Path, test_filter: str) -> Dict[str, Any]:
        """Parses config and instantiates test objects."""
        full_test_dict = config_to_dict(str(config_file))
        filtered_config = parse_test_list(full_test_dict, test_filter)
        
        # Discover available test classes dynamically
        # Note: This relies on the side-effect of importing load_functional_tests
        available_subclasses = TestDiscovery._get_subclasses(FunctionalTest, FunctionalTestWithDrivers)
        
        instantiated_tests = {}
        for name, params in filtered_config.items():
            # Find the class that has the matching 'name' attribute
            match = next((cls for cls in available_subclasses if getattr(cls, 'name', None) == name), None)
            
            if match:
                instantiated_tests[name] = match(params)
            else:
                logger.warning(f"Test class for '{name}' not found. Skipping.")
                
        return instantiated_tests

    @staticmethod
    def _get_subclasses(*base_classes) -> List[Type]:
        """Recursive subclass finder."""
        subclasses = set()
        for base in base_classes:
            for sub in base.__subclasses__():
                subclasses.add(sub)
                # If you need to go deeper than one level, recurse here
        return list(subclasses)