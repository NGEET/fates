import logging
import subprocess
from pathlib import Path
from abc import ABC, abstractmethod
from typing import Optional
from framework.utils.path import add_cime_lib_to_path

# initialze CIME path
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

# constants
_TEST_SUB_DIR = "testing"


class FatesTest(ABC):
    """Base class for all FATES tests"""

    def __init__(self, name: str, test_dir: str, test_exe: str):
        self.name = name
        self.test_dir = Path(test_dir)
        self.test_exe = test_exe
    
    @abstractmethod
    def run(self, **kwargs):
        """Each category of test will define its own requirements"""
        raise NotImplementedError

    def find_build(self, build_dir: Path):
        """Check to see if the required binary exists and return it

        Args:
            build_dir (Path): path build directory

        Raises:
            FileNotFoundError: Could not find executable
        """
        exe_path = build_dir / _TEST_SUB_DIR / self.test_dir / self.test_exe
        if not exe_path.exists():
            raise FileNotFoundError(f"[{self.name}] Executable not found at {exe_path}")
        
        return exe_path

    def execute_shell(self, cmd: list[str], run_dir: Path) -> str:
        """Run executable or test

        Args:
            run_dir (Path): path to run directory

        Raises:
            subprocess.CalledProcessError: error from subprocess call

        Returns:
            str: output from subprocess call
        """
        
        logging.info("--> Running Test: %s", self.name)
        
        stat, out, _ = run_cmd(
            " ".join(cmd), from_dir=str(run_dir), combine_output=True
        )
        if stat:
            logging.error(out)
            raise subprocess.CalledProcessError(stat, cmd, out)
        return out
