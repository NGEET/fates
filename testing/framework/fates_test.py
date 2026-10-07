"""Abstract class for FATES tests"""

import logging
import subprocess
from pathlib import Path
from abc import ABC, abstractmethod
from framework.utils.path import get_cime_module

# setup logging
logging.basicConfig(
    level=logging.WARNING, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


class FatesTest(ABC):
    """Base class for all FATES tests"""

    def __init__(self, name: str, test_dir: str):
        self.name = name
        self.test_dir = Path(test_dir)

    @abstractmethod
    def run(self, build_dir, run_dir, param_file):
        """Each category of test will define its own requirements"""
        raise NotImplementedError

    def execute_shell(self, cmd: list[str], run_dir: Path) -> str:
        """Run executable or test

        Args:
            run_dir (Path): path to directory to run test

        Raises:
            subprocess.CalledProcessError: error from subprocess call

        Returns:
            str: output from subprocess call
        """
        cime_utils = get_cime_module("CIME.utils")

        logging.info("--> Running Test: %s", self.name)

        stat, out, _ = cime_utils.run_cmd(
            " ".join(cmd), from_dir=str(run_dir), combine_output=True
        )
        if stat:
            logging.error(out)
            raise subprocess.CalledProcessError(stat, cmd, out)
        return out
