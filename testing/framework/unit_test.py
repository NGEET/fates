from pathlib import Path
from framework.fates_test import FatesTest

class UnitTest(FatesTest):
    """Class for running FATES unit tests via CTest"""

    def __init__(self, name: str, config: dict):
        super().__init__(
            name=name, 
            test_dir=config["test_dir"]
        )
        # unit tests don't produce plots
        self.plot = False

    def run(self, build_dir: Path):
        """
        Executes a unit test. Note: Unit tests only require the build_dir
        because they run inside the build tree using ctest.
        """
        # locate the test directory in the build tree
        # unit tests typically run from: build/testing/<test_dir>
        test_path = build_dir / "testing" / self.test_dir
        
        if not test_path.exists():
            raise FileNotFoundError(f"Unit test path not found: {test_path}")

        # build the CTest command
        # --output-on-failure ensures we see the log if it crashes
        cmd = ["ctest", "--output-on-failure"]
        return self.execute_shell(cmd, test_path)