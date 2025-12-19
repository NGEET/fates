"""Class for FATES functional tests"""

from abc import abstractmethod
from pathlib import Path
from framework.fates_test import FatesTest
from framework.utils.general import str_to_bool, str_to_list, copy_file

# constants
_TEST_SUB_DIR = "testing"


class FunctionalTest(FatesTest):
    """Class for running FATES functional tests"""

    def __init__(self, name: str, config: dict):
        super().__init__(
            name=name,
            test_dir=config["test_dir"],
        )

        self.test_exe = config["test_exe"]
        self.out_file = config["out_file"]
        self.use_param_file = str_to_bool(config["use_param_file"])
        self.other_args = str_to_list(config["other_args"])
        self.plot = True

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

    def run_command(self, param_file: Path | None = None) -> list[str]:
        """Builds run command for executing binary

        Args:
            param_file (str | None, optional): input parameter file path. Defaults to None.

        Returns:
            list[str]: list of arguments
        """

        cmd = [f"./{self.test_exe}"]

        if self.use_param_file:
            cmd.append(str(param_file))

        if self.other_args:
            cmd.extend(self.other_args)

        return cmd

    def run(
        self,
        build_dir: Path,
        run_dir: Path,
        param_file: Path | None = None,
    ) -> str:
        """Runs a functional test

        Args:
            build_dir (Path): path to build directory
            run_dir (Path): path to run directory
            param_file (Path | None, optional): path to parameter file. Defaults to None.

        Returns:
            str: any output from subprocess call
        """

        # find executable
        exe_path = self.find_build(build_dir)

        # copy file to run directory
        copy_file(exe_path, run_dir)

        # get specific run command and execute
        cmd = self.run_command(param_file)
        return self.execute_shell(cmd, run_dir)

    @abstractmethod
    def plot_output(self, run_dir: str, save_figs: bool, plot_dir: str):
        """Every functional test must implement its own plotting logic"""
        raise NotImplementedError


class FunctionalTestWithDrivers(FunctionalTest):
    """Class for running FATES functional tests with driver files"""

    def __init__(self, name: str, config: dict):
        super().__init__(name=name, config=config)
        # driver tests specifically look for this key in the config
        self.datm_file = config.get("datm_file")

    def run_command(self, param_file: Path | None = None) -> list[str]:
        """Builds run command for executing binary

        Args:
            param_file (str | None, optional): input parameter file path. Defaults to None.

        Returns:
            list[str]: list of arguments
        """

        cmd = [f"./{self.test_exe}"]

        if self.use_param_file and param_file:
            cmd.append(str(param_file))

        if hasattr(self, "datm_file") and self.datm_file:
            cmd.append(str(self.datm_file))

        if self.other_args:
            cmd.extend(self.other_args)

        return cmd
