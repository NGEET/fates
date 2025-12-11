"""Classes to define test layout and directory policy patterns for FATES tests"""

from pathlib import Path
from abc import ABC, abstractmethod


class TestLayoutProtocol(ABC):
    """Abstract interface for test layout classes."""

    def __init__(self, root: Path):
        self.root = root

    @property
    def templates_dir(self) -> Path:
        """Defines the templates directory"""
        return self.root / "templates"

    @property
    def main_cmakelists(self) -> Path:
        """Defines the main testing/CMakeLists.txt directory"""
        return self.root / "CMakeLists.txt"

    @property
    @abstractmethod
    def base_dir(self) -> Path:
        """Base directory for this type of test"""
        raise NotImplementedError

    @property
    @abstractmethod
    def config_file(self) -> Path:
        """Defines the config file path for a type of test"""
        raise NotImplementedError

    @property
    @abstractmethod
    def program_template(self) -> Path:
        """Defines the program file template file path"""
        raise NotImplementedError

    @property
    @abstractmethod
    def cmake_template(self) -> Path:
        """Defines the CMakeLists.txt template file path"""
        raise NotImplementedError

    @abstractmethod
    def test_dir(self, test_name: str, sub_dir: str | None = None) -> Path:
        """Defines a test directory path"""
        raise NotImplementedError

    @abstractmethod
    def build_name(self, test_name: str) -> str:
        """Defines a build directory name"""
        raise NotImplementedError

    @abstractmethod
    def filename(self, module_name: str) -> str:
        """Defines a file name"""
        raise NotImplementedError


class FunctionalLayoutProtocol(TestLayoutProtocol, ABC):
    """Extra properties for functional test layouts."""

    @property
    @abstractmethod
    def class_template(self) -> Path:
        """Defines the class template"""
        raise NotImplementedError

    @property
    @abstractmethod
    def load_class_file(self) -> Path:
        """Defines the load class file path"""
        raise NotImplementedError

    @abstractmethod
    def executable_name(self, module_name) -> str:
        """Defines the name of the executable"""
        raise NotImplementedError


class UnitTestLayout(TestLayoutProtocol):
    """Class for defining directories and filenames for unit tests"""

    @property
    def base_dir(self) -> Path:
        return self.root / "unit_testing"

    @property
    def config_file(self) -> Path:
        return self.root / "unit_tests.cfg"

    @property
    def cmake_template(self) -> Path:
        return self.templates_dir / "cmake_utest_template.txt"

    @property
    def program_template(self) -> Path:
        return self.templates_dir / "pfunit_template.txt"

    def test_dir(self, test_name: str, sub_dir: str | None = None) -> Path:
        """Defines a test directory path

        Args:
            test_name (str): name of test
            sub_dir (str | None, optional): subdirectory name. Defaults to None.

        Returns:
            Path: full path to directory
        """
        name = f"{test_name}_test"
        return self.base_dir / sub_dir / name if sub_dir else self.base_dir / name

    def build_name(self, test_name: str) -> str:
        """Defines a build directory name

        Args:
            test_name (str): name of test

        Returns:
            str: build directory name
        """
        return f"fates_{test_name}_utest"

    def filename(self, module_name: str) -> str:
        """Defines a test file name

        Args:
            module_name (str): name of module

        Returns:
            str: file name
        """
        return f"test_{module_name}.pf"


class FunctionalTestLayout(FunctionalLayoutProtocol):
    """Class for defining directories and filenames for unit tests"""

    @property
    def base_dir(self) -> Path:
        return self.root / "functional_testing"

    @property
    def config_file(self) -> Path:
        return self.root / "functional_tests.cfg"

    @property
    def cmake_template(self) -> Path:
        return self.templates_dir / "cmake_ftest_template.txt"

    @property
    def program_template(self) -> Path:
        return self.templates_dir / "fortran_test_template.txt"

    @property
    def class_template(self) -> Path:
        return self.templates_dir / "test_class_template.txt"

    @property
    def load_class_file(self) -> Path:
        return self.root / "load_functional_tests.py"

    def test_dir(self, test_name: str, sub_dir: str | None = None) -> Path:
        """Defines a test directory path

        Args:
            test_name (str): name of test
            sub_dir (str | None, optional): subdirectory name. Defaults to None.

        Returns:
            Path: full path to directory
        """
        return (
            self.base_dir / sub_dir / test_name
            if sub_dir
            else self.base_dir / test_name
        )

    def build_name(self, test_name: str) -> str:
        """Defines a build directory name

        Args:
            test_name (str): name of test

        Returns:
            str: build directory name
        """
        return f"fates_{test_name}_ftest"

    def filename(self, module_name: str) -> str:
        """Defines a test file name

        Args:
            module_name (str): name of module

        Returns:
            str: file name
        """
        return f"test_{module_name}.F90"

    def executable_name(self, module_name: str) -> str:
        """Defines the name of the executable

        Args:
            module_name (str): name of module

        Returns:
            str: file name
        """
        return f"{module_name}_exe"
