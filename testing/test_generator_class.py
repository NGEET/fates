"""Module for generating empty test classes"""

import re
from pathlib import Path
from abc import ABC, abstractmethod
from file_system import FileSystem
from utils import snake_to_camel

# root directory of this test package
_TEST_ROOT = Path(__file__).resolve().parent

# top-level test directories
_UNIT_TESTING_DIR = _TEST_ROOT / "unit_testing"
_FUNCTIONAL_TESTING_DIR = _TEST_ROOT / "functional_testing"

# template directory
_TEMPLATE_DIR = _TEST_ROOT / "templates"

# template file names
_UNIT_CMAKE_TEMPLATE = "cmake_utest_template.txt"
_FUNCTIONAL_CMAKE_TEMPLATE = "cmake_ftest_template.txt"
_PFUNIT_TEMPLATE = "pfunit_template.txt"
_FUNCTIONAL_TEST_TEMPLATE = "fortran_test_template.txt"
_CLASS_TEMPLATE = "test_class_template.txt"

# configuration files
_TEST_CMAKELISTS = _TEST_ROOT / "CMakeLists.txt"
_UNIT_TESTS_CONFIG = _TEST_ROOT / "unit_tests.cfg"
_FUNCTIONAL_TESTS_CONFIG = _TEST_ROOT / "functional_tests.cfg"

# dynamic import loader
_LOAD_CLASS_FILE = _TEST_ROOT / "load_functional_tests.py"


class GenerateTestClass(ABC):
    """Abstract base class for creating empty FATES tests"""

    def __init__(
        self, test_name: str, sub_dir: str | None = None, fs: FileSystem | None = None
    ):
        # normalize test name
        self.test_name = re.sub(r"_test$", "", test_name.lower())
        self.sub_dir = sub_dir
        self.fs = fs or FileSystem()

        # values derived from test_name
        self.test_dir = self.generate_test_dir()
        self.build_dir = self.generate_test_build_dir()
        self.module_name = snake_to_camel(self.test_name)
        self.file_name = self.generate_file_name()

    # ----------------------------
    # Abstract interface contract
    # ----------------------------

    @property
    @abstractmethod
    def section_header(self) -> str:
        """String that marks where this test's add_subdirectory should go in the
        testing/CMakeLists.txt"""
        raise NotImplementedError

    @property
    @abstractmethod
    def config_file(self) -> Path:
        """Path to the config file that should be modified during test creation."""
        raise NotImplementedError

    @abstractmethod
    def generate_test_dir(self) -> Path:
        """Generates a test directory path (must be implemented by subclasses)."""
        raise NotImplementedError

    @abstractmethod
    def generate_test_build_dir(self) -> str:
        """Generates a test build directory name (must be implemented by subclasses)."""
        raise NotImplementedError

    @abstractmethod
    def generate_file_name(self) -> str:
        """ "Generates a file name for the test program (must be implemented by
        subclases)."""
        raise NotImplementedError

    @abstractmethod
    def append_config_lines(self, lines):
        """Appends lines to the config file for the test (must be implemented by
        subclasses)."""
        raise NotImplementedError

    @abstractmethod
    def create_cmake_file(self):
        """Creates the new CMakeLists.txt file for the test (must be implemented by
        subclasses)."""
        raise NotImplementedError

    @abstractmethod
    def create_template_program(self):
        """Creates a template program file (must be implemented by subclasses)."""
        raise NotImplementedError

    # --------------------------
    # Concrete helpers
    # --------------------------

    def create_test_directory(self):
        """Creates a new test directory

        Raises:
            FileExistsError: Directory already exists
        """
        if self.fs.exists(self.test_dir):
            raise FileExistsError(f"Error: Directory '{self.test_dir}' already exists.")
        self.fs.makedirs(self.test_dir)
        print(f"Created directory: {self.test_dir}")

    def update_main_cmake(self):
        """Updates the main testing/CMakeLists.txt file by adding a new test directory in
        the correct section.

        Raises:
            FileNotFoundError: testing/CMakeLists.txt file not found
        """
        if not self.fs.exists(_TEST_CMAKELISTS):
            raise FileNotFoundError(f"ERROR: {_TEST_CMAKELISTS} not found.")
        lines = self.fs.read_lines(_TEST_CMAKELISTS)

        new_entry = f"add_subdirectory({self.test_dir} {self.build_dir})\n"
        updated_lines = []
        inside_block = False

        for line in lines:
            updated_lines.append(line)

            # insert new entry in the correct section
            if line.strip() == self.section_header:
                inside_block = True
                continue

            # if we're in the correct section and hit another section header, insert it
            if inside_block and line.startswith("## "):
                updated_lines.insert(-2, new_entry)  # insert before new section header
                inside_block = False  # stop inserting

        # if we never found a new section header, append at end
        if inside_block:
            updated_lines.append(new_entry)

        # write modified contents back to file
        self.fs.write_lines(_TEST_CMAKELISTS, updated_lines)
        print(f"Updated {_TEST_CMAKELISTS} to include {self.test_dir}.")

    def update_config_file(self):
        """Loads config file lines"""

        if not self.fs.exists(self.config_file):
            raise FileNotFoundError(f"ERROR: {self.config_file} not found.")
        lines = self.fs.read_lines(self.config_file)

        # strip trailing newlines and whitespace-only lines
        while lines and lines[-1].strip() == "":
            lines.pop()

        # add one blank line before new section
        lines.append("\n")

        # new config block
        lines = self.append_config_lines(lines)

        self.fs.write_lines(self.config_file, lines)
        print(f"Updated {self.config_file}.")

    def load_template(self, template_name: str) -> str:
        """Load a template file from the templates directory

        Args:
            template_name (str): name of template

        Raises:
            FileNotFoundError: Template file not found

        Returns:
            str: template string
        """

        template_path = _TEMPLATE_DIR / template_name
        if not self.fs.exists(template_path):
            raise FileNotFoundError(f"ERROR: Template '{template_path}' not found.")
        return self.fs.read_text(template_path)
    
    def render_template(self, template_name: str, out_path: Path, **kwargs):
        """Write text to a file based on an input template and keyword arguments

        Args:
            template_name (str): name of template
            out_path (Path): path to file to write
        """
        template = self.load_template(template_name)
        rendered = template.format(**kwargs)
        self.fs.write_file(out_path, rendered)

    def setup_test(self):
        """Runs all setup steps: creates directory, updates CMake, and config."""

        self.create_test_directory()
        self.create_cmake_file()
        self.update_main_cmake()
        self.update_config_file()
        self.create_template_program()


class GenerateUnitTest(GenerateTestClass):
    """Concrete generator for unit test boilerplate"""

    @property
    def section_header(self) -> str:
        return "## Unit tests"

    @property
    def config_file(self) -> Path:
        return _UNIT_TESTS_CONFIG

    # ---------------------------------------------------------
    # Directory + filename generation
    # ---------------------------------------------------------

    def generate_test_dir(self) -> Path:
        """Returns the directory path where this unit test will be created

        Returns:
            Path: test directory path
        """
        test_name = f"{self.test_name}_test"  # append "_test"
        if self.sub_dir:
            return _UNIT_TESTING_DIR / self.sub_dir / test_name
        return _UNIT_TESTING_DIR / test_name

    def generate_test_build_dir(self) -> str:
        """Return the build-directory name

        Returns:
            str: test build directory path
        """
        return f"fates_{self.test_name}_utest"

    def generate_file_name(self) -> str:
        """Returns the filename for the test program.

        Returns:
            str: file name
        Raises:
            RuntimeError: self.module name not set yet
        """
        if not self.module_name:
            raise RuntimeError(
                "self.module_name must be set before file name generated!"
            )
        return f"test_{self.module_name}.pf"

    # ---------------------------------------------------------
    # Files that get written
    # ---------------------------------------------------------

    def create_cmake_file(self):
        """Create the test's CMakeLists.txt."""
        cmake_path = self.test_dir / "CMakeLists.txt"
        self.render_template(_UNIT_CMAKE_TEMPLATE, cmake_path, file_name=self.file_name, 
                             module_name=self.module_name)
        print(f"Created {cmake_path}.")

    def append_config_lines(self, lines: list[str]) -> list[str]:
        """Appends this test's config section to the config file contents

        Args:
            lines (list[str]): lines to append to

        Returns:
            list[str]: updated lines
        """

        lines.append(f"[{self.test_name}]\n")
        lines.append(f"test_dir = {self.build_dir}\n")
        return lines

    def create_template_program(self):
        """Generate the .pf unit test boilerplate file."""
        pfunit_path = self.test_dir / self.file_name
        self.render_template(_PFUNIT_TEMPLATE, pfunit_path, module_name=self.module_name)
        print(f"Added template test files in {self.test_dir}.")

class GenerateFunctionalTest(GenerateTestClass):
    """Concrete generator for functional test boilerplate"""

    @property
    def section_header(self) -> str:
        return "## Functional tests"

    @property
    def config_file(self) -> Path:
        return _FUNCTIONAL_TESTS_CONFIG

    def __init__(self, test_name: str, sub_dir: str | None = None):
        super().__init__(test_name, sub_dir)
        self.executable_name = f"{self.module_name}_exe"

    # ---------------------------------------------------------
    # Directory + filename generation
    # ---------------------------------------------------------

    def generate_test_dir(self) -> Path:
        """Return the functional test directory path.

        Returns:
            Path: test directory path
        """
        if self.sub_dir:
            return _FUNCTIONAL_TESTING_DIR / self.sub_dir / self.test_name
        return _FUNCTIONAL_TESTING_DIR / self.test_name

    def generate_file_name(self) -> str:
        """Return the filename for the functional test's Fortran program.

        Returns:
            str: file name
        Raises:
            RuntimeError: self.module name not set yet
        """
        if not self.module_name:
            raise RuntimeError(
                "self.module_name must be set before file name generated!"
            )
        return f"Test{self.module_name}.F90"

    def generate_test_build_dir(self) -> str:
        """Return the build-directory name used in add_subdirectory().

        Returns:
            str: test build directory path
        """
        return f"fates_{self.test_name}_ftest"

    # ---------------------------------------------------------
    # Files that get written
    # ---------------------------------------------------------

    def create_cmake_file(self):
        """Creates the test's CMakeLists.txt file for the test"""
        cmake_path = self.test_dir / "CMakeLists.txt"
        self.render_template(_FUNCTIONAL_CMAKE_TEMPLATE, cmake_path,
                             file_name=self.file_name,
                             executable_name=self.executable_name)
        print(f"Created {cmake_path}.")

    def append_config_lines(self, lines: list[str]) -> list[str]:
        """Append this functional test's config section to the config file

        Args:
            lines (list[str]): lines to append to

        Returns:
            str: updated lines
        """
        lines.append(f"[{self.test_name}]\n")
        lines.append(f"test_dir = {self.build_dir}\n")
        lines.append(f"test_exe = {self.executable_name}\n")
        lines.append("out_file = None\n")
        lines.append("use_param_file = False\n")
        lines.append("other_args = []\n")
        return lines

    def create_template_program(self):
        """Generate the functional test Fortran program and Python class.

        Raises:
            FileNotFoundError: Can't find the class file
        """

        # create a fortran template file
        test_path = self.test_dir / self.file_name
        self.render_template(_FUNCTIONAL_TEST_TEMPLATE, test_path,
                             module_name=self.module_name)
        
        # create the python template file
        class_path = self.test_dir / f"{self.test_name}.py"
        self.render_template(_CLASS_TEMPLATE, class_path, module_name=self.module_name)
        
        print(f"Added template test files in {self.test_dir}.")

        # add import to class loader file
        if not self.fs.exists(_LOAD_CLASS_FILE):
            raise FileNotFoundError(f"ERROR: {_LOAD_CLASS_FILE} not found.")

        if self.sub_dir:
            name = f"{self.sub_dir}.{self.test_name}"
        else:
            name = self.test_name

        content = f"from functional_testing.{name} import {self.module_name}\n"
        self.fs.append_file(_LOAD_CLASS_FILE, content)
        print(f"Updated {_LOAD_CLASS_FILE}.")


def generate_test(
    test_type: str, test_name: str, sub_dir: str | None = None
) -> GenerateTestClass:
    """Generates all boilerplate associated with a new test

    Args:
        test_type (str): test type ('unit' or 'functional')
        test_name (str): name of test
        sub_dir (str, optional): optional subdirectory to place. Defaults to None.

    Raises:
        RuntimeError: _description_

    Returns:
        GenerateTestClass: _description_
    """
    if test_type == "unit":
        return GenerateUnitTest(test_name, sub_dir)
    if test_type == "functional":
        return GenerateFunctionalTest(test_name, sub_dir)
    raise RuntimeError("test_type must be one of ['unit', 'functional']")
