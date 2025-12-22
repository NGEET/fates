"""Module for generating empty test classes"""

import re
import logging
from pathlib import Path
from abc import ABC, abstractmethod
from framework.utils.test_layout_protocol import (
    TestLayoutProtocol,
    FunctionalLayoutProtocol,
    UnitTestLayout,
    FunctionalTestLayout,
)
from framework.utils.general import snake_to_camel

# configure logger
logger = logging.getLogger(__name__)

# root directory of this test package
_TEST_ROOT = Path(__file__).resolve().parents[1]


class GenerateTestClass(ABC):
    """Abstract base class for creating boilerplate for FATES tests"""

    def __init__(
        self,
        test_name: str,
        sub_dir: str | None = None,
        layout: TestLayoutProtocol | None = None,
        verbose: bool = False,
    ):
        self.test_name = re.sub(r"_test$", "", test_name.lower())
        self.sub_dir = sub_dir
        self.layout = layout or self.default_layout()
        self.verbose = verbose

        # values derived from test_name
        self.module_name = snake_to_camel(self.test_name)

        # set logging level
        if self.verbose:
            logging.getLogger().setLevel(logging.INFO)
        else:
            logging.getLogger().setLevel(logging.WARNING)

    # ----------------------------
    # Computed properties
    # ----------------------------

    @property
    def test_dir(self) -> Path:
        """Returns the test directory"""
        return self.layout.test_dir(self.test_name, self.sub_dir)

    @property
    def build_dir(self) -> str:
        """Returns the build name"""
        return self.layout.build_name(self.test_name)

    @property
    def file_name(self) -> str:
        """Returns the program filename"""
        return self.layout.filename(self.module_name)

    @property
    def config_file(self) -> Path:
        """Returns a config file path name"""
        return self.layout.config_file

    @property
    def cmake_template(self) -> Path:
        """Returns a cmake template file path"""
        return self.layout.cmake_template

    @property
    def program_template(self) -> Path:
        """Returns a program template file path"""
        return self.layout.program_template

    @property
    def templates_dir(self) -> Path:
        """Returns path to templates directory"""
        return self.layout.templates_dir

    @property
    def main_cmakelists(self) -> Path:
        """Returns path to main testing CMakeLists.txt"""
        return self.layout.main_cmakelists

    # ----------------------------
    # Abstract interface contract
    # ----------------------------

    @abstractmethod
    def default_layout(self) -> TestLayoutProtocol:
        """Returns the default layout for this test type"""
        raise NotImplementedError

    @property
    @abstractmethod
    def section_header(self) -> str:
        """String that marks where this test's add_subdirectory should go in the
        testing/CMakeLists.txt"""
        raise NotImplementedError

    @abstractmethod
    def append_config_lines(self):
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
        if self.test_dir.exists():
            raise FileExistsError(f"Directory '{self.test_dir}' already exists.")
        self.test_dir.mkdir(parents=True, exist_ok=True)
        if self.verbose:
            logger.info("Created directory: %d", self.test_dir)

    def update_main_cmake(self):
        """Updates the main testing/CMakeLists.txt file by adding a new test directory in
        the correct section.

        Raises:
            FileNotFoundError: testing/CMakeLists.txt file not found
        """
        if not self.main_cmakelists.exists():
            raise FileNotFoundError(f"{self.main_cmakelists} not found.")

        lines = self.main_cmakelists.read_text(encoding="utf-8").splitlines()

        test_dir_rel = self.test_dir.relative_to(_TEST_ROOT)
        new_entry = f"add_subdirectory({test_dir_rel} {self.build_dir})"

        # avoid duplicates
        if any(new_entry in line for line in lines):
            logger.warning("Test already exists in CMakeLists. Skipping.")
            return

        updated_lines = []
        found_section = False
        inserted = False

        for line in lines:
            # if we hit the next section header while in our target section,
            # or if we are at the very end of the file and haven't inserted yet
            is_next_header = line.startswith("##") and found_section

            if is_next_header and not inserted:
                # if there's a blank line before the next header,
                # we want to insert BEFORE that blank line.
                if len(updated_lines) > 0 and updated_lines[-1].strip() == "":
                    updated_lines.insert(-1, new_entry.strip())
                else:
                    updated_lines.append(new_entry.strip())
                inserted = True

            updated_lines.append(line)

            if line.strip() == self.section_header:
                found_section = True

        # if we reached the end of the file and haven't inserted (last section)
        if found_section and not inserted:
            # check for trailing empty lines to keep it clean
            while updated_lines and updated_lines[-1].strip() == "":
                updated_lines.pop()
            updated_lines.append(new_entry.strip())

        # write modified contents back to file
        self.main_cmakelists.write_text("\n".join(updated_lines) + "\n")
        logger.info("Updated %s to include %s.", self.main_cmakelists, self.test_dir)

    def update_config_file(self):
        """Loads config file lines"""

        if not self.config_file.exists():
            raise FileNotFoundError(f"{self.config_file} not found.")

        content = self.config_file.read_text(encoding="utf-8").rstrip()

        # get the new lines from the subclass
        new_lines = self.append_config_lines()

        # add a double new line before the new block
        new_block = "\n" + "".join(new_lines)

        self.config_file.write_text(content + new_block + "\n")
        logger.info("Updated %s.", self.config_file)

    def load_template(self, template_path: Path) -> str:
        """Load a template file from the templates directory

        Args:
            template_path (str): name of template

        Raises:
            FileNotFoundError: Template file not found

        Returns:
            Path: template path
        """
        if not template_path.exists():
            raise FileNotFoundError(f"Template '{template_path}' not found.")
        return template_path.read_text(encoding="utf-8")

    def render_template(self, template_path: Path, out_path: Path, **kwargs):
        """Write text to a file based on an input template and keyword arguments

        Args:
            template_path (Path): path of template
            out_path (Path): path to file to write
        """
        template = self.load_template(template_path)
        out_path.write_text(template.format(**kwargs), encoding="utf-8")

    # --------------------------
    # Template Method
    # --------------------------

    def setup_test(self):
        """Template method for setting up a test"""
        self.create_test_directory()
        self.create_cmake_file()
        self.update_main_cmake()
        self.update_config_file()
        self.create_template_program()


# ---------------------------------------------------------
# Concrete Unit Test Generator
# ---------------------------------------------------------
class GenerateUnitTest(GenerateTestClass):
    """Concrete generator for unit test boilerplate"""

    layout: TestLayoutProtocol

    def default_layout(self) -> TestLayoutProtocol:
        """Sets the default layout for file naming"""
        return UnitTestLayout(_TEST_ROOT)

    @property
    def section_header(self) -> str:
        """Section header for the CMakeLists.txt file"""
        return "## Unit tests"

    def create_cmake_file(self):
        """Create the test's CMakeLists.txt."""
        cmake_path = self.test_dir / "CMakeLists.txt"
        self.render_template(
            self.cmake_template,
            cmake_path,
            file_name=self.file_name,
            module_name=self.module_name,
        )
        logger.info("Created %s.", cmake_path)

    def append_config_lines(self) -> list[str]:
        """Appends this test's config section to the config file contents

        Args:
            lines (list[str]): lines to append to

        Returns:
            list[str]: updated lines
        """
        lines = ["\n", f"[{self.test_name}]\n", f"test_dir = {self.build_dir}\n"]
        return lines

    def create_template_program(self):
        """Generate the .pf unit test boilerplate file."""
        pfunit_path = self.test_dir / self.file_name
        self.render_template(
            self.program_template, pfunit_path, module_name=self.module_name
        )
        logger.info("Added template test files in %s.", self.test_dir)


# ---------------------------------------------------------
# Concrete Functional Test Generator
# ---------------------------------------------------------
class GenerateFunctionalTest(GenerateTestClass):
    """Concrete generator for functional test boilerplate"""

    layout: FunctionalLayoutProtocol

    def default_layout(self) -> FunctionalLayoutProtocol:
        return FunctionalTestLayout(_TEST_ROOT)

    @property
    def section_header(self) -> str:
        return "## Functional tests"

    @property
    def class_template(self) -> Path:
        """Returns a class template path"""
        return self.layout.class_template

    @property
    def executable_name(self) -> str:
        """Defines the name of the executable"""
        return self.layout.executable_name(self.module_name)

    def create_cmake_file(self):
        """Creates the test's CMakeLists.txt file for the test"""
        cmake_path = self.test_dir / "CMakeLists.txt"
        self.render_template(
            self.cmake_template,
            cmake_path,
            file_name=self.file_name,
            executable_name=self.executable_name,
        )
        logger.info("Created %s.", cmake_path)

    def append_config_lines(self) -> list[str]:
        """Append this functional test's config section to the config file

        Args:
            lines (list[str]): lines to append to

        Returns:
            str: updated lines
        """
        lines = [
            "\n",
            f"[{self.test_name}]\n",
            "test_dir = {self.build_dir}\n",
            f"test_exe = {self.executable_name}\n",
            "out_file = None\n",
            "use_param_file = False\n",
            "other_args = []\n",
        ]
        return lines

    def create_template_program(self):
        """Generate the functional test Fortran program and Python class.

        Raises:
            FileNotFoundError: Can't find the class file
        """

        # create a fortran template file
        test_path = self.test_dir / self.file_name
        self.render_template(
            self.program_template, test_path, module_name=self.module_name
        )

        # create the python template file
        class_path = self.test_dir / f"{self.test_name}_test.py"
        self.render_template(
            self.class_template, class_path, module_name=self.module_name
        )
        logger.info("Added template test files in %s.", self.test_dir)
        
        # create the __init__.py
        init_path = self.test_dir / "__init__.py"
        init_path.touch()
        logger.info("Added __init__.py file in %s.", self.test_dir)


# ---------------------------------------------------------
# Factory
# ---------------------------------------------------------
TEST_TYPE_REGISTRY = {"unit": GenerateUnitTest, "functional": GenerateFunctionalTest}


def generate_test(
    test_type: str, test_name: str, sub_dir: str | None = None, verbose: bool = True
) -> GenerateTestClass:
    """Generates all boilerplate associated with a new test

    Args:
        test_type (str): test type ('unit' or 'functional')
        test_name (str): name of test
        sub_dir (str, optional): optional subdirectory to place. Defaults to None.
        verbose (bool, optional): whether or not to print log messages to screen

    Raises:
        ValueError: incorrect test type

    Returns:
        GenerateTestClass: test generator instance
    """
    cls = TEST_TYPE_REGISTRY.get(test_type)
    if cls is None:
        raise ValueError(
            f"Invalid test_type '{test_type}'. Must be one of {list(TEST_TYPE_REGISTRY.keys())}"
        )
    return cls(test_name, sub_dir, layout=None, verbose=verbose)
