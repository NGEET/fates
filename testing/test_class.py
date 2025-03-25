import os
import re
from abc import ABC, abstractmethod
from utils import snake_to_camel

_TEST_SUB_DIR = os.path.dirname(os.path.abspath(__file__))
_UNIT_TESTING_DIR =  "unit_testing"
_FUNCTIONAL_TESTING_DIR = "functional_testing"
_TEMPLATE_DIR = os.path.join(_TEST_SUB_DIR, 'templates')
_UNIT_CMAKE_TEMPLATE = 'cmake_utest_template.txt'
_FUNCTIONAL_CMAKE_TEMPLATE = 'cmake_ftest_template.txt'
_TEST_CMAKELISTS = os.path.join(_TEST_SUB_DIR, 'CMakeLists.txt')
_UNIT_TESTS_CONFIG = os.path.join(_TEST_SUB_DIR, 'unit_tests.cfg')
_FUNCTIONAL_TESTS_CONFIG = os.path.join(_TEST_SUB_DIR, 'functional_tests.cfg')

class TestClass(ABC):
    """Abstract base class for FATES tests"""

    def __init__(self, test_name: str, sub_dir: str=None):
        self.test_name = re.sub(r"_test$", "", test_name.lower())
        self.sub_dir = sub_dir
        self.test_dir = self.generate_test_dir()
        self.build_dir = self.generate_test_build_dir()
        self.module_name = snake_to_camel(self.test_name)
        self.file_name = self.generate_file_name()
        
        if not hasattr(self, "section_header"):
            raise ValueError("Subclass must define 'section_header'")
    
    @abstractmethod
    def generate_test_dir(self) -> str:
        """Generates a test directory path (must be implemented by subclasses)."""
        pass
    
    @abstractmethod
    def generate_test_build_dir(self) -> str:
        """Generates a test build directory path (must be implemented by subclasses)."""
        pass
    
    @abstractmethod
    def generate_file_name(self) -> str:
        """"Generates a file name for the test program (must be implemented by subclases)."""
        pass
    
    @abstractmethod
    def create_cmake_file(self):
        """Creates the new CMakeLists.txt file for the test (must be implemented by subclasses)."""
        pass
    
    @abstractmethod
    def update_config_file(self):
        """Appends information to the correct config file (must be implemented by subclasses)."""
        pass
    
    @abstractmethod
    def create_template_program(self):
        """Creates a template program file (must be implemented by subclasses)."""
        pass
        
    def create_test_directory(self):
        """Creates a new test directory

        Raises:
            FileExistsError: Directory already exists
        """
        if os.path.exists(self.test_dir):
            raise FileExistsError(f"Error: Directory '{self.test_dir}' already exists.")
        os.makedirs(self.test_dir)
        print(f"Created directory: {self.test_dir}")
    
    def update_main_cmake(self):
        """Updates the main testing/CMakeLists.txt file by adding a new test directory in
        the correct section.            

        Raises:
            FileNotFoundError: testing/CMakeLists.txt file not found
        """
        if not os.path.exists(_TEST_CMAKELISTS):
            raise FileNotFoundError(f"Warning: {_TEST_CMAKELISTS} not found.")
        
        with open(_TEST_CMAKELISTS, "r") as f:
            lines = f.readlines()
            
        new_entry = f"add_subdirectory({self.test_dir} {self.build_dir})\n"
        section_found = False 
        updated_lines = []
        
        for line in lines:
            updated_lines.append(line)
            
            # insert new entry in the correct section
            if line.strip() == self.section_header:
                section_found = True 
                
            # if we're in the correct section and hit another section header, insert it
            elif section_found and line.startswith("## "):
                updated_lines.insert(-1, new_entry) # insert before new section header
                section_found = False # stop inserting
            
        # if we never found a new section header, append at end
        if section_found:
            updated_lines.append(new_entry)
        
        # write modified contents back to file 
        with open(_TEST_CMAKELISTS, "w") as f:
            f.writelines(updated_lines)
        
        print(f"Updated {_TEST_CMAKELISTS} to include {self.test_dir}")
        
    def setup_test(self):
        """Runs all setup steps: creates directory, updates CMake, and config."""
        
        self.create_test_directory()
        self.create_cmake_file()
        self.update_main_cmake()
        self.update_config_file()
    
    @staticmethod
    def load_template(template_name: str) -> str:
        """Load a template file from the templates directory

        Args:
            template_name (str): name of template

        Raises:
            FileNotFoundError: Template file not found

        Returns:
            str: template string
        """
        
        template_path = os.path.join(_TEMPLATE_DIR, template_name)
        if not os.path.exists(template_path):
            raise FileNotFoundError(f"Error: Template '{template_path}' not found.")
        with open(template_path, "r") as f:
            return f.read()
        
class UnitTest(TestClass):
    """Class for unit tests"""
    
    section_header = "## Unit tests"
    
    def generate_test_dir(self) -> str:
        """Generates a test directory path

        Returns:
            str: test directory path
        """
        test_name = f"{self.test_name}_test"  # append "_test"
        if self.sub_dir:
            return os.path.join(_UNIT_TESTING_DIR, self.sub_dir, test_name)
        return os.path.join(_UNIT_TESTING_DIR, test_name)
    
    def generate_test_build_dir(self) -> str:
        """Generates a test build directory path

        Returns:
            str: test build directory path
        """
        return f"fates_{self.test_name}_utest"
    
    def generate_file_name(self) -> str:
        """Generates a file name for the test program.

        Returns:
            str: file name
        Raises:
            RuntimeError: self.module name not set yet
        """
        if not self.module_name:
            raise RuntimeError("self.module_name must be set before file name generated!")
        return f"test_{self.module_name}.pf"
    
    def create_cmake_file(self):
        """Creates the new CMakeLists.txt file for the test
        """
        cmake_template = self.load_template(_UNIT_CMAKE_TEMPLATE)
        cmake_path = os.path.join(self.test_dir, 'CMakeLists.txt')
        
        with open(cmake_path, 'w') as f:
            f.write(cmake_template.format(file_name=self.file_name,
                                          module_name=self.module_name))
        print(f"Created {cmake_path}")
        
    def update_config_file(self):
        
        if not os.path.exists(_UNIT_TESTS_CONFIG):
            raise FileNotFoundError(f"Warning: {_UNIT_TESTS_CONFIG} not found.")
    
        with open(_UNIT_TESTS_CONFIG, "a") as f:
            f.write(f"[{self.test_name}]\n")
            f.write(f"test_dir = {self.build_dir}\n")
        print(f"Updated {_UNIT_TESTS_CONFIG} with {_UNIT_TESTS_CONFIG}")
        
    def create_template_program(self):
        pass
                
class FunctionalTest(TestClass):
    """Class for functional tests"""
    
    section_header = "## Functional tests"
    
    def __init__(self, test_name, sub_dir = None):
        super().__init__(test_name, sub_dir)
        
        self.executable_name = f"{self.module_name}_exe"
    
    def generate_test_dir(self) -> str:
        """Generates a test directory path

        Returns:
            str: test directory path
        """
        if self.sub_dir:
            return os.path.join(_FUNCTIONAL_TESTING_DIR, self.sub_dir, self.test_name)
        return os.path.join(_FUNCTIONAL_TESTING_DIR, self.test_name)
    
    def generate_file_name(self) -> str:
        """Generates a file name for the test program.

        Returns:
            str: file name
        Raises:
            RuntimeError: self.module name not set yet
        """
        if not self.module_name:
            raise RuntimeError("self.module_name must be set before file name generated!")
        return f"Test{self.module_name}.F90"
    
    def generate_test_build_dir(self) -> str:
        """Generates a test build directory path

        Returns:
            str: test build directory path
        """
        return f"fates_{self.test_name}_ftest"

    def create_cmake_file(self):
        """Creates the new CMakeLists.txt file for the test
        """
        cmake_template = self.load_template(_FUNCTIONAL_CMAKE_TEMPLATE)
        cmake_path = os.path.join(self.test_dir, 'CMakeLists.txt')
        
        with open(cmake_path, 'w') as f:
            f.write(cmake_template.format(file_name=self.file_name,
                                          executable_name=self.executable_name))
        print(f"Created {cmake_path}")
        
    def update_config_file(self):
        
        if not os.path.exists(_FUNCTIONAL_TESTS_CONFIG):
            raise FileNotFoundError(f"Warning: {_FUNCTIONAL_TESTS_CONFIG} not found.")
    
        with open(_FUNCTIONAL_TESTS_CONFIG, "a") as f:
            f.write(f"\n[{self.test_name}]\n")
            f.write(f"test_dir = {self.build_dir}\n")
            f.write(f"test_exe = {self.executable_name}\n")
            f.write(f"out_file = None\n")
            f.write(f"use_param_file = False\n")
            f.write(f"other_args = []\n")
        print(f"Updated {_FUNCTIONAL_TESTS_CONFIG} with {self.test_name}")
                
def generate_test(test_type: str, test_name: str, sub_dir: str = None) -> TestClass:
    if test_type == "unit":
        return UnitTest(test_name, sub_dir)
    elif test_type == "functional":
        return FunctionalTest(test_name, sub_dir)
    else:
        raise RuntimeError("test_type must be one of ['unit', 'functional']")
    