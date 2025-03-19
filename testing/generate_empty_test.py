#!/usr/bin/env python

"""
|------------------------------------------------------------------|
|---------------------  Instructions  -----------------------------|
|------------------------------------------------------------------|

This script creates an 'empty' test (unit or functional) that can then be modified.

Specifically, this script peforms the following steps:

1. Creates a new test directory in unit_testing or functional_testing, depending on the
    script arguments
2. Generates the CMakeLists.txt file for the new test.
3. Updates the existing CMakeLists.txt file to include the new test.
4. Appends an entry for this test to the configuration file.
5. If it's a functional test, generates a new Python TestCase file.

Note that you may need to 

"""

import os 
import argparse

_TEST_SUB_DIR = os.path.dirname(os.path.abspath(__file__))
_TEMPLATE_DIR = os.path.join(_TEST_SUB_DIR, 'templates')
_UNIT_TESTING_DIR = os.path.join(_TEST_SUB_DIR, "unit_testing")
_FUNCTIONAL_TESTING_DIR = os.path.join(_TEST_SUB_DIR, "functional_testing")
_UNIT_CMAKE_TEMPLATE = 'cmake_utest_template.txt'
_TEST_CMAKELISTS = os.path.join(_TEST_SUB_DIR, 'CMakeLists.txt')

def commandline_args():
    """Parse and return command-line arguments"""

    description = """
    Driver for creating a new empty test

    Typical usage:

    ./generate_empty_test.py 

    """
    parser = argparse.ArgumentParser(
        description=description, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('test_name', help="Name of the new test")
    parser.add_argument('test_type', choices=['unit', 'functional'], help='Type of test')
    parser.add_argument('--module', default="NewTest", help="Module name for test")

    args = parser.parse_args()

    return args

def load_template(template_name):
    """Load a template file from the templates directory."""
    template_path = os.path.join(_TEMPLATE_DIR, template_name)
    if not os.path.exists(template_path):
        raise FileNotFoundError(f"Error: Template '{template_path}' not found.")
    with open(template_path, "r") as f:
        return f.read()

def generate_test_dir_name(test_type, test_name):
    if test_type == 'unit':
        return os.path.join(_UNIT_TESTING_DIR, f"fates_{test_name}_utest")
    elif test_type == 'functional':
        return os.path.join(_FUNCTIONAL_TESTING_DIR, f"fates_{test_name}_ftest")
    else:
        raise RuntimeError("test_type must be one of [unit, functional]")

def create_directory(path):
    if os.path.exists(path):
        raise FileExistsError(f"Error: Directory '{path}' already exists.")
    os.makedirs(path)
    print(f"Created directory: {path}")
    
def create_cmake_file(test_dir, test_type, module_name):
    if test_type == 'unit':
        cmake_template = load_template(_UNIT_CMAKE_TEMPLATE)
    else:
        return
        
    cmake_path = os.path.join(test_dir, "CMakeLists.txt")
    with open(cmake_path, "w") as f:
        f.write(cmake_template.format(module_name=module_name))
    print(f"Created {cmake_path}")
    
def update_main_cmake(test_dir):
    if os.path.exists(_TEST_CMAKELISTS):
        with open(_TEST_CMAKELISTS, "a") as f:
            f.write(f"\nadd_subdirectory({test_dir})\n")
        print(f"Updated {_TEST_CMAKELISTS} to include {test_dir}")
    else:
        raise FileNotFoundError(f"Warning: {_TEST_CMAKELISTS} not found.")

def main():
  
  args = commandline_args()
  
  test_dir = generate_test_dir_name(args.test_type, args.test_name)
  create_directory(test_dir)
  create_cmake_file(test_dir, args.test_type, args.module)
  update_main_cmake(test_dir)

if __name__ == "__main__":
  main()
  