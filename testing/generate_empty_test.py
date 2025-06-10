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

"""

import os 
import argparse
import textwrap
from utils import snake_to_camel
from test_class import generate_test

def commandline_args():
    """Parse and return command-line arguments"""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.print_usage = parser.print_help
    subparsers = parser.add_subparsers(
        help="Two types of tests, either:", dest="test_type"
    )
    unit_parser = subparsers.add_parser('unit', help="Create an empty unit test")
    functional_parser = subparsers.add_parser('functional', help="Create an empty functional test")
    
    # Unit test parser options
    
    # Functional test parser options
    
    # common options between both subparsers
    for subparser in [unit_parser, functional_parser]:
        subparser.add_argument('--test-name', default="hello_world", help="Name of the new test")
        subparser.add_argument('--test-sub-dir', default=None, help='Optional test subdirectory path')
        
    # print help for both subparsers
    parser.epilog = textwrap.dedent(
        f"""\
         {unit_parser.format_help()}
         {functional_parser.format_help()}
         """
    )
    
    args = parser.parse_args()
    
    return args

def get_config_path(test_type: str) -> str:
    """Returns path to config file based on input test type

    Args:
        test_type (str): test type [unit, functional]

    Raises:
        RuntimeError: test_type must be either unit or functional

    Returns:
        str: path to config file
    """
    if test_type == 'unit':
        return _UNIT_TESTS_CONFIG
    elif test_type == 'functional':
        return _FUNCTIONAL_TESTS_CONFIG
    else:
        raise RuntimeError("test_type must be one of [unit, functional]")
        
def update_config(config_path: str, test_name: str, build_dir: str):
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Warning: {config_path} not found.")
    
    if os.path.exists(config_path):
        with open(config_path, "a") as f:
            f.write(f"[{test_name}]\n")
            f.write(f"test_dir = {build_dir}\n")
        print(f"Updated {config_path} with {test_name}")

        

def main():
  
    #args = commandline_args()

    # Example - create a unit test
    unit_test = generate_test("unit", "hello_world")
    unit_test.setup_test()

    # Example: Creating a functional test
    functional_test = generate_test("functional", "hello_world")
    functional_test.setup_test()

if __name__ == "__main__":
  main()
  