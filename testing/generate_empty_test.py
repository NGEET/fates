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
5. Creates a template test file in the new directory.
6. If it's a functional test, generates a new Python TestCase file and updates load_functional_tests.py

"""

import os
import argparse
import textwrap
from test_generator_class import generate_test


def commandline_args():
    """Parse and return command-line arguments"""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.print_usage = parser.print_help
    subparsers = parser.add_subparsers(
        help="Two types of tests, either:", dest="test_type"
    )
    unit_parser = subparsers.add_parser("unit", help="Create an empty unit test")
    functional_parser = subparsers.add_parser(
        "functional", help="Create an empty functional test"
    )

    # Unit test parser options

    # Functional test parser options

    # common options between both subparsers
    for subparser in [unit_parser, functional_parser]:
        subparser.add_argument(
            "--test-name", default="hello_world", help="Name of the new test"
        )
        subparser.add_argument(
            "--test-sub-dir", default=None, help="Optional test subdirectory path"
        )

    # print help for both subparsers
    parser.epilog = textwrap.dedent(
        f"""\
         {unit_parser.format_help()}
         {functional_parser.format_help()}
         """
    )

    args = parser.parse_args()

    return args

def main():

    # args = commandline_args()

    # Example - create a unit test
    unit_test = generate_test("unit", "hello_world")
    unit_test.setup_test()

    # Example: Creating a functional test
    # functional_test = generate_test("functional", "hello_world")
    # functional_test.setup_test()


if __name__ == "__main__":
    main()
