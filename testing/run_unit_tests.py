#!/usr/bin/env python

"""
|------------------------------------------------------------------|
|---------------------  Instructions  -----------------------------|
|------------------------------------------------------------------|

Though this script does not require any host land model code, it does require some CIME
and shr code, so you should still get these repositories as you normally would
(i.e., manage_externals, etc.)

Additionally, this requires netcdf and netcdff as well as a fortran compiler.

You must also have a .cime folder in your home directory which specifies machine
configurations for CIME.

This script builds and runs FATES units tests.

"""
import os
import argparse

from build_fortran_tests import build_tests
from path_utils import add_cime_lib_to_path
from utils import config_to_dict, parse_test_list

add_cime_lib_to_path()

from CIME.utils import run_cmd_no_fail  # pylint: disable=wrong-import-position,import-error,wrong-import-order

# constants for this script
_CMAKE_BASE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../")
_DEFAULT_CONFIG_FILE = "unit_tests.cfg"
_TEST_SUB_DIR = "testing"


def commandline_args():
    """Parse and return command-line arguments"""

    description = """
    Driver for running FATES unit tests

    Typical usage:

    ./run_fates_unit_tests -t fire_weather

    """
    parser = argparse.ArgumentParser(
        description=description, formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(
        "-b",
        "--build-dir",
        type=str,
        default=os.path.join(_CMAKE_BASE_DIR, "_build"),
        help="Directory where tests are built.\n"
        "Will be created if it does not exist.\n",
    )

    parser.add_argument(
        "--make-j",
        type=int,
        default=8,
        help="Number of processes to use for build.",
    )

    parser.add_argument(
        "-c",
        "--clean",
        action="store_true",
        help="Clean build directory before building.\n"
        "Removes CMake cache and runs 'make clean'.\n",
    )

    parser.add_argument(
        "--verbose-make", action="store_true", help="Run make with verbose output."
    )

    parser.add_argument(
        "-t",
        "--test-list",
        action="store",
        dest="test_list",
        type=str,
        default="all",
        help="Test(s) to run. Comma-separated list of test names, or 'all'\n"
        "for all tests. If not supplied, will run all tests.",
    )

    args = parser.parse_args()

    return args


def run_unit_tests(clean, verbose_make, build_dir, make_j, test_dict):
    """Builds and runs the fates unit tests

    Args:
        clean (bool): whether or not to clean the build directory
        verbose_make (bool): whether or not to run make with verbose output
        build_dir (str): build directory
        make_j (int): number of processors for the build
        test_dict (dict): dictionary of test classes to run
    """

    # absolute path to desired build directory
    build_dir_path = os.path.abspath(build_dir)

    # compile code
    build_tests(
        build_dir_path, _CMAKE_BASE_DIR, make_j, clean=clean, verbose=verbose_make
    )

    # run unit tests
    print("Running unit tests...")
    for _, attributes in test_dict.items():

        test_dir = os.path.join(build_dir_path, _TEST_SUB_DIR, attributes["test_dir"])
        ctest_command = ["ctest", "--output-on-failure"]
        output = run_cmd_no_fail(
            " ".join(ctest_command), from_dir=test_dir, combine_output=True
        )
        print(output)


def main():
    """Main script
    Reads in command-line arguments and then runs the tests.
    """

    full_test_dict = config_to_dict(_DEFAULT_CONFIG_FILE)

    args = commandline_args()
    test_dict = parse_test_list(full_test_dict, args.test_list)

    run_unit_tests(
        args.clean, args.verbose_make, args.build_dir, args.make_j, test_dict
    )


if __name__ == "__main__":
    main()
