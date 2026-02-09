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
from pathlib import Path
import argparse

import framework.utils.env_check
from framework.unit_test import UnitTest
from framework.builder import build_tests
from framework.utils.general import config_to_dict, parse_test_list

# constants for this script
_CMAKE_BASE_DIR = Path(__file__).resolve().parents[1]
_DEFAULT_CONFIG_FILE = Path(__file__).resolve().parents[0] / "config" / "unit.cfg"


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
        default=_CMAKE_BASE_DIR / "_build",
        help="Directory where tests are built.\n"
        "Will be created if it does not exist.\n",
    )

    parser.add_argument(
        "--config-file",
        type=str,
        default=_DEFAULT_CONFIG_FILE,
        help=f"Configuration file where test list is defined. Default: '{_DEFAULT_CONFIG_FILE}'",
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
        "--verbose", action="store_true", help="Run make with verbose output."
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


def main():
    """Main script
    Reads in command-line arguments and then runs the tests.
    """

    args = commandline_args()
    build_dir = Path(args.build_dir)

    full_test_dict = config_to_dict(args.config_file)
    test_dict = parse_test_list(full_test_dict, args.test_list)

    # build tests
    build_tests(
        build_dir, _CMAKE_BASE_DIR, args.make_j, clean=args.clean, verbose=args.verbose
    )

    # run unit tests
    for name, attributes in test_dict.items():
        test = UnitTest(name, attributes)
        out = test.run(build_dir)
        
        print(out)


if __name__ == "__main__":
    main()
