#!/usr/bin/env python

"""
|------------------------------------------------------------------|
|---------------------  Instructions  -----------------------------|
|------------------------------------------------------------------|
To run this script the following python packages are required:
        - numpy
        - xarray
        - matplotlib
        - pandas

Though this script does not require any host land model code, it does require some CIME
and shr code, so you should still get these repositories as you normally would
(i.e., manage_externals, etc.)

Additionally, this requires netcdf and netcdff as well as a fortran compiler.

You must also have a .cime folder in your home directory which specifies machine
configurations for CIME.

This script builds and runs FATES functional tests, and plots any relevant output from
those tests.

You can supply your own parameter file (a json formatted file), or if you do not
specify anything, the script will use the default FATES parameter json file.

"""
import os
import sys
import argparse
import subprocess
import matplotlib.pyplot as plt

from framework.builder import build_tests, build_exists
from framework.utils.path import add_cime_lib_to_path
from framework.utils.general import copy_file, create_nc_from_cdl, config_to_dict, parse_test_list

# root_dir = os.path.dirname(os.path.abspath(__file__))
# if root_dir not in sys.path:
#     sys.path.insert(0, root_dir)

# load the functional test classes
from framework.load_functional_tests import *

# initialze CIME path
add_cime_lib_to_path()
try:
    from CIME.utils import run_cmd
except ImportError as e:
    # fail fast if environment isn't set up
    raise ImportError(
        f"CIME dependencies missing. Ensure environment is configured. Error: {e}"
    ) from e

# constants for this script
_FILE_DIR = os.path.dirname(__file__)
_DEFAULT_CONFIG_FILE = os.path.join(_FILE_DIR, "config", "functional.cfg")
_DEFAULT_CDL_PATH = os.path.abspath(
    os.path.join(
        _FILE_DIR,
        os.pardir,
        "parameter_files",
        "fates_params_default.json",
    )
)
_CMAKE_BASE_DIR = os.path.join(_FILE_DIR, os.pardir)
_TEST_SUB_DIR = "testing"


def commandline_args():
    """Parse and return command-line arguments"""

    description = """
    Driver for running FATES functional tests

    Typical usage:

    ./run_fates_functional_tests -t allometry

    """
    parser = argparse.ArgumentParser(
        description=description, formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(
        "-f",
        "--parameter-file",
        type=str,
        default=_DEFAULT_CDL_PATH,
        help="Parameter file to run the FATES tests with.\n"
        "This should be JSON formatted.\n"
        "If no file is specified the script will use the default .json file in the\n"
        "parameter_files directory.\n",
    )

    parser.add_argument(
        "--config-file",
        type=str,
        default=_DEFAULT_CONFIG_FILE,
        help=f"Configuration file where test list is defined. Default: '{_DEFAULT_CONFIG_FILE}'",
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
        "-r",
        "--run-dir",
        type=str,
        default=os.path.join(_CMAKE_BASE_DIR, "_run"),
        help="Directory where tests are run.\n"
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
        "--skip-build",
        action="store_true",
        help="Skip building and compiling the test code.\n"
        "Only do this if you already have run build.\n"
        "Script will check to make sure executables are present.\n",
    )

    parser.add_argument(
        "--skip-run-executables",
        action="store_true",
        help="Skip running test code executables.\n"
        "Only do this if you already have run the code previously.\n"
        "Script will check to make sure required output files are present.\n",
    )

    parser.add_argument(
        "--save-figs",
        action="store_true",
        help="Write out generated figures to files.\n"
        "Will be placed in run_dir/plots.\n"
        "Should probably do this on remote machines.\n",
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

    check_arg_validity(args)

    return args


def check_arg_validity(args):
    """Checks validity of input script arguments

    Args:
        args (parse_args): input arguments
    """
    # check to make sure parameter file exists and is one of the correct forms
    check_param_file(args.parameter_file)

    # make sure relevant output files exist:
    if args.skip_run_executables:
        # if you skip the run we assume you want to skip the build
        print("--skip-run specified, assuming --skip-build")
        args.skip_build = True
        check_out_files(args.run_dir, args.test_dict)

    # make sure build directory exists
    if args.skip_build:
        if args.verbose_make:
            raise argparse.ArgumentError(
                None,
                "Can't run verbose make and skip build.\n"
                "Re-run script without --skip-build",
            )
        check_build_dir(args.build_dir, args.test_dict)

    # Check that config file exists and is a file
    if not os.path.exists(args.config_file):
        raise FileNotFoundError(args.config_file)
    if not os.path.isfile(args.config_file):
        raise RuntimeError(f"config 'file' is a directory: '{args.config_file}'")


def check_param_file(param_file):
    """Checks to see if param_file exists and is of the correct form (.json)

    Args:
        param_file (str): path to parameter file

    Raises:
        argparse.ArgumentError: Parameter file is not of the correct form (.json)
        argparse.ArgumentError: Can't find parameter file
    """
    file_suffix = os.path.basename(param_file).split(".")[-1]
    if not file_suffix in ["json"]:
        raise argparse.ArgumentError(
            None, "Must supply parameter file with .json. Ending."
        )
    if not os.path.isfile(param_file):
        raise FileNotFoundError(param_file)


def check_build_dir(build_dir, test_dict):
    """Checks to see if all required build directories and executables are present

    Args:
        build_dir (str): build directory
        test_list (list, str): list of test names

    Raises:
        argparse.ArgumentError: Can't find a required build directory or executable
    """
    for attributes in test_dict.values():
        if not build_exists(build_dir, attributes["test_dir"], attributes["test_exe"]):
            raise argparse.ArgumentError(
                None,
                "Build directory or executable does not exist.\n"
                "Re-run script without --skip-build.",
            )


def check_out_files(run_dir, test_dict):
    """Checks to see that required output files are present in the run directory

    Args:
        run_dir (str): run directory
        test_dict (dict): dictionary of tests to run

    Raises:
        argparse.ArgumentError: Can't find a required output file
    """
    for test, attributes in dict(
        filter(lambda pair: pair[1]["out_file"] is not None, test_dict.items())
    ).items():
        if not os.path.isfile(
            os.path.join(os.path.abspath(run_dir), attributes["out_file"])
        ):
            raise argparse.ArgumentError(
                None,
                f"Required file for {test} test does not exist.\n"
                "Re-run script without --skip-run.",
            )


def run_functional_tests(
    clean,
    verbose_make,
    build,
    run_executables,
    build_dir,
    run_dir,
    make_j,
    param_file,
    save_figs,
    test_dict,
):
    """Builds and runs the fates functional tests

    Args:
        clean (bool): whether or not to clean the build directory
        verbose_make (bool): whether or not to run make with verbose output
        build (bool): whether or not to build the exectuables
        run_executables (bool): whether or not to run the executables
        build_dir (str): build directory
        run_dir (str): run directory
        make_j (int): number of processors for the build
        param_file (str): input FATES parameter file
        save_figs (bool): whether or not to write figures to file
        test_dict (dict): dictionary of test classes to run
    """

    # absolute path to desired build directory
    build_dir_path = os.path.abspath(build_dir)

    # absolute path to desired run directory
    run_dir_path = os.path.abspath(run_dir)

    # make run directory if it doesn't already exist
    if not os.path.isdir(run_dir_path):
        os.mkdir(run_dir_path)

    # create plot directories if we need to
    if save_figs:
        make_plotdirs(os.path.abspath(run_dir), test_dict)

    # move parameter file to correct location
    param_file = create_param_file(param_file, run_dir)
    
    # compile code
    if build:
        build_tests(
            build_dir, _CMAKE_BASE_DIR, make_j, clean=clean, verbose=verbose_make
        )

    # run executables for each test in test list
    if run_executables:
        print("Running executables")
        for _, test in test_dict.items():
            args = test.other_args
            # prepend datm file (if required) to argument list
            if isinstance(test, FunctionalTestWithDrivers) and test.datm_file:
                args.insert(0, test.datm_file)
            # prepend parameter file (if required) to argument list
            if test.use_param_file:
                args.insert(0, param_file)
            # run
            run_fortran_exectuables(
                build_dir_path, test.test_dir, test.test_exe, run_dir_path, args
            )

    # plot output for relevant tests
    for name, test in dict(
        filter(lambda pair: pair[1].plot, test_dict.items())
    ).items():
        test.plot_output(
            run_dir_path, save_figs, os.path.join(run_dir_path, "plots", name)
        )
    # show plots
    plt.show()


def make_plotdirs(run_dir, test_dict):
    """Create plotting directories if they don't already exist

    Args:
        run_dir (str): full path to run directory
        test_dict (dict): dictionary of test to run
    """
    # make main plot directory
    plot_dir = os.path.join(run_dir, "plots")
    if not os.path.isdir(plot_dir):
        os.mkdir(plot_dir)

    # make sub-plot directories
    for test in dict(filter(lambda pair: pair[1].plot, test_dict.items())):
        sub_dir = os.path.join(plot_dir, test)
        if not os.path.isdir(sub_dir):
            os.mkdir(sub_dir)


def create_param_file(param_file, run_dir):
    """Moves the input parameter file to the run directory

    Args:
        param_file (str): path to parmaeter file
        run_dir (str): full path to run directory

    Raises:
        RuntimeError: Supplied parameter file is not JSON

    Returns:
        str: full path to new parameter file name/location
    """
    if param_file is None:
        print("Using default parameter file.")
        param_file = _DEFAULT_CDL_PATH
        param_file_update = create_nc_from_cdl(param_file, run_dir)
    else:
        print(f"Using parameter file {param_file}.")
        file_suffix = os.path.basename(param_file).split(".")[-1]
        if file_suffix == "json":
            param_file_update = copy_file(param_file, run_dir)
        else:
            raise RuntimeError("Must supply parameter file with .json. Ending.")

    return param_file_update


def run_fortran_exectuables(build_dir, test_dir, test_exe, run_dir, args):
    """Run the generated Fortran executables

    Args:
        build_dir (str): full path to build directory
        run_dir (str): full path to run directory
        test_dir (str): test directory within the run directory
        test_exe (str): test executable to run
        args ([str]):   arguments for executable
    """

    # move executable to run directory
    exe_path = os.path.join(build_dir, _TEST_SUB_DIR, test_dir, test_exe)
    copy_file(exe_path, run_dir)

    # run the executable
    new_exe_path = os.path.join(run_dir, test_exe)
    run_command = [new_exe_path]
    run_command.extend(args)

    os.chdir(run_dir)
    cmd = " ".join(run_command)
    stat, out, _ = run_cmd(cmd, combine_output=True)
    if stat:
        print(out)
        raise subprocess.CalledProcessError(stat, cmd, out)
    print(out)


def get_test_subclasses(*argv):
    """
    Given a FunctionalTest* class, find all its test subclasses. Do not include child
    FunctionalTest* classes.
    """
    test_subclasses = []
    for ftest_class in argv:
        test_subclasses += [x for x in ftest_class.__subclasses__() if hasattr(x, "name")]
    return test_subclasses


def main():
    """Main script
    Reads in command-line arguments and then runs the tests.
    """

    args = commandline_args()

    full_test_dict = config_to_dict(args.config_file)
    config_dict = parse_test_list(full_test_dict, args.test_list)

    test_dict = {}

    # Get all the possible test subclasses.
    subclasses = get_test_subclasses(FunctionalTest, FunctionalTestWithDrivers)

    # Associate each test in the config file with the appropriate test subclass
    for name in config_dict.keys():
        test_class = list(filter(lambda subclass: subclass.name == name, subclasses))[
            0
        ](config_dict[name])
        test_dict[name] = test_class

    build = not args.skip_build
    run = not args.skip_run_executables

    run_functional_tests(
        args.clean,
        args.verbose_make,
        build,
        run,
        args.build_dir,
        args.run_dir,
        args.make_j,
        args.parameter_file,
        args.save_figs,
        test_dict,
    )


if __name__ == "__main__":
    main()
