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

This script builds and runs various FATES unit and functional tests, and plots any
relevant output from those tests.

You can supply your own parameter file (either a .cdl or a .nc file), or if you do not
specify anything, the script will use the default FATES parameter cdl file.

"""
import os
import argparse
import matplotlib.pyplot as plt
from build_fortran_tests import build_unit_tests, build_exists
from path_utils import add_cime_lib_to_path
from utils import copy_file, create_nc_file
from functional_testing.allometry.allometry_plotting import plot_allometry_dat
from functional_testing.math_utils.math_plotting import plot_quadratic_dat

add_cime_lib_to_path()

from CIME.utils import run_cmd_no_fail # pylint: disable=wrong-import-position,import-error,wrong-import-order

# Constants for this script
_DEFAULT_CDL_PATH = os.path.abspath("../parameter_files/fates_params_default.cdl")
_CMAKE_BASE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../")
_TEST_NAME = "fates_tests"
_TEST_SUB_DIR = "testing"

# Dictionary with needed constants for running the executables and reading in the
# output files - developers who add tests should add things here.

# NOTE: if the functional test you write requires a parameter file read in as a
# command-line argument, this should be the *first* (or only) argument in the
# command-line argument list
_ALL_TESTS_DICT = {
        "allometry": {
          "test_dir": "fates_allom_ftest",
          "test_exe": "FATES_allom_exe",
          "out_file": "allometry_out.nc",
          "has_unit_test": False,
          "use_param_file": True,
          "other_args": [],
          "plotting_function": plot_allometry_dat,
        },
        "quadratic": {
          "test_dir": "fates_math_ftest",
          "test_exe": "FATES_math_exe",
          "out_file": "quad_out.nc",
          "has_unit_test": False,
          "use_param_file": False,
          "other_args": [],
          "plotting_function": plot_quadratic_dat,
        },
        "fire_weather":{
            "test_dir": "fates_fire_weather_utest",
            "test_exe": None,
            "out_file": None,
            "has_unit_test": True,
            "use_param_file": False,
            "other_args": [],
            "plotting_function": None,
        }
    }

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
    out = run_cmd_no_fail(" ".join(run_command), combine_output=True)
    print(out)

def make_plotdirs(run_dir, test_dict):
    """Create plotting directories if they don't already exist

    Args:
        run_dir (str): full path to run directory
        test_dict (dict): dictionary of test to run
    """
    # make main plot directory
    plot_dir = os.path.join(run_dir, 'plots')
    if not os.path.isdir(plot_dir):
        os.mkdir(plot_dir)

    # make sub-plot directories
    for test in dict(filter(lambda pair: pair[1]['plotting_function'] is not None,
                                            test_dict.items())):
        sub_dir = os.path.join(plot_dir, test)
        if not os.path.isdir(sub_dir):
            os.mkdir(sub_dir)

def create_param_file(param_file, run_dir):
    """Creates and/or move the default or input parameter file to the run directory
    Creates a netcdf file from a cdl file if a cdl file is supplied

    Args:
        param_file (str): path to parmaeter file
        run_dir (str): full path to run directory

    Raises:
        RuntimeError: Supplied parameter file is not netcdf (.cd) or cdl (.cdl)

    Returns:
        str: full path to new parameter file name/location
    """
    if param_file is None:
        print("Using default parameter file.")
        param_file = _DEFAULT_CDL_PATH
        param_file_update = create_nc_file(param_file, run_dir)
    else:
        print(f"Using parameter file {param_file}.")
        file_suffix = os.path.basename(param_file).split(".")[-1]
        if file_suffix == 'cdl':
            param_file_update = create_nc_file(param_file, run_dir)
        elif file_suffix == "nc":
            param_file_update = copy_file(param_file, run_dir)
        else:
            raise RuntimeError("Must supply parameter file with .cdl or .nc ending.")

    return param_file_update

def run_tests(clean, verbose_make, build_tests, run_executables, build_dir, run_dir,
              make_j, param_file, save_figs, test_dict):
    """Builds and runs the fates tests

    Args:
        clean (bool): whether or not to clean the build directory
        verbose_make (bool): whether or not to run make with verbose output
        build_tests (bool): whether or not to build the exectuables
        run_executables (bool): whether or not to run the executables
        build_dir (str): build directory
        run_dir (str): run directory
        make_j (int): number of processors for the build
        param_file (str): input FATES parameter file
        save_figs (bool): whether or not to write figures to file
        test_dict (dict): dictionary of tests to run
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

    # move parameter file to correct location (creates nc file if cdl supplied)
    param_file = create_param_file(param_file, run_dir)

    # compile code
    if build_tests:
        build_unit_tests(build_dir, _TEST_NAME, _CMAKE_BASE_DIR, make_j, clean=clean,
                         verbose=verbose_make)

    # run executables for each test in test list
    if run_executables:
        print("Running executables")
        # we don't run executables for only pfunit tests
        for attributes in dict(filter(lambda pair: pair[1]['test_exe'] is not None,
                                            test_dict.items())).values():
            # prepend parameter file (if required) to argument list
            args = attributes['other_args']
            if attributes['use_param_file']:
                args.insert(0, param_file)
            # run
            run_fortran_exectuables(build_dir_path, attributes['test_dir'],
                            attributes['test_exe'], run_dir_path, args)
            
    # run unit tests 
    for test, attributes in dict(filter(lambda pair: pair[1]['has_unit_test'],
                                        test_dict.items())).items():
        print(f"Running unit tests for {test}.")
        
        test_dir = os.path.join(build_dir_path, _TEST_SUB_DIR, attributes['test_dir'])
        ctest_command = ["ctest", "--output-on-failure"]
        output = run_cmd_no_fail(" ".join(ctest_command), from_dir=test_dir, 
                                 combine_output=True)
        print(output)

    # plot output for relevant tests
    for test, attributes in dict(filter(lambda pair: pair[1]['plotting_function'] is not None,
                                            test_dict.items())).items():
        attributes['plotting_function'](run_dir_path,
                                        attributes['out_file'], save_figs,
                                        os.path.join(run_dir_path, 'plots', test))
    # show plots
    plt.show()

def out_file_exists(run_dir, out_file):
    """Checks to see if the file out_file exists in the run_dir

    Args:
        run_dir (str): full path to run directory
        out_file (str): output file name

    Returns:
        bool: yes/no file exists in correct location
    """
    return os.path.isfile(os.path.join(run_dir, out_file))

def parse_test_list(test_string):
    """Parses the input test list and checks for errors

    Args:
        test (str): user-supplied comma-separated list of test names

    Returns:
        dictionary: filtered dictionary of tests to run

    Raises:
        RuntimeError: Invalid test name supplied
    """
    valid_test_names = _ALL_TESTS_DICT.keys()

    if test_string != "all":
        test_list = test_string.split(',')
        for test in test_list:
            if test not in valid_test_names:
                raise argparse.ArgumentTypeError("Invalid test supplied, \n"
                                                 "must supply one of:\n"
                                  f"{', '.join(valid_test_names)}\n"
                                  "or do not supply a test name to run all tests.")
        test_dict = {key: _ALL_TESTS_DICT[key] for key in test_list}
    else:
        test_dict = _ALL_TESTS_DICT

    return test_dict

def commandline_args():
    """Parse and return command-line arguments"""

    description = """
    Driver for running FATES unit and functional tests

    Typical usage:

    ./run_fates_tests -f parameter_file.nc

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
        "Can be a netcdf (.nc) or cdl (.cdl) file.\n"
        "If no file is specified the script will use the default .cdl file in the\n"
        "parameter_files directory.\n",
    )

    parser.add_argument(
        "-b",
        "--build-dir",
        type=str,
        default="../_build",
        help="Directory where tests are built.\n"
        "Will be created if it does not exist.\n",
    )

    parser.add_argument(
        "-r",
        "--run-dir",
        type=str,
        default="../_run",
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
      "--verbose-make",
      action="store_true",
      help="Run make with verbose output."
    )

    parser.add_argument(
      "-t",
      "--test-list",
      action="store",
      dest="test_dict",
      type=parse_test_list,
      default="all",
      help="Test(s) to run. Comma-separated list of test names, or 'all'\n"
      "for all tests. If not supplied, will run all tests."
    )

    args = parser.parse_args()

    check_arg_validity(args)

    return args

def check_param_file(param_file):
    """Checks to see if param_file exists and is of the correct form (.nc or .cdl)

    Args:
        param_file (str): path to parameter file

    Raises:
        argparse.ArgumentError: Parameter file is not of the correct form (.nc or .cdl)
        argparse.ArgumentError: Can't find parameter file
    """
    file_suffix = os.path.basename(param_file).split(".")[-1]
    if not file_suffix in ['cdl', 'nc']:
        raise argparse.ArgumentError(None, "Must supply parameter file with .cdl or .nc ending.")
    if not os.path.isfile(param_file):
        raise argparse.ArgumentError(None, f"Cannot find file {param_file}.")

def check_build_dir(build_dir, test_dict):
    """Checks to see if all required build directories and executables are present

    Args:
        build_dir (str): build directory
        test_list (list, str): list of test names

    Raises:
        argparse.ArgumentError: Can't find a required build directory or executable
    """
    for attributes in test_dict.values():
        if not build_exists(build_dir, attributes['test_dir'], attributes['test_exe']):
            raise argparse.ArgumentError(None, "Build directory or executable does not exist.\n"
                                "Re-run script without --skip-build.")

def check_out_files(run_dir, test_dict):
    """Checks to see that required output files are present in the run directory

    Args:
        run_dir (str): run directory
        test_dict (dict): dictionary of tests to run

    Raises:
        argparse.ArgumentError: Can't find a required output file
    """
    for test, attributes in dict(filter(lambda pair: pair[1]['out_file'] is not None,
                                            test_dict.items())).items():
        if not out_file_exists(os.path.abspath(run_dir), attributes['out_file']):
            raise argparse.ArgumentError(None, f"Required file for {test} test does not exist.\n"
                                "Re-run script without --skip-run.")

def check_arg_validity(args):
    """Checks validity of input script arguments

    Args:
        args (parse_args): input arguments
    """
    # check to make sure parameter file exists and is one of the correct forms
    check_param_file(args.parameter_file)

    # make sure build directory exists
    if args.skip_build:
        if args.verbose_make:
            raise argparse.ArgumentError(None, "Can't run verbose make and skip build.\n" 
                                         "Re-run script without --skip-build")
        check_build_dir(args.build_dir, args.test_dict)

    # make sure relevant output files exist:
    if args.skip_run_executables:
        check_out_files(args.run_dir, args.test_dict)

def main():
    """Main script
      Reads in command-line arguments and then runs the tests.
    """

    args = commandline_args()

    build = not args.skip_build
    run = not args.skip_run_executables

    run_tests(args.clean, args.verbose_make, build, run, args.build_dir, args.run_dir, 
              args.make_j, args.parameter_file, args.save_figs, args.test_dict)

if __name__ == "__main__":
    main()
