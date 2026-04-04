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
import argparse
import logging
from pathlib import Path
import matplotlib.pyplot as plt
import framework.utils.env_check
from framework.loader import get_test_instances, validate_test_configs
from framework.utils.general import config_to_dict, parse_test_list, copy_file
from framework.builder import build_tests

# constants
_DEFAULT_CONFIG_FILE = Path(__file__).resolve().parents[0] / "config" / "functional.cfg"
_DEFAULT_PARAM_FILE = (
    Path(__file__).resolve().parents[1]
    / "parameter_files"
    / "fates_params_default.json"
)
_CMAKE_BASE_DIR = Path(__file__).resolve().parents[1]

# setup logging
logging.basicConfig(
    level=logging.WARNING, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


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
        default=_DEFAULT_PARAM_FILE,
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
        default=_CMAKE_BASE_DIR / "_build",
        help="Directory where tests are built.\n"
        "Will be created if it does not exist.\n",
    )

    parser.add_argument(
        "-r",
        "--run-dir",
        type=str,
        default=_CMAKE_BASE_DIR / "_run",
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
        "--skip-run",
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
        "--verbose", action="store_true", help="Run with verbose output."
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

    return parser.parse_args()


def make_plotdirs(run_dir: Path, test_instances: dict):
    """Create plotting directories

    Args:
        run_dir (Path): path to run directory
        test_instances (dict): dictionary of test instances
    """

    # create top-level plot directory
    plot_dir = run_dir / "plots"
    if not plot_dir.exists():
        plot_dir.mkdir(parents=True, exist_ok=True)

    # make sub-plot directories
    for test_name, test in test_instances.items():
        if test.plot:
            sub_dir = plot_dir / test_name
            if not sub_dir.exists():
                sub_dir.mkdir(parents=True, exist_ok=True)


def copy_param_file(param_file: Path, run_dir: Path) -> Path:
    """Moves the input parameter file to the run directory

    Args:
        param_file (Path): path to parameter file
        run_dir (Path): path to run directory

    Raises:
        ValueError: parameter file not .json file

    Returns:
        Path: path to new parameter file in run directory
    """

    logging.info("Using parameter file, %s", param_file)
    file_suffix = param_file.name.split(".")[-1]
    if file_suffix == "json":
        return copy_file(param_file, run_dir)

    raise ValueError("Must supply parameter file with .json ending.")


def prep_directories(run_dir: Path, test_instances: dict, save_figs: bool):
    """Creates directories (run_dir and plotting directories) for the tests

    Args:
        run_dir (Path): path to run directory
        test_instances (dict): dictionary of test class instances
        save_figs (bool): whether or not to save figures
    """

    # create run directory
    if not run_dir.exists():
        run_dir.mkdir(parents=True, exist_ok=True)

    # create plot directories
    if save_figs:
        make_plotdirs(run_dir, test_instances)


def main():
    """Main script
    Reads in command-line arguments and then runs tests
    """

    args = commandline_args()

    build_dir = Path(args.build_dir)
    run_dir = Path(args.run_dir)

    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)
    else:
        logging.getLogger().setLevel(logging.WARNING)

    # get config of specific tests to run
    full_test_dict = config_to_dict(args.config_file)
    config_dict = parse_test_list(full_test_dict, args.test_list)
    
    # get instances of tests to run
    test_instances = get_test_instances(config_dict)
    validate_test_configs(test_instances)

    # build tests
    build_tests(
        build_dir,
        _CMAKE_BASE_DIR,
        args.make_j,
        clean=args.clean,
        verbose=args.verbose,
    )

    # create run and plotting directories
    prep_directories(run_dir, test_instances, args.save_figs)

    # move parameter file into run directory
    param_file_moved = copy_param_file(Path(args.parameter_file), run_dir)

    # run tests
    if not args.skip_run:
        for test_name, test in test_instances.items():
            out = test.run(build_dir, run_dir, param_file_moved)
            print(out)

    # plot output
    for test_name, test in test_instances.items():
        if test.plot:
            test.plot_output(run_dir, args.save_figs, run_dir / "plots" / test_name)

    # show plots
    plt.show()


if __name__ == "__main__":
    main()
