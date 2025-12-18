"""Utility functions for file checking, math equations, etc.
Do not include any third-party modules here.
"""

import math
import os
import re
import configparser
import argparse
from framework.utils.path import add_cime_lib_to_path

add_cime_lib_to_path()

from CIME.utils import (
    run_cmd_no_fail,
)  # pylint: disable=wrong-import-position,import-error,wrong-import-order


def round_up(num: float, decimals: int = 0) -> float:
    """Rounds a number up

    Args:
        num (float): number to round
        decimals (int, optional): number of decimals to round to. Defaults to 0.

    Returns:
        float: input number rounded up
    """
    multiplier = 10**decimals
    return math.ceil(num * multiplier) / multiplier


def truncate(num: float, decimals: int = 0) -> float:
    """Rounds a number down

    Args:
        num (float): number to round
        decimals (int, optional): Decimals to round down to. Defaults to 0.

    Returns:
        float: number rounded down
    """
    multiplier = 10**decimals
    return int(num * multiplier) / multiplier


def create_nc_from_cdl(cdl_path: str, run_dir: str) -> str:
    """Creates a netcdf file from a cdl file and return path to new file.

    Args:
        cdl_path (str): full path to desired cdl file
        run_dir (str): where the file should be written to
    """
    file_basename = os.path.basename(cdl_path).split(".")[-2]
    file_nc_name = f"{file_basename}.nc"

    file_gen_command = ["ncgen -o", os.path.join(run_dir, file_nc_name), cdl_path]
    out = run_cmd_no_fail(" ".join(file_gen_command), combine_output=True)
    print(out)

    return file_nc_name


def copy_file(file_path: str, directory) -> str:
    """Copies a file file to a desired directory and returns path to file.

    Args:
        file_path (str): full path to file
        dir (str): where the file should be copied to
    """
    file_basename = os.path.basename(file_path)

    file_copy_command = ["cp", os.path.abspath(file_path), os.path.abspath(directory)]
    out = run_cmd_no_fail(" ".join(file_copy_command), combine_output=True)
    print(out)

    return file_basename


def get_abspath_from_config_file(relative_path, config_file):
    """
    Gets the absolute path of a file relative to the config file where it was defined.

    Args:
      relative_path: The path to the target file, relative to the base file.
      config_file: The path to the config file.

    Returns:
      The absolute path of the target file.
    """

    # Do nothing if it's already a absolute path
    if os.path.isabs(relative_path):
        return relative_path

    base_dir = os.path.dirname(os.path.abspath(config_file))
    absolute_path = os.path.abspath(os.path.join(base_dir, relative_path))
    return absolute_path


def config_to_dict(config_file: str) -> dict:
    """Convert a config file to a python dictionary

    Args:
        config_file (str): full path to config file

    Returns:
        dictionary: dictionary of config file
    """

    # Define list of config file options that we expect to be paths
    options_that_are_paths = ["datm_file"]

    config = configparser.ConfigParser()
    config.read(config_file)

    dictionary = {}
    for section in config.sections():
        dictionary[section] = {}
        for option in config.options(section):
            value = config.get(section, option)

            # If the option is one that we expect to be a path, ensure it's an absolute path.
            if option in options_that_are_paths:
                value = get_abspath_from_config_file(value, config_file)

            # Save value to dictionary
            dictionary[section][option] = value

    return dictionary


def parse_test_list(full_test_dict, test_string):
    """Parses the input test list and checks for errors

    Args:
        test (str): user-supplied comma-separated list of test names

    Returns:
        dictionary: filtered dictionary of tests to run

    Raises:
        RuntimeError: Invalid test name supplied
    """
    valid_test_names = full_test_dict.keys()

    if test_string != "all":
        test_list = test_string.split(",")
        for test in test_list:
            if test not in valid_test_names:
                raise argparse.ArgumentTypeError(
                    "Invalid test supplied, \n"
                    "must supply one of:\n"
                    f"{', '.join(valid_test_names)}\n"
                    "or do not supply a test name to run all tests."
                )
        test_dict = {key: full_test_dict[key] for key in test_list}
    else:
        test_dict = full_test_dict

    return test_dict


def str_to_bool(val: str) -> bool:
    """Convert a string representation of truth to True or False.

    Args:
        val (str): input string

    Raises:
        ValueError: can't figure out what the string should be converted to

    Returns:
        bool: True or False
    """
    if val.lower() in ("y", "yes", "t", "true", "on", "1"):
        return True
    if val.lower() in ("n", "no", "f", "false", "off", "0"):
        return False
    raise ValueError(f"invalid truth value {val}")

def snake_to_camel(snake_str: str) -> str:
    """Convert a snake_case string to CamelCase.

    Args:
        snake_str (str): input snake case

    Returns:
        str: output CamelCase
    """
    return "".join(word.capitalize() for word in snake_str.split("_"))

def camel_to_snake(camel_str: str) -> str:
    """Convert a CamelCase string to snake_case.

    Args:
        camel_str (str): input camel case

    Returns:
        str: output snake_case
    """
    return re.sub(r'(?<!^)(?=[A-Z])', '_', camel_str).lower()

def str_to_list(val: str) -> list:
    """converts string representation of list to actual list

    Args:
        val (str): string representation of list

    Returns:
        list: actual list
    """
    if val in ("", "[]"):
        # empty list
        return []
    res = val.strip("][").split(",")
    return [n.strip() for n in res]
