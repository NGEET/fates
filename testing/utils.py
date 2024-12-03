"""Utility functions for plotting, file checking, math equations, etc.
"""

import math
import os
import configparser
import argparse
import matplotlib.pyplot as plt
from path_utils import add_cime_lib_to_path

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


def get_color_palette(number: int) -> list:
    """_summary_

    Args:
        number (int): number of colors to get - must be <= 20

    Raises:
        ValueError: number must be less than hard-coded list

    Returns:
        list[tuple]: list of colors to use in plotting
    """

    # hard-coded list of colors, can add more here if necessary
    all_colors = [
        (31, 119, 180),
        (174, 199, 232),
        (255, 127, 14),
        (255, 187, 120),
        (44, 160, 44),
        (152, 223, 138),
        (214, 39, 40),
        (255, 152, 150),
        (148, 103, 189),
        (197, 176, 213),
        (140, 86, 75),
        (196, 156, 148),
        (227, 119, 194),
        (247, 182, 210),
        (127, 127, 127),
        (199, 199, 199),
        (188, 189, 34),
        (219, 219, 141),
        (23, 190, 207),
        (158, 218, 229),
    ]

    if number > len(all_colors):
        raise ValueError(f"get_color_palette: number must be <= {len(all_colors)}")

    colors = [
        (red / 255.0, green / 255.0, blue / 255.0) for red, green, blue in all_colors
    ]

    return colors[:number]


def config_to_dict(config_file: str) -> dict:
    """Convert a config file to a python dictionary

    Args:
        config_file (str): full path to config file

    Returns:
        dictionary: dictionary of config file
    """
    config = configparser.ConfigParser()
    config.read(config_file)

    dictionary = {}
    for section in config.sections():
        dictionary[section] = {}
        for option in config.options(section):
            dictionary[section][option] = config.get(section, option)

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


def blank_plot(
    x_max: float,
    x_min: float,
    y_max: float,
    y_min: float,
    draw_horizontal_lines: bool = False,
):
    """Generate a blank plot with set attributes

    Args:
        x_max (float): maximum x value
        x_min (float): minimum x value
        y_max (float): maximum y value
        y_min (float): minimum y value
        draw_horizontal_lines (bool, optional): whether or not to draw horizontal
        lines across plot. Defaults to False.
    """

    plt.figure(figsize=(7, 5))
    axis = plt.subplot(111)
    axis.spines["top"].set_visible(False)
    axis.spines["bottom"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["left"].set_visible(False)

    axis.get_xaxis().tick_bottom()
    axis.get_yaxis().tick_left()

    plt.xlim(0.0, x_max)
    plt.ylim(0.0, y_max)

    plt.yticks(fontsize=10)
    plt.xticks(fontsize=10)

    if draw_horizontal_lines:
        inc = (int(y_max) - y_min) / 20
        for i in range(0, 20):
            plt.plot(
                range(math.floor(x_min), math.ceil(x_max)),
                [0.0 + i * inc] * len(range(math.floor(x_min), math.ceil(x_max))),
                "--",
                lw=0.5,
                color="black",
                alpha=0.3,
            )

    plt.tick_params(bottom=False, top=False, left=False, right=False)
