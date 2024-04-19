"""Utility functions for file checking, math equations, etc.
"""

import math
import os
from path_utils import add_cime_lib_to_path

add_cime_lib_to_path()

from CIME.utils import run_cmd_no_fail

def round_up(num, decimals=0):
    multiplier = 10**decimals
    return math.ceil(num * multiplier)/multiplier

def truncate(num, decimals=0):
    multiplier = 10**decimals
    return int(num * multiplier)/multiplier

def create_nc_file(cdl_path, run_dir):
    """Creates a netcdf file from a cdl file

    Args:
        cdl_path (str): full path to desired cdl file
        run_dir (str): where the file should be written to
    """
    file_basename = os.path.basename(cdl_path).split(".")[-2]
    file_nc_name = f"{file_basename}.nc"

    file_gen_command = [
            "ncgen -o",
            os.path.join(run_dir, file_nc_name),
            cdl_path
    ]
    out = run_cmd_no_fail(" ".join(file_gen_command), combine_output=True)
    print(out)

    return file_nc_name

def copy_file(file_path, dir):
    """Copies a file file to a desired directory

    Args:
        file_path (str): full path to file
        dir (str): where the file should be copied to
    """
    file_basename = os.path.basename(file_path)

    file_copy_command = [
            "cp",
            os.path.abspath(file_path),
            os.path.abspath(dir)
    ]
    run_cmd_no_fail(" ".join(file_copy_command), combine_output=True)

    return file_basename

def get_color_pallete():
    """Generate a color pallete

    Returns:
        real: array of colors to use in plotting
    """

    colors = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
            (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
            (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
            (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
            (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]

    for i in range(len(colors)):
        r, g, b = colors[i]
        colors[i] = (r/255., g/255., b/255.)

    return colors