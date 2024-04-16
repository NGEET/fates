#!/usr/bin/env python

import os
import sys
from build_fortran_tests import build_unit_tests

import math
import pandas as pd
import numpy as np
import xarray as xr
import matplotlib
import matplotlib.pyplot as plt

from utils import add_cime_lib_to_path, round_up
add_cime_lib_to_path()

from CIME.utils import run_cmd_no_fail

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


def plot_allometry_var(data, var, varname, units):
    """Plot an allometry variable

    Args:
        data (xarray DataArray): the data array of the variable to plot
        var (str): variable name (for data structure)
        varname (str): variable name for plot labels
        units (str): variable units for plot labels
    """
    df = pd.DataFrame({'dbh': np.tile(data.dbh, len(data.pft)),
                       'pft': np.repeat(data.pft, len(data.dbh)),
                       var: data.values.flatten()})
    
    maxdbh = df['dbh'].max()
    maxvar = round_up(df[var].max())
    
    colors = get_color_pallete()
    
    plt.figure(figsize=(7, 5))
    ax = plt.subplot(111)
    ax.spines["top"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

    plt.xlim(0.0, maxdbh)
    plt.ylim(0.0, maxvar)

    plt.yticks(fontsize=10)
    plt.xticks(fontsize=10)
    
    for y in range(0, int(maxvar), 5):
        plt.plot(range(math.floor(0), math.ceil(maxdbh)),
                  [y] * len(range(math.floor(0), math.ceil(maxdbh))),
                  "--", lw=0.5, color="black", alpha=0.3)
        
    plt.tick_params(bottom=False, top=False, left=False, right=False)
    
    pfts = np.unique(df.pft.values)
    for rank, pft in enumerate(pfts):
        data = df[df.pft == pft]
        plt.plot(data.dbh.values, data[var].values, lw=2, color=colors[rank],
                 label=pft)
        
    plt.xlabel('DBH (cm)', fontsize=11)
    plt.ylabel(f'{varname} ({units})', fontsize=11)
    plt.title(f"Simulated {varname} for input parameter file", fontsize=11)
    plt.legend(loc='upper left', title='PFT')
    
    plt.show()
    
def create_nc_file(cdl_path):
    file_basename = os.path.basename(cdl_path).split(".")[-2]
    file_nc_name = f"{file_basename}.nc"
    
    file_gen_command = [
            "ncgen -o",
            os.path.join(file_nc_name),
            cdl_path
    ]
    run_cmd_no_fail(" ".join(file_gen_command), combine_output=True)

def main(clean, build, run, build_dir, make_j, param_file):
    
    # Constants for now
    out_file = "allometry_out.nc"
    test_dir = "fates_allom_test"
    test_exe = "FATES_allom_exe"
    name = "fates_unit_tests"
    default_cdl_path = "../parameter_files/fates_params_default.cdl"
    
    # absolute path to desired build directory
    build_dir_path = os.path.abspath(build_dir)
    
    if param_file is None:
      print("Using default parameter file.")
      param_file = default_cdl_path
      create_nc_file(param_file)
    else:
      print("Using parameter file {param_file}.")

    if build:
        build_unit_tests(build_dir, name, os.path.abspath("../"), make_j, clean=clean)
        
    if run:
        exe_path = os.path.join(build_dir_path, test_dir, test_exe)
        run_command = [exe_path, os.path.abspath(param_file)]
        out = run_cmd_no_fail(" ".join(run_command), combine_output=True)
        print(out)
        
    # read in allometry data
    allometry_dat = xr.open_dataset(out_file)
    
    # plot allometry data
    plot_allometry_var(allometry_dat.height, 'height', 'height', 'm')
    
    
if __name__ == "__main__":
    
    ## Arguments
    clean = False
    build = False
    run = True
    build_dir = "../_build"
    make_j = 8
    param_file = None
    
    main(clean, build, run, build_dir, make_j, param_file)