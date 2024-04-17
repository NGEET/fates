#!/usr/bin/env python

import os
import sys
from build_fortran_tests import build_unit_tests

import argparse
import math
import pandas as pd
import numpy as np
import xarray as xr
import matplotlib
import matplotlib.pyplot as plt

from utils import add_cime_lib_to_path, round_up
add_cime_lib_to_path()

from CIME.utils import run_cmd_no_fail

DEFAULT_CDL_PATH = "../parameter_files/fates_params_default.cdl"
CMAKE_BASE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../")

# Constants for now
out_file = "allometry_out.nc"
test_dir = "fates_allom_test"
test_exe = "FATES_allom_exe"
name = "fates_unit_tests"

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
    
    inc = (int(maxvar) - 0)/20
    for i in range(0, 20):
        y = 0.0 + i*inc
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
    
    
def plot_total_biomass(data):
    """Plot two calculations of total biomass against each other

    Args:
        data (xarray DataSet): the allometry dataset
    """
    df = pd.DataFrame({'dbh': np.tile(data.dbh, len(data.pft)),
                       'pft': np.repeat(data.pft, len(data.dbh)),
                       'total_biomass_parts': data.total_biomass_parts.values.flatten(),
                       'total_biomass_tissues': data.total_biomass_tissues.values.flatten()})
    
    colors = get_color_pallete()
        
    plt.figure(figsize=(7, 5))
    ax = plt.subplot(111)
    ax.spines["top"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    
    maxbiomass = np.maximum(df['total_biomass_parts'].max(), df['total_biomass_tissues'].max())
    
    plt.xlim(0.0, maxbiomass)
    plt.ylim(0.0, maxbiomass)

    plt.yticks(fontsize=10)
    plt.xticks(fontsize=10)
    plt.tick_params(bottom=False, top=False, left=False, right=False)
    
    pfts = np.unique(df.pft.values)
    for rank, pft in enumerate(pfts):
        data = df[df.pft == pft]
        plt.scatter(data.total_biomass_parts.values, data.total_biomass_parts.values,
                 color=colors[rank], label=pft)
        
    plt.xlabel('Total biomass (kgC) from parts', fontsize=11)
    plt.ylabel('Total biomass (kgC) from tissues', fontsize=11)
    plt.title("Simulated total biomass for input parameter file", fontsize=11)
    plt.legend(loc='upper left', title='PFT')
  
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
    
    
def run_exectuables(build_dir, test_dir, test_exe, run_dir, args):
    """Run the generated executables

    Args:
        build_dir (str): full path to build directory
        run_dir (str): full path to run directory
        test_dir (str): test directory within the run directory
        test_exe (str): test executable to run
        args ([str]):   arguments for executable
    """
    
    # move executable to run directory
    exe_path = os.path.join(build_dir, test_dir, test_exe)
    copy_file(exe_path, run_dir)
    
    # run the executable
    new_exe_path = os.path.join(run_dir, test_exe)
    run_command = [new_exe_path]
    run_command.extend(args)
    
    os.chdir(run_dir)
    print("Running exectuables")
    out = run_cmd_no_fail(" ".join(run_command), combine_output=True)
    print(out)
    

def run_tests(clean, build, run, build_dir, run_dir, make_j, param_file):
    """Builds and runs the fates tests

    Args:
        clean (bool): whether or not to clean the build directory
        build (bool): whether or not to build the exectuables
        run (bool): whether or not to run the executables
        build_dir (str): build directory
        run_dir (str): run directory
        make_j (int): number of processors for the build
        param_file (str): input FATES parameter file

    Raises:
        RuntimeError: Parameter file is not the correct file type
    """
        
    # absolute path to desired build directory
    build_dir_path = os.path.abspath(build_dir)
    
    # absolute path to desired run directory
    run_dir_path = os.path.abspath(run_dir)
    
    if not os.path.isdir(run_dir_path):
        os.mkdir(run_dir_path)
    
    if param_file is None:
        print("Using default parameter file.")
        param_file = DEFAULT_CDL_PATH
        param_file = create_nc_file(param_file, run_dir_path)
    else:
        print(f"Using parameter file {param_file}.")
        file_suffix = os.path.basename(param_file).split(".")[-1]
        if file_suffix == 'cdl':
            param_file = create_nc_file(param_file, run_dir_path)
        elif file_suffix == "nc":
            param_file = copy_file(param_file, run_dir_path)
        else:
            raise RuntimeError("Must supply file with .cdl or .nc ending.")

    if build:
        build_unit_tests(build_dir, name, CMAKE_BASE_DIR, make_j, clean=clean)
    
    if run:
        run_exectuables(build_dir_path, test_dir, test_exe, run_dir_path, [param_file])

    # read in allometry data
    allometry_dat = xr.open_dataset(os.path.join(run_dir_path, out_file))
    
    # plot allometry data
    plot_allometry_var(allometry_dat.height, 'height', 'height', 'm')
    plot_allometry_var(allometry_dat.bagw, 'bagw', 'aboveground biomass', 'kgC')
    plot_allometry_var(allometry_dat.blmax, 'blmax', 'maximum leaf biomass', 'kgC')
    plot_allometry_var(allometry_dat.crown_area, 'crown_area', 'crown area', 'm$^2$')
    plot_allometry_var(allometry_dat.sapwood_area, 'sapwood_area', 'sapwood area', 'm$^2$')
    plot_allometry_var(allometry_dat.bsap, 'bsap', 'sapwood biomass', 'kgC')
    plot_allometry_var(allometry_dat.bbgw, 'bbgw', 'belowground biomass', 'kgC')
    plot_allometry_var(allometry_dat.fineroot_biomass, 'fineroot_biomass', 'fineroot biomass', 'kgC')
    plot_allometry_var(allometry_dat.bstore, 'bstore', 'storage biomass', 'kgC')
    plot_allometry_var(allometry_dat.bdead, 'bdead', 'deadwood biomass', 'kgC')
    plot_allometry_var(allometry_dat.total_biomass_parts, 'total_biomass_parts', 'total biomass (calculated from parts)', 'kgC')
    plot_allometry_var(allometry_dat.total_biomass_tissues, 'total_biomass_tissues', 'total biomass (calculated from tissues)', 'kgC')
    plot_total_biomass(allometry_dat)
    plt.show()


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
        default=DEFAULT_CDL_PATH,
        help="Parameter file to run the FATES tests with.\n"
        "Can be a netcdf (.nc) or cdl (.cdl) file.\n"
        "If no file is specified the script will use the default .cdl file in the\n"
        "parameter_files directory.\n",
    )
    
    parser.add_argument(
        "-b",
        "--build-dir",
        default="../_build",
        help="Directory where tests are built.\n"
        "Will be created if it does not exist.\n",
    )
    
    parser.add_argument(
        "-r",
        "--run-dir",
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
      "--skip-run",
      action="store_true",
      help="Skip running test code executables.\n"
      "Only do this if you already have run the code previously.\n"
      "Script will check to make sure required output files are present.\n",
    )
    
    args = parser.parse_args()
    
    check_arg_validity(args)
    
    return args


def check_build_exists(build_dir):
    """Checks to see if the build directory and associated executables exist.
    
        Args:
          build_dir (str): build directory
    """
   
    build_path = os.path.abspath(build_dir)
    if not os.path.isdir(build_path):
        return False
     
    exe_path = os.path.join(build_path, test_dir, test_exe)
    if not os.path.isfile(exe_path):
        return False
   
    return True


def check_out_file_exists(out_file):
    """Checks to see if the required output files exist.
    
        Args:
          out_file (str): required output file
    """
   
    full_path = os.path.abspath(out_file)
    if not os.path.isfile(full_path):
        return False
   
    return True


def check_arg_validity(args):
    """Checks validity of input script arguments

    Args:
        args (parse_args): input arguments

    Raises:
        RuntimeError: Can't find input parameter file
        RuntimeError: Can't find build directory or required executables
        RuntimeError: Can't find required output files for plotting
    """
    if args.parameter_file is not None:
        if not os.path.isfile(args.parameter_file):
            raise RuntimeError(f"Cannot find file {args.parameter_file}.")
    if args.skip_build:
        if not check_build_exists(os.path.abspath(args.build_dir)):
            raise RuntimeError("Can't find build directory or executables, run again without --skip-build")
    if args.skip_run:
        if not check_out_file_exists(os.path.join(os.path.abspath(args.run_dir), out_file)):
            raise RuntimeError(f"Can't find output file {out_file}, run again without --skip-run")

def main():
    """Main script
    """
    
    args = commandline_args()
    
    build = not args.skip_build
    run = not args.skip_run
    
    run_tests(args.clean, build, run, args.build_dir, args.run_dir, args.make_j, args.parameter_file)

if __name__ == "__main__":
    
    main()