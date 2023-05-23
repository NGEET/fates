#!/usr/bin/env python3

# LUH2 python script
# Usage: python luh2.py -l <luh2file> -s <luh2staticfile>

import argparse, os
from cftime import date2num
from luh2mod import ImportData, SetMaskLUH2, SetMaskSurfData
from luh2mod import RegridConservative, RegridLoop, CorrectStateSum

# Add version checking here in case environment.yml not used
def main():

    # Add argument parser - subfunction? Seperate common module?
    # input_files and range should be the only arguments
    # Allow variable input files (state and/or transitions and/or management)
    args = CommandLineArgs()

    # Import and prep the LUH2 datasets and regrid target
    ds_luh2 = ImportData(args.luh2_file,args.begin,args.end)
    ds_regrid_target = ImportData(args.regridder_target_file,args.begin,args.end)

    # Import the LUH2 static data to use for masking
    ds_luh2_static = ImportData(args.luh2_static_file)

    # Create new variable where the ice water fraction is inverted w
    ds_luh2_static["landfrac"] = 1 - ds_luh2_static.icwtr

    # Mask all LUH2 input data using the ice/water fraction for the LUH2 static data
    ds_luh2 = SetMaskLUH2(ds_luh2, ds_luh2_static)
    ds_luh2_static = SetMaskLUH2(ds_luh2_static, ds_luh2_static)

    # Mask the regrid target
    ds_regrid_target = SetMaskSurfData(ds_regrid_target)

    # Determine if we are saving a new regridder or using an old one
    # TO DO: add check to handle if the user enters the full path
    # TO DO: check if its possible to enter nothing with the argument
    regrid_reuse = False
    # If we are merging files together, we assume that the weights file
    # being supplied exists on file
    if (not isinstance(args.luh2_merge_file,type(None))):
        regrid_reuse = True

    # Regrid the luh2 data to the target grid
    # TO DO: provide a check for the save argument based on the input arguments
    regrid_luh2,regridder_luh2 = RegridConservative(ds_luh2, ds_regrid_target,
                                                    args.regridder_weights, regrid_reuse)

    # Regrid the inverted ice/water fraction data to the target grid
    regrid_land_fraction = regridder_luh2(ds_luh2_static)

    # Adjust the luh2 data by the land fraction
    # TO DO: determine if this is necessary for the transitions and management data
    regrid_luh2 = regrid_luh2 / regrid_land_fraction.landfrac

    # Correct the state sum (checks if argument passed is state file in the function)
    regrid_luh2 = CorrectStateSum(regrid_luh2)

    # Add additional required variables for the host land model
    # Add 'YEAR' as a variable.
    # This is an old requirement of the HLM and should simply be a copy of the `time` dimension
    # If we are merging, we might not need to do this, so check to see if its there already
    if (not "YEAR" in list(regrid_luh2.variables)):
        regrid_luh2["YEAR"] = date2num(regrid_luh2.time,'common_years since 0-00-00 00:00:00')
        regrid_luh2["LONGXY"] = ds_regrid_target["LONGXY"] # TO DO: double check if this is strictly necessary
        regrid_luh2["LATIXY"] = ds_regrid_target["LATIXY"] # TO DO: double check if this is strictly necessary

    # Rename the dimensions for the output.  This needs to happen after the "LONGXY/LATIXY" assignment
    if (not 'lsmlat' in list(regrid_luh2.dims)):
        regrid_luh2 = regrid_luh2.rename_dims({'lat':'lsmlat','lon':'lsmlon'})

    # Merge existing regrided luh2 file with merge input target
    # TO DO: check that the grid resolution 
    # We could do this with an append during the write phase instead of the merge
    if (not(isinstance(args.luh2_merge_file,type(None)))):
        ds_luh2_merge = ImportData(args.luh2_merge_file,args.begin,args.end,merge_flag=True)
        #ds_luh2_merge = ds_luh2_merge.merge(regrid_luh2)
        regrid_luh2 = regrid_luh2.merge(ds_luh2_merge)

    # Write the files
    # TO DO: add check to handle if the user enters the full path
    # TO DO: add check if user sets the output the same as the input
    if (args.output == None):
        output_filename = 'LUH2_timeseries.nc'
    else:
        output_filename = args.output


    output_file = os.path.join(os.getcwd(),output_filename)
    print("generating output: {}".format(output_file))
    regrid_luh2.to_netcdf(output_file)

def CommandLineArgs():

    parser = argparse.ArgumentParser(description="placeholder desc")

    # Required input luh2 datafile
    # TO DO: using the checking function to report back if invalid file input
    parser.add_argument("-l","--luh2_file", required=True)

    # Required static luh2 data to get the ice/water fraction for masking
    parser.add_argument("-s", "--luh2_static_file", required=True)

    # File to use as regridder target (e.g. a surface dataset)
    parser.add_argument("-r","--regridder_target_file", required=True)

    # Filename to use or save for the regridder weights
    parser.add_argument("-w", "--regridder_weights", required=True)

    # Optional input to subset the time range of the data
    parser.add_argument("-b","--begin")
    parser.add_argument("-e","--end")

    # Optional output argument
    parser.add_argument("-o","--output")

    # Optional merge argument to enable merging of other files
    parser.add_argument("-m", "--luh2_merge_file")

    args = parser.parse_args()

    return(args)


if __name__ == "__main__":
    main()
