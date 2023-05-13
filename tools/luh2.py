#!/usr/bin/env python3

# LUH2 python script
# Usage: python luh2.py -l <luh2file> -s <luh2staticfile>

import argparse, os
from luh2mod import PrepDataSet, ImportData, SetMaskLUH2, SetMaskSurfData
from luh2mod import RegridConservative, CorrectStateSum

# Add version checking here in case environment.yml not used
def main():

    # Add argument parser - subfunction? Seperate common module?
    # input_files and range should be the only arguments
    # Allow variable input files (state and/or transitions and/or management)
    args = CommandLineArgs()


    # Prep the LUH2 datasets and regrid target
    ds_luh2 = PrepDataSet(args.luh2_file,args.begin,args.end)


    if (args.regridder_file == None):
        ds_regrid_target = PrepDataSet(args.regridder_target_file,args.begin,args.end)

        # Import the LUH2 static data to use for masking
        ds_luh2_static = ImportData(args.luh2_static_file)

        # Create new variable where the ice water fraction is inverted w
        ds_luh2_static["landfrac"] = 1 - ds_luh2_static.icwtr

        # Mask all LUH2 input data using the ice/water fraction for the LUH2 static data
        ds_luh2 = SetMaskLUH2(ds_luh2, ds_luh2_static)
        ds_luh2_static = SetMaskLUH2(ds_luh2_static, ds_luh2_static)

        # Mask the regrid target
        ds_regrid_target = SetMaskSurfData(ds_regrid_target)

        # Handle regridder file save name
        # TO DO: add check to handle if the user enters the full path
        # TO DO: check if its possible to enter nothing with the argument
        if (args.regridder_save_name == None):
            regridder_save_file = None
            print("Warning: Regridder will not be saved to file")
        else:
            output_filename = regridder_save_name

        # Regrid the luh2 data to the target grid
        # TO DO: provide a check for the save argument based on the input arguments
        regrid_luh2,regridder_luh2 = RegridConservative(ds_luh2, ds_regrid_target, regridder_save_file)

    elif (args.regridder_target_file == None):
        regridder_luh2 = ImportData(args.regridder_file)
        # TO DO: check that the time bounds match the argument bounds
        # TO DO: create bypass option to regridder function

    # # Regrid the inverted ice/water fraction data to the target grid
    #regridder_land_fraction = RegridConservative(ds_luh2_static, ds_regrid_target)
    #regrid_land_fraction = regridder_land_fraction(ds_luh2_static)
    regrid_land_fraction = regridder_luh2(ds_luh2_static)

    # Adjust the luh2 data by the land fraction
    regrid_luh2 = regrid_luh2 / regrid_land_fraction.landfrac

    # Correct the state sum (will be returned as if in not state values)
    regrid_luh2 = CorrectStateSum(regrid_luh2)

    # Add additional required variables for the host land model
    # Add 'YEAR' as a variable.
    # This is an old requirement of the HLM and should simply be a copy of the `time` dimension
    regrid_luh2["YEAR"] = regrid_luh2.time
    regrid_luh2["LONGXY"] = ds_regrid_target["LONGXY"] # TO DO: double check if this is strictly necessary
    regrid_luh2["LATIXY"] = ds_regrid_target["LATIXY"] # TO DO: double check if this is strictly necessary


    # Rename the dimensions for the output
    regrid_luh2 = regrid_luh2.rename_dims({'lat':'lsmlat','lon':'lsmlon'})

    # Merge existing regrided luh2 file with merge input target
    # TO DO: check that the grid resolution and time bounds match
    if (args.luh2_merge_file != None):
        ds_luh2_merge = ImportData(args.luh2_merge_file)
        regrid_luh2 = regrid_luhs.merge(ds_luh2_merge)

    # Write the files
    # TO DO: add check to handle if the user enters the full path
    if (args.output == None):
        output_filename = 'LUH2_timeseries.nc'
    else:
        output_filename = args.output

    output_file = os.path.join(os.getcwd(),output_filename)
    print("generating output: {}".format(output_file))
    regrid_luh2.to_netcdf(output_file)

    # Example of file naming scheme
    # finb_luh2_all_regrid.to_netcdf('LUH2_historical_1850_2015_4x5_cdk_220302.nc')

def CommandLineArgs():

    parser = argparse.ArgumentParser(description="placeholder desc")

    # Required input luh2 datafile
    # TO DO: using the checking function to report back if invalid file input
    parser.add_argument("-l","--luh2_file", required=True)

    # Provide mutually exlusive arguments for regridding input selection
    # Currently assuming that if a target is provided that a regridder file will be saved
    regrid_target = parser.add_mutually_exclusive_group(required=True)
    regrid_target.add_argument("-rf","--regridder_file") # use previously save regridder file
    regrid_target.add_argument("-rt","--regridder_target_file")  # use a dataset to regrid to

    # TO DO: static file is required if regridder file argument is not used
    # Required static luh2 data to get the ice/water fraction
    parser.add_argument("-s", "--luh2_static_file")

    # Optional argument for defining the regridder file name
    parser.add_argument("-rs", "--regridder_save_name")

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
