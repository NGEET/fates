#!/usr/bin/env python3

import argparse, os
from luh2mod import PrepDataSet, ImportData, SetMaskLUH2, SetMaskSurfData
from luh2mod import RegridConservative

# Add version checking here in case environment.yml not used
def main():

    # Add argument parser - subfunction? Seperate common module?
    # input_files and range should be the only arguments
    # Allow variable input files (state and/or transitions and/or management)
    args = CommandLineArgs()

    # Prep the LUH2 datasets and regrid target
    # To Do: need logic to handle taking in a saved regrid file
    ds_luh2 = PrepDataSet(args.luh2_file,args.begin,args.end)
    ds_regrid_target= PrepDataSet(args.regridder_target_file,args.begin,args.end)

    # Import the LUH2 static data to use for masking
    ds_luh2_static = ImportData(args.luh2_static_file)

    # Create new variable where the ice water fraction is inverted
    ds_luh2_static["landfrac"] = 1 - ds_luh2_static.icwtr

    # Mask all LUH2 input data using the ice/water fraction for the LUH2 static data
    ds_luh2 = SetMaskLUH2(ds_luh2, ds_luh2_static)
    ds_luh2_static = SetMaskLUH2(ds_luh2_static, ds_luh2_static)

    # Mask the regrid target
    ds_regrid_target = SetMaskSurfData(ds_regrid_target)

    # Regrid the luh2 data to the target grid
    regridder_luh2 = RegridConservative(ds_luh2, ds_regrid_target, save=True)
    regrid_luh2 = regridder_luh2(ds_luh2)

    # # Regrid the inverted ice/water fraction data to the target grid
    # # regridder_land_fraction = RegridConservative(ds_luh2_static, ds_regrid_target)
    # # regrid_land_fraction = regridder_land_fraction(ds_luh2_static)
    # regrid_land_fraction = regridder_luh2(ds_luh2_static)

    # # Adjust the luh2 data by the land fraction
    # # To Do: check if we need to do this for transition and management data as well
    # regrid_luh2 = regrid_luh2 / regrid_land_fraction.landfrac

    # # Rename the dimensions for the output
    # regrid_luh2 = regrid_luh2.rename_dims(dims_dict={'latitude':'lsmlat','longitude':'lsmlon'})
    # regrid_luh2["LONGXY"] = ds_regrid_target["LONGXY"]
    # regrid_luh2["LATIXY"] = ds_regrid_target["LATIXY"]

    # # Add 'YEAR' as a variable.  This is an old requirement of the HLM and should simply be a copy of the `time` dimension
    # regrid_luh2["YEAR"] = regrid_luh2.time

    # # Write the files
    # # TO DO: add check to handle if the user enters the full path
    # output_filename = args.output
    # if (args.output == None):
    #     output_filename = 'LUH2_timeseries.nc'

    # output_file = os.path.join(os.getcwd(),output_filename)
    # regrid_luh2.to_netcdf(output_file)

    # Example of file naming scheme
    # finb_luh2_all_regrid.to_netcdf('LUH2_historical_1850_2015_4x5_cdk_220302.nc')

def CommandLineArgs():

    parser = argparse.ArgumentParser(description="placeholder desc")

    # Required input luh2 datafile
    # TO DO: using the checking function to report back if invalid file input
    parser.add_argument("-l","--luh2_file", required=True)

    # Required static luh2 data to get the ice/water fraction
    parser.add_argument("-s", "--luh2_static_file", required=True)

    # Provide mutually exlusive arguments for regridding input selection
    # Currently assuming that if a target is provided that a regridder file will be saved
    regrid_target = parser.add_mutually_exclusive_group(required=True)
    regrid_target.add_argument("-rf","--regridder_file") # use previously save regridder file
    regrid_target.add_argument("-rt","--regridder_target_file")  # use a dataset to regrid to

    # Optional input to subset the time range of the data
    parser.add_argument("-b","--begin")
    parser.add_argument("-e","--end")

    # Optional output argument
    parser.add_argument("-o","--output")

    args = parser.parse_args()

    return(args)


if __name__ == "__main__":
    main()
