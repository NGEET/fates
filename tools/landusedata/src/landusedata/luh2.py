import argparse, os, sys

import xarray as xr

from landusedata.luh2mod import ImportLUH2StaticFile
from landusedata.luh2mod import ImportLUH2TimeSeries

from landusedata.utils import ImportRegridTarget

from landusedata.luh2mod import SetMaskLUH2, SetMaskSurfData
from landusedata.luh2mod import RegridConservative, RegridLoop, CorrectStateSum

# Add version checking here in case environment.yml not used
def main(args):

    # Import and prep the LUH2 datasets and regrid target
    ds_luh2 = ImportLUH2TimeSeries(args.luh2_file,args.begin,args.end)

    # Import the LUH2 static data to use for masking
    ds_luh2_static = ImportLUH2StaticFile(args.luh2_static_file)

    # Create new variable where the ice water fraction is inverted w
    ds_luh2_static["landfrac"] = 1 - ds_luh2_static.icwtr

    # Mask all LUH2 input data using the ice/water fraction for the LUH2 static data
    ds_luh2 = SetMaskLUH2(ds_luh2, ds_luh2_static)
    ds_luh2_static = SetMaskLUH2(ds_luh2_static, ds_luh2_static)

    # Import and prep the regrid target surface dataset
    ds_regrid_target = ImportRegridTarget(args.regrid_target_file)

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
    # If we are merging, we might not need to do this, so check to see if its there already
    # This is a requirement of the HLM dyn_subgrid module and should be the actual year.
    # Note that the time variable from the LUH2 data is 'years since ...' so we need to
    # add the input data year
    if (not "YEAR" in list(regrid_luh2.variables)):
        regrid_luh2["YEAR"] = regrid_luh2.time + ds_luh2.timesince
        regrid_luh2["LONGXY"] = ds_regrid_target["LONGXY"] # TO DO: double check if this is strictly necessary
        regrid_luh2["LATIXY"] = ds_regrid_target["LATIXY"] # TO DO: double check if this is strictly necessary

    # Rename the dimensions for the output.  This needs to happen after the "LONGXY/LATIXY" assignment
    if (not 'lsmlat' in list(regrid_luh2.dims)):
        regrid_luh2 = regrid_luh2.rename_dims({'lat':'lsmlat','lon':'lsmlon'})

    # Reapply the coordinate attributes.  This is a workaround for an xarray bug (#8047)
    # Currently only need time
    regrid_luh2.time.attrs = ds_luh2.time.attrs
    regrid_luh2.lat.attrs = ds_luh2.lat.attrs
    regrid_luh2.lon.attrs = ds_luh2.lon.attrs

    # Merge existing regrided luh2 file with merge input target
    # TO DO: check that the grid resolution 
    # We could do this with an append during the write phase instead of the merge
    if (not(isinstance(args.luh2_merge_file,type(None)))):
        ds_luh2_merge = ImportData(args.luh2_merge_file,args.begin,args.end,merge_flag=True)
        #ds_luh2_merge = ds_luh2_merge.merge(regrid_luh2)
        regrid_luh2 = regrid_luh2.merge(ds_luh2_merge)

    # Write the files
    # TO DO: add check to handle if the user enters the full path
    output_file = os.path.join(os.getcwd(),args.output)
    print("generating output: {}".format(output_file))
    regrid_luh2.to_netcdf(output_file)

if __name__ == "__main__":
    main()
