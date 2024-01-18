import argparse, os, sys

import xarray as xr

from landusedata.luh2mod import ImportLUH2TimeSeries, CorrectStateSum
from landusedata.regrid import RegridConservative, RegridLoop
from landusedata.utils import ImportLUH2StaticFile, ImportRegridTarget
from landusedata.utils import SetMaskRegridTarget, DefineStaticMask


# Add version checking here in case environment.yml not used
def main(args):

    # Import and prep the LUH2 datasets and regrid target
    filelist = [args.luh2_states_file,
                args.luh2_transitions_file,
                args.luh2_management_file]
    ds_luh2 = []
    for filename in filelist:
        ds_luh2.append(ImportLUH2TimeSeries(filename,start=args.begin,stop=args.end))

    # Import the LUH2 static data to use for masking
    ds_luh2_static = ImportLUH2StaticFile(args.luh2_static_file)

    # Create new variable where the ice water fraction is inverted w
    ds_luh2_static["landfrac"] = 1 - ds_luh2_static.icwtr

    # Define static luh2 ice/water mask
    mask_icwtr = DefineStaticMask(ds_luh2_static)

    # Mask all LUH2 input data using the ice/water fraction for the LUH2 static data
    ds_luh2_static['mask'] = mask_icwtr
    for dataset in ds_luh2:
        dataset['mask'] = mask_icwtr

    # Import and prep the regrid target surface dataset
    ds_regrid_target = ImportRegridTarget(args.regrid_target_file)

    # Mask the regrid target
    ds_regrid_target = SetMaskRegridTarget(ds_regrid_target)

    # Create an output dataset to contain individually regridded landuse percent datasets
    ds_output = xr.Dataset()

    # Loop through the data set list, regrid and merge the sets.
    # At the beginning of the loop, there is no regrid weights file to reuse
    regrid_reuse = False
    for dataset in ds_luh2:

        # Regrid the luh2 data to the target grid
        # TO DO: provide a check for the save argument based on the input arguments
        regrid_luh2,regridder_luh2 = RegridConservative(dataset, ds_regrid_target,
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
            regrid_luh2["YEAR"] = regrid_luh2.time + dataset.timesince
            regrid_luh2["LONGXY"] = ds_regrid_target["LONGXY"] # TO DO: double check if this is strictly necessary
            regrid_luh2["LATIXY"] = ds_regrid_target["LATIXY"] # TO DO: double check if this is strictly necessary

        # Rename the dimensions for the output.  This needs to happen after the "LONGXY/LATIXY" assignment
        if (not 'lsmlat' in list(regrid_luh2.dims)):
            regrid_luh2 = regrid_luh2.rename_dims({'lat':'lsmlat','lon':'lsmlon'})

        # Reapply the coordinate attributes.  This is a workaround for an xarray bug (#8047)
        # Currently only need time
        regrid_luh2.time.attrs = dataset.time.attrs
        regrid_luh2.lat.attrs = dataset.lat.attrs
        regrid_luh2.lon.attrs = dataset.lon.attrs

        # Merge previous regrided luh2 file with merge input target
        # TO DO: check that the grid resolution matches
        # We could do this with an append during the write phase instead of the merge
        # What is the performance impact of writing to disk between dataset regrids
        # versus holding the output in memory?
        ds_output = ds_output.merge(regrid_luh2)

        # If regrid_reuse is False, change it to True for the next loop so that
        # the previously generated weights file is reused
        if (not(regrid_reuse)):
            regrid_reuse = True

    # Write the files
    # TO DO: add check to handle if the user enters the full path
    output_file = os.path.join(os.getcwd(),args.output)
    print("generating output: {}".format(output_file))
    ds_output.to_netcdf(output_file)

if __name__ == "__main__":
    main()
