#!/usr/bin/env python3

import re
import numpy as np
import xarray as xr
import xesmf as xe
from nco import Nco
from nco.custom import Atted

# Add version checking here in case environment.yml not used

# Import luh2 data
def importdata(inputfile):

    # Open files
    # Check to see if a ValueError is raised which is likely due
    # to the LUH2 time units being undecodable by cftime module
    try:
       datasetout = xr.open_dataset(inputfile)
    except ValueError as err:
       print("ValueError:", err)
       errmsg = "User direction: If error is due to units being 'years since ...' " \
                "update the input data file to change to 'common_years since...'. " \
                "This can be done using the luh2.attribupdate function."
       print()
       print(errmsg)

    return(datasetout)

# Modify the luh2 metadata to enable xarray to read in data
# This issue here is that the luh2 time units start prior to
# year 1672, which cftime should be able to handle, but it
# appears to need a specific unit name convention "common_years"
def attribupdate(inputfile,output_append="modified"):

    # Define the output filename
    index = inputfile.find(".nc")
    outputfile = inputfile[:index] + "_" + output_append + inputfile[index:]

    nco = Nco()

    # Get the 'time:units' string from the input using ncks
    timeunitstr = nco.ncks(input=inputfile,variable="time",options=["-m"]).decode()

    # Grab the units string and replace "years" with "common_years"
    substr = re.search('time:units.*".*"',timeunitstr)
    newstr = substr.group().replace("time:units = \"years","\"common_years")

    # Use ncatted to update the time units
    att = "units"
    var = "time"
    att_type = "c"
    opts = [" -a {0},{1},o,{2},{3}".format(att, var, att_type, newstr)]
    nco.ncatted(input=inputfile, output=outputfile, options=opts)

    print("Generated modified output file: {}".format(outputfile))

    return(outputfile)

    # The following is fixed with PR #62 for pynco but isn't in that latest update yet
    # on conda
    # nco.ncatted(input=inputfile,output=outputfile,options=[
    #     Atted(mode="overwrite",
    #           att_name="units",
    #           var_name="time",
    #           value=newstr
    #           stype="c"
    #           ),
    # ])

# Fix the boundaries of the LUH2 data
# Each lat/lon boundary array is a 2D array corresponding to the bounds of each
# coordinate position (e.g. lat_boundary would be 90.0 and 89.75 for lat coordinate
# of 89.875).
def BoundsFixLUH2():


# Make necessary (?) changes to metadata for XESMF conservative regrid
# Any LUH2 data set should work as the input dataset, but
# we should have some sort of check to make sure that the
# data sets being used are consistent in the final calling script
def MetadataUpdateLUH2(inputdataset):

    # if this is LUH2 data
    # Drop the invalid lat, lon variable labels and replace with "lat_b" and "lon_b"
    # This is only necessary for the conservative method (?) per xESMF docs

    # Is dropping these variables really necessary?  Will xESMF get confused or is this
    # for the users's sake?
    # outputdatatset = inputdataset.drop(labels=['lat_bounds','lon_bounds'])

    # Create lat and lon bounds as a single dimension array out of the LUH2 two dimensional
    # _bounds array.
    # xESMF needs these variable names for bounding the conservative regridding
    outputdatatset["lat_b"] = np.insert(inputdataset.lat_bounds[:,1].data,0,inputdataset.lat_bounds[0,0].data)
    outputdatatset["lon_b"] = np.insert(inputdataset.lon_bounds[:,1].data,0,inputdataset.lon_bounds[0,0].data)
    # outputdatatset["time"] = np.arange(len(fin["time"]), dtype=np.int16) + 850
    # outputdatatset["YEAR"] = xr.DataArray(np.arange(len(fin["time"]), dtype=np.int16) + 850, dims=("time"))

def definemask(inputdataset,outputdataset,label_to_mask):
    # make mask of LUH2 data
    outputdatatset["mask"] = xr.where(~np.isnan(inputdataset["primf_to_range"].isel(time=0)), 1, 0)
    # make mask of surface data file
    fin2b["mask"] = fin2b["PCT_NATVEG"]> 0.

# load CLM surface data file
def MetadataUpdateSurfDS(inputdataset,outputdataset):

    # Move this out to be handled by opening function
    # fin2 = xr.open_dataset('surfdata_4x5_hist_16pfts_Irrig_CMIP6_simyr2000_c190214.nc')

    # Rename the surface dataset dimensions to something recognizable by xESMF.
    # That said, since regridder just needs the values, I'm not sure this is necessary.
    # We might just be able to rename t
    fin2b = fin2.rename_dims(dims_dict={'lsmlat':'latitude','lsmlon':'longitude'})

    # Populate the new surface dataset with the actual lat/lon values
    fin2b['longitude'] = fin2b.LONGXY.isel(latitude=0)
    fin2b['latitude'] = fin2b.LATIXY.isel(longitude=0)

def RegridConservative(datasets_in,datasets_out,regridder):
    # define the regridder transformation
    regridder = xe.Regridder(dataset_in, dataset_out, "conservative")

    #regrid the various LUH2 datasets, from smallest to largest
    fin_states_regrid = regridder(finb_states)
    fin_management_regrid = regridder(finb_management)

    ### memory crashes on the transition data
    #fin_transitions_regrid = regridder(finb)

# Update the formatting to meet HLM needs
# def clmformatter():

# General functionality needed
# - collect data for specific user-defined time period
# - collect subset of the data variables (e.g. pasture, rangeland, etc)
# - write the subset to the necessary format to pass to regridding tools
#   - This may need to be specific to the particular hlm tooling
# - call regridding tooling in hlm model
#   - this will need to point at existing mapping data files, which are external to hlm
#   - this may need a new mksurf_xxxx module if it can't conform to existing modules due
#     to needing to use new fields
# future functionality:
# - importdata function for aforestation and other data to add to luh2 data
