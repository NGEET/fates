#!/usr/bin/env python3

import re
import numpy as np
import xarray as xr
import xesmf as xe
from nco import Nco
from nco.custom import Atted

# Add version checking here in case environment.yml not used

# Import luh2 data
def importdata(inputarg):

    # Open files
    # Check to see if a ValueError is raised which is likely due
    # to the LUH2 time units being undecodable by cftime module
    try:
       datasetout = xr.open_dataset(inputarg)
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
def attribupdate(inputarg,outputarg):

    nco = Nco()

    # Get the 'time:units' string from the input using ncks
    timeunitstr = nco.ncks(input=inputarg,variable="time",options=["-m"]).decode()

    # Grab the units string and replace "years" with "common_years"
    substr = re.search('time:units.*".*"',timeunitstr)
    newstr = substr.group().replace("time:units = \"years","\"common_years")

    # Use ncatted to update the time units
    att = "units"
    var = "time"
    att_type = "c"
    opts = [" -a {0},{1},o,{2},{3}".format(att, var, att_type, newstr)]
    nco.ncatted(input=inputarg, output=outputarg, options=opts)

    # The following is fixed with PR #62 for pynco but isn't in that latest update yet
    # on conda
    # nco.ncatted(input=inputarg,output=outputarg,options=[
    #     Atted(mode="overwrite",
    #           att_name="units",
    #           var_name="time",
    #           value=newstr
    #           stype="c"
    #           ),
    # ])

# Make necessary (?) changes to metadata for XESMF conservative regrid
# Any LUH2 data set should work as the input dataset, but
# we should have some sort of check to make sure that the
# data sets being used are consistent in the final calling script
def metadatachange(inputdataset,outputdataset):

    # Drop the invalid lat, lon label names
    outputdatatset = inputdataset.drop(labels=['lat_bounds','lon_bounds'])
    outputdatatset["lat_b"] = np.insert(inputdataset.lat_bounds[:,1].data,0,inputdataset.lat_bounds[0,0].data)
    outputdatatset["lon_b"] = np.insert(inputdataset.lon_bounds[:,1].data,0,inputdataset.lon_bounds[0,0].data)
    # outputdatatset["time"] = np.arange(len(fin["time"]), dtype=np.int16) + 850
    # outputdatatset["YEAR"] = xr.DataArray(np.arange(len(fin["time"]), dtype=np.int16) + 850, dims=("time"))

    # make mask of LUH2 data
    outputdatatset["mask"] = xr.where(~np.isnan(inputdataset["primf_to_range"].isel(time=0)), 1, 0)

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
