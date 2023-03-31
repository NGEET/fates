#!/usr/bin/env python3

import re
import numpy as np
import xarray as xr
import xesmf as xe
from nco import Nco
from nco.custom import Atted

# Add version checking here in case environment.yml not used



# Primary function to regrid luh2 data
#
# Prepare the
def PrepDataSet(inputfile_luh2,inputfile_surface):

    # Import the data
    ds_luh2 = ImportData(inputfile_luh2)
    ds_surfdata = ImportData(inputfile_surface)

    # Correct the necessary variables for both datasets
    ds_luh2 = BoundsVariableFixLUH2(ds_luh2)
    ds_surfdata = DimensionFixSurfData(ds_surfdata)

    # Set dataset masks
    ds_luh2 = SetMask(ds_luh2)
    ds_surfdata = SetMask(ds_surfdata)

    # Define the xESMF regridder if necessary
    regridder = RegridConservative(ds_luh2,ds_surfdata)

    return(ds_luh2,ds_surfdata,regridder)

# Import luh2 or surface data sets
def ImportData(inputfile):

    # Open files
    # Check to see if a ValueError is raised which is likely due
    # to the LUH2 time units being undecodable by cftime module
    try:
       datasetout = xr.open_dataset(inputfile)
       print("Input file dataset opened: {}".format(inputfile))
       return(datasetout)
    except ValueError as err:
       print("ValueError:", err)
       errmsg = "User direction: If error is due to units being 'years since ...' " \
                "update the input data file to change to 'common_years since...'. " \
                "This can be done using the luh2.attribupdate function."
       print()
       print(errmsg)


# Modify the luh2 metadata to enable xarray to read in data
# This issue here is that the luh2 time units start prior to
# year 1672, which cftime should be able to handle, but it
# appears to need a specific unit name convention "common_years"
def AttribUpdateLUH2(inputfile,output_append="modified"):

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

# Create the necessary variable "lat_b" and "lon_b" for xESMF conservative regridding
# Each lat/lon boundary array is a 2D array corresponding to the bounds of each
# coordinate position (e.g. lat_boundary would be 90.0 and 89.75 for lat coordinate
# of 89.875).
def BoundsVariableFixLUH2(inputdataset):

    # Drop the old boundary names to avoid confusion
    outputdataset = inputdataset.drop(labels=['lat_bounds','lon_bounds'])

    # Create lat and lon bounds as a single dimension array out of the LUH2 two dimensional_bounds array.
    # Future todo: is it possible to have xESMF recognize and use the original 2D array?
    outputdataset["lat_b"] = np.insert(inputdataset.lat_bounds[:,1].data,0,inputdataset.lat_bounds[0,0].data)
    outputdataset["lon_b"] = np.insert(inputdataset.lon_bounds[:,1].data,0,inputdataset.lon_bounds[0,0].data)

    print("LUH2 dataset lat/lon boundary variables formatted and added as new variable for xESMF")

    return(outputdataset)

# The user will need to use a surface data set to regrid from, but the surface datasets
# need to have their dimensions renamed to something recognizable by xESMF
def DimensionFixSurfData(surfdataset):

    # Rename the surface dataset dimensions to something recognizable by xESMF.
    outputdataset = surfdataset.rename_dims(dims_dict={'lsmlat':'latitude','lsmlon':'longitude'})

    # Populate the new surface dataset with the actual lat/lon values
    outputdataset['longitude'] = outputdataset.LONGXY.isel(latitude=0)
    outputdataset['latitude'] = outputdataset.LATIXY.isel(longitude=0)

    print("Surface dataset dimensions renamed for xESMF")

    return(outputdataset)

def SetMask(inputdataset):

    # check what sort of inputdata is being provided; surface dataset or luh2
    # LUH2 data will need to be masked based on the variable input to mask
    dsflag,dstype = CheckDataSet(inputdataset)
    if (dsflag):
        if(dstype == "LUH2"):
            SetMaskLUH2(inputdataset,'primf') # temporary test
        elif(dstype == "Surface"):
            SetMaskSurfData(inputdataset)
        print("mask added")

# Check which dataset we're working with
def CheckDataSet(inputdataset):

    dsflag = False
    if('primf' in list(inputdataset.variables)):
        dstype = 'LUH2'
        dsflag = True
        print("LUH2")
    elif('natpft' in list(inputdataset.variables)):
        dstype = 'Surface'
        dsflag = True
        print("Surface")
    else:
        dstype = 'Unknown'
        print("Unrecognize data set")

    return(dsflag,dstype)

# LUH2 specific masking sub-function
def SetMaskLUH2(inputdataset,label_to_mask):
    # Instead of passing the label_to_mask, loop through this for all labels?
    inputdataset["mask"] = xr.where(~np.isnan(inputdataset[label_to_mask].isel(time=0)), 1, 0)
    # return(outputdataset)

# Surface dataset specific masking sub-function
def SetMaskSurfData(inputdataset):
    # Instead of passing the label_to_mask, loop through this for all labels?
    inputdataset["mask"] = inputdataset["PCT_NATVEG"]> 0
    # return(outputdataset)

def RegridConservative(dataset_from,dataset_to):
    # define the regridder transformation
    regridder = xe.Regridder(dataset_from, dataset_to, "conservative")

    return(regridder)


    # Apply regridder
    # regrid_states = regridder(finb_states)
    # regrid_management= regridder(finb_management)
    # regrid_transitions= regridder(finb_management)

    ### memory crashes on the transition data
    #fin_transitions_regrid = regridder(finb)

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
