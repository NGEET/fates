#!/usr/bin/env python3

import re, sys
import numpy as np
import xarray as xr
import xesmf as xe

# Import luh2 or surface data sets
def ImportData(input_file,start=None,stop=None,merge_flag=False):

    # Open files
    # Set decode_times to false as the luh2 raw data is outside the range
    # of the standard NetCDF datetime format.
    datasetout = xr.open_dataset(input_file, cache=False, decode_times=False)
    print("Input file dataset opened: {}".format(input_file))

    # Prep the input data for use
    datasetout = PrepDataset(datasetout,start,stop,merge_flag)

    return(datasetout)

# Prepare the input_file to be used for regridding
def PrepDataset(input_dataset,start=None,stop=None,merge_flag=False):

    # Check that the input dataset is a valid type
    dsflag, dstype = CheckDataset(input_dataset)

    # Use the maximum span if start and stop are not present
    # This assumes that the luh2 raw data will always use a
    # 'years since' style format.
    if(not(dstype in ('static','regrid'))):

        if (dstype == 'LUH2'):
            # Get the units to determine the file time
            # It is expected that the units of time is 'years since ...'
            time_since_array = input_dataset.time.units.split()
            if (time_since_array[0] != 'years'):
                sys.exit("FileTimeUnitsError: input file units of time is not 'years since ...'")

            # Note that datetime package is not used as the date range might
            # be beyond the bounds of the packages applicable bounds
            time_since = int(time_since_array[2].split('-')[0])

            # Get the time bounds of the input file
            start_bound = input_dataset.time.values[0]
            stop_bound = input_dataset.time.values[-1]

            # If no input provided, simply get the bounds of the time
            if (isinstance(start,type(None))):
               start = start_bound + time_since

            if (isinstance(stop,type(None))):
               stop = stop_bound + time_since

            # Convert the input dates to years since 0850
            years_since_start = start - time_since
            years_since_stop = stop - time_since

            # Abort if the times provided are outside the applicable range
            if (years_since_start < start_bound or years_since_stop < start_bound or
                years_since_start > stop_bound or years_since_stop > stop_bound):
                sys.exit("StartStopBoundError: the input start or stop date is outside the applicable range of {} to {}".format(time_since+start_bound,time_since+stop_bound))

            # Truncate the data to the user defined range
            # This might need some more error handling for when
            # the start/stop is out of range
            input_dataset = input_dataset.sel(time=slice(years_since_start,years_since_stop))

        # Correct the necessary variables for both datasets
        # We don't need to Prep the incoming dataset if it's being opened to merge
        if(not merge_flag):
            input_dataset = PrepDataset_ESMF(input_dataset,dsflag,dstype)

    return(input_dataset)

# Updating datasets to work with xESMF
def PrepDataset_ESMF(input_dataset,dsflag,dstype):

    if (dsflag):
        if(dstype == "LUH2"):
            print("PrepDataset: LUH2")
            input_dataset = BoundsVariableFixLUH2(input_dataset)
        elif(dstype == "surface"):
            print("PrepDataset: SurfData")
            input_dataset = DimensionFixSurfData(input_dataset)
        print("data set updated for xESMF\n")

    return(input_dataset)

# Create the necessary variable "lat_b" and "lon_b" for xESMF conservative regridding
# Each lat/lon boundary array is a 2D array corresponding to the bounds of each
# coordinate position (e.g. lat_boundary would be 90.0 and 89.75 for lat coordinate
# of 89.875).
def BoundsVariableFixLUH2(input_dataset):

    # Create lat and lon bounds as a single dimension array out of the LUH2 two dimensional_bounds array.
    # Future todo: is it possible to have xESMF recognize and use the original 2D array?
    input_dataset["lat_b"] = np.insert(input_dataset.lat_bounds[:,1].data,0,input_dataset.lat_bounds[0,0].data)
    input_dataset["lon_b"] = np.insert(input_dataset.lon_bounds[:,1].data,0,input_dataset.lon_bounds[0,0].data)

    # Drop the old boundary names to avoid confusion
    input_dataset = input_dataset.drop(labels=['lat_bounds','lon_bounds'])

    print("LUH2 dataset lat/lon boundary variables formatted and added as new variable for xESMF")

    return(input_dataset)

# The user will need to use a surface data set to regrid from, but the surface datasets
# need to have their dimensions renamed to something recognizable by xESMF
def DimensionFixSurfData(input_dataset):

    # Rename the surface dataset dimensions to something recognizable by xESMF.
    input_dataset = input_dataset.rename_dims(dims_dict={'lsmlat':'lat','lsmlon':'lon'})

    # Populate the new surface dataset with the actual lat/lon values
    input_dataset['lon'] = input_dataset.LONGXY.isel(lat=0)
    input_dataset['lat'] = input_dataset.LATIXY.isel(lon=0)

    print("Surface dataset dimensions renamed for xESMF")

    return(input_dataset)

# LUH2 specific masking sub-function
def SetMaskLUH2(input_dataset,static_data_set):

    # Mask the luh2 data where the ice/water fraction is unity (i.e. fully ice covered gridcell)
    input_dataset["mask"] = (static_data_set.icwtr != 1)
    return(input_dataset)

# Surface dataset specific masking sub-function
def SetMaskSurfData(input_dataset):
    # Instead of passing the label_to_mask, loop through this for all labels?
    input_dataset["mask"] = input_dataset["PCT_NATVEG"] > 0
    return(input_dataset)

# Check which dataset we're working with
def CheckDataset(input_dataset):

    dsflag = False
    dsvars = list(input_dataset.variables)
    if('primf' in dsvars or
       'primf_to_secdn' in dsvars or
       any('irrig' in subname for subname in dsvars)):
        dstype = 'LUH2'
        dsflag = True
        print("LUH2")
    elif('natpft' in dsvars):
        dstype = 'surface'
        dsflag = True
        print("Surface")
    elif('icwtr' in dsvars):
        dstype = 'static'
        dsflag = True
    elif('col' in dsvars):
        dstype = 'regrid'
        dsflag = True
    else:
        dstype = 'Unknown'
        sys.exit("CheckDataSetError: Unrecognize data set")

    return(dsflag,dstype)

def RegridConservative(ds_to_regrid, ds_regrid_target, regridder_weights, regrid_reuse):

    # define the regridder transformation
    regridder = GenerateRegridder(ds_to_regrid, ds_regrid_target, regridder_weights, regrid_reuse)

    # Loop through the variables to regrid
    ds_regrid = RegridLoop(ds_to_regrid, regridder)

    return (ds_regrid, regridder)

def GenerateRegridder(ds_to_regrid, ds_regrid_target, regridder_weights_file, regrid_reuse):

    regrid_method = "conservative"
    print("\nDefining regridder, method: ", regrid_method)

    if (regrid_reuse):
        regridder = xe.Regridder(ds_to_regrid, ds_regrid_target,
                                 regrid_method, weights=regridder_weights_file)
    else:
        regridder = xe.Regridder(ds_to_regrid, ds_regrid_target, regrid_method)

        # If we are not reusing the regridder weights file, then save the regridder
        filename = regridder.to_netcdf(regridder_weights_file)
        print("regridder saved to file: ", filename)

    return(regridder)

def RegridLoop(ds_to_regrid, regridder):

    # To Do: implement this with dask
    print("\nRegridding")

    # Loop through the variables one at a time to conserve memory
    ds_varnames = list(ds_to_regrid.variables.keys())
    varlen = len(ds_to_regrid.variables)
    first_var = False
    for i in range(varlen-1):

        # Skip time variable
        if (ds_varnames[i] != "time"):

            # Only regrid variables that match the lat/lon shape.
            if (ds_to_regrid[ds_varnames[i]][0].shape == (ds_to_regrid.lat.shape[0], ds_to_regrid.lon.shape[0])):
                print("regridding variable {}/{}: {}".format(i+1, varlen, ds_varnames[i]))

                # For the first non-coordinate variable, copy and regrid the dataset as a whole.
                # This makes sure to correctly include the lat/lon in the regridding.
                if (not(first_var)):
                    ds_regrid = ds_to_regrid[ds_varnames[i]].to_dataset() # convert data array to dataset
                    ds_regrid = regridder(ds_regrid)
                    first_var = True

                # Once the first variable has been included, then we can regrid by variable
                else:
                    ds_regrid[ds_varnames[i]] = regridder(ds_to_regrid[ds_varnames[i]])
            else:
                print("skipping variable {}/{}: {}".format(i+1, varlen, ds_varnames[i]))
        else:
            print("skipping variable {}/{}: {}".format(i+1, varlen, ds_varnames[i]))

    print("\n")
    return(ds_regrid)

# Temporary: Add minor correction factor to assure states sum to one
def CorrectStateSum(input_dataset):

    # Only calculate the state sum to unity correction for the appropiate dataset
    # TO DO: Update this to use the check function
    if (not(any('irrig' in var for var in input_dataset) or
            any('_to_' in var for var in input_dataset))):

        # Drop the secma and secmb variables temporarily
        temp_dataset = input_dataset.drop({'secma','secmb'})

        # Sum the remaining state variables and normalize
        state_sum = temp_dataset.to_array().sum(dim='variable')
        state_sum = state_sum.where(state_sum != 0)
        temp_dataset = temp_dataset / state_sum

        # Update dataset with new scaled values
        input_dataset.update(temp_dataset)

        # Save the correction value
        input_dataset["stscf"] = 1.0 / state_sum

    return(input_dataset)
