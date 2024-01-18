import re, sys
import numpy as np
import xarray as xr
import xesmf as xe

# Import luh2 or surface data sets
def ImportLUH2TimeSeries(input_file,start=None,stop=None):

    # Open files
    # Set decode_times to false as the luh2 raw data is outside the range
    # of the standard NetCDF datetime format.
    datasetout = xr.open_dataset(input_file, cache=False, decode_times=False)
    print("Input file dataset opened: {}".format(input_file))

    # Prep the input data for use
    datasetout = PrepDataset(datasetout,start,stop)

    return(datasetout)

# Prepare the input_file to be used for regridding
def PrepDataset(input_dataset,start=None,stop=None):

    # Use the maximum span if start and stop are not present
    # This assumes that the luh2 raw data will always use a
    # 'years since' style format.

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

    # Save the timesince as a variable for future use
    input_dataset["timesince"] = time_since

    # Correct the necessary variables for both datasets
    input_dataset = _BoundsVariableFixLUH2(input_dataset)

    return(input_dataset)

# Create the necessary variable "lat_b" and "lon_b" for xESMF conservative regridding
# Each lat/lon boundary array is a 2D array corresponding to the bounds of each
# coordinate position (e.g. lat_boundary would be 90.0 and 89.75 for lat coordinate
# of 89.875).
def _BoundsVariableFixLUH2(input_dataset):

    # Create lat and lon bounds as a single dimension array out of the LUH2 two dimensional_bounds array.
    # Future todo: is it possible to have xESMF recognize and use the original 2D array?
    input_dataset["lat_b"] = np.insert(input_dataset.lat_bounds[:,1].data,0,input_dataset.lat_bounds[0,0].data)
    input_dataset["lon_b"] = np.insert(input_dataset.lon_bounds[:,1].data,0,input_dataset.lon_bounds[0,0].data)

    # Drop the old boundary names to avoid confusion
    input_dataset = input_dataset.drop(labels=['lat_bounds','lon_bounds'])

    print("LUH2 dataset lat/lon boundary variables formatted and added as new variable for xESMF")

    return(input_dataset)

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
