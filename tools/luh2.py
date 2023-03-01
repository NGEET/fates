#!/usr/bin/env python3
#
import xarray as xr
import pynco

def importdata(inputarg):

    # Open files
    # Check to see if a ValueError is raised which is likely due
    # to the LUH2 time units being undecodable by cftime module
    try:
       ds = xr.open_dataset(inputarg)
    except ValueError as err:
       print("ValueError:", err)
       errmsg = "User direction: If error is due to units being 'years since ...' " \
                "update the input data file to change to 'common_years since...'. " \
                "This can be done using the NCO tool `ncatted`."
       print()
       print(errmsg)
       # exit(2)

# Close files
# management_ds.close()
