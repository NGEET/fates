#!/usr/bin/env python3

import luh2mod as luh2

# Location of luh2 data - move this to input
# file_states = "/home/glemieux/Data/luh2/orig/states.nc"
# file_management = "/home/glemieux/Data/luh2/orig/management.nc"
# file_transitions = "/home/glemieux/Data/luh2/orig/transitions.nc"

# # Modify the files
# mod_file_states = luh2.attribupdate(file_states)
# mod_file_management = luh2.attribupdate(file_management)
# mod_file_transitions = luh2.attribupdate(file_transitions)
mod_file_states = "/home/glemieux/Data/luh2/states_modified.nc"
mod_file_management = "/home/glemieux/Data/luh2/management_modified.nc"
mod_file_transitions = "/home/glemieux/Data/luh2/transitions_modified.nc"

# Open modified files
ds_states = luh2.importdata(mod_file_states)
ds_management = luh2.importdata(mod_file_management)
ds_transitions = luh2.importdata(mod_file_transitions)

# Fix the bounds
ds_states = luh2.BoundsFixLUH2(ds_states)
ds_management= luh2.BoundsFixLUH2(ds_management)
ds_transitions = luh2.BoundsFixLUH2(ds_transitions)

# Grab surface data set to use as regrid "target"
# Import the surface dataset
file_surfdata = "/home/glemieux/Data/luh2/surfdata_4x5_16pfts_Irrig_CMIP6_simyr2000_c170824.nc"

# modify the surface data set to enable xesmf regridder to find necessary data
ds_surfdata = luh2.DimensionFixSurfData(file_surfdata)

print("done")


# def main():
# if __name__ == "__main__":
#     main()
