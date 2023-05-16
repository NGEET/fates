#!/bin/bash

# Note that this script must be run with the luh2 conda environment

# Regrid the luh2 data against a target surface data set
python luh2.py  -l ~/Data/luh2/states.nc -s ~/Data/luh2/staticData_quarterdeg.nc \
                -rt ~/Data/luh2/surfdata_4x5_16pfts_Irrig_CMIP6_simyr2000_c170824.nc \
                -rs regridder.nc \
                -o LUH2_historical_0850_2015_4x5.nc

# Regrid the luh2 transitions data using the saved regridder file and merge into previous regrid output
python luh2.py  -l ~/Data/luh2/transitions.nc \
                -rf regridder.nc \
                -m LUH2_historical_0850_2015_4x5.nc \
                -o LUH2_historical_0850_2015_4x5.nc

# Regrid the luh2 management data using the saved regridder file and merge into previous regrid output
python luh2.py  -l ~/Data/luh2/management.nc \
                -rf regridder.nc \
                -m LUH2_historical_0850_2015_4x5.nc \
                -o LUH2_historical_0850_2015_4x5.nc
