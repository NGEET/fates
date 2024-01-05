#!/bin/bash

data_loc=/home/glemieux/Data/luh2
static=${data_loc}/staticData_quarterdeg.nc
forest=${data_loc}/lu-pft/CLM5_current_luhforest_deg025.nc
pasture=${data_loc}/lu-pft/CLM5_current_luhpasture_deg025.nc
other=${data_loc}/lu-pft/CLM5_current_luhother_deg025.nc
surface=${data_loc}/lu-pft/CLM5_current_surf_deg025.nc
regrid=${data_loc}/surfdata_4x5_16pfts_Irrig_CMIP6_simyr2000_c170824.nc

python landusepft.py -c ${static} -f ${forest} -p ${pasture} -o ${other} -s ${surface} -r ${regrid} -O 'fates_landuse_pft_map_4x5.nc'
