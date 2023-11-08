#!/bin/bash
# WARNING: This script generates intermediate copies of the LUH2
# data which at its peak takes up approximately 42G of space.
#
# Note that this script must be run with the luh2 conda environment
# It requires a single argument that points to the full path location
# of the luh2 data and the dataset to regrid against

# LUH2 data names
DATA_LOC=$1
STATIC_LOC=$2
TARGET_LOC=$3
OUTPUT_LOC=$4
STATES_FILE=states.nc
TRANSITIONS_FILE=transitions.nc
MANAGE_FILE=management.nc
STATIC_FILE=staticData_quarterdeg.nc
REGRID_TARGET_FILE=surfdata_4x5_16pfts_Irrig_CMIP6_simyr2000_c170824.nc

START=1850
END=2015

# Save files
REGRID_SAVE=regridder.nc
OUTPUT_FILE=LUH2_historical_0850_2015_4x5.nc

# Combine strings
STATES=${DATA_LOC}/${STATES_FILE}
TRANSITIONS=${DATA_LOC}/${TRANSITIONS_FILE}
MANAGE=${DATA_LOC}/${MANAGE_FILE}
STATIC=${STATIC_LOC}/${STATIC_FILE}
REGRID_TARGET=${TARGET_LOC}/${REGRID_TARGET_FILE}
REGRIDDER=${OUTPUT_LOC}/${REGRID_SAVE}

# Comment this out if the user already has the modified datasets available

# Regrid the luh2 data against a target surface data set and then remove the states_modified file
echo "starting storage"
du -h ${OUTPUT_LOC}
python luh2.py  -b ${START} -e ${END} -l ${STATES} -s ${STATIC} -r ${REGRID_TARGET} -w ${REGRIDDER} -o ${OUTPUT_LOC}/states_regrid.nc
echo -e"storage status:\n"
du -h ${OUTPUT_LOC}

# Regrid the luh2 transitions data using the saved regridder weights file and merge into previous regrid output
python luh2.py  -b ${START} -e ${END} -l ${TRANSITIONS} -s ${STATIC} -r ${REGRID_TARGET} -w ${REGRIDDER} \
                -m ${OUTPUT_LOC}/states_regrid.nc -o ${OUTPUT_LOC}/states_trans_regrid.nc
echo -e"storage status:\n"
du -h ${OUTPUT_LOC}
rm ${DATA_LOC}/states_regrid.nc

# Regrid the luh2 management data using the saved regridder file and merge into previous regrid output
python luh2.py  -b ${START} -e ${END} -l ${MANAGE} -s ${STATIC} -r ${REGRID_TARGET} -w ${REGRIDDER} \
                -m ${OUTPUT_LOC}/states_trans_regrid.nc -o ${OUTPUT_LOC}/${OUTPUT_FILE}
echo -e"storage status:\n"
du -h ${OUTPUT_LOC}
rm ${OUTPUT_LOC}/states_trans_regrid.nc
rm ${REGRIDDER}
