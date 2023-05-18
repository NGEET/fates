#!/bin/bash
# WARNING: This script generates intermediate copies of the LUH2
# data which at its peak takes up approximately 42G of space.
#
# Note that this script must be run with the luh2 conda environment
# It requires a single argument that points to the full path location
# of the luh2 data and the dataset to regrid against

# LUH2 data names
export DATA_LOC=$1
export STATES_FILE=states.nc
export TRANSITIONS_FILE=transitions.nc
export MANAGE_FILE=management.nc
export STATIC_FILE=staticData_quarterdeg.nc
export REGRID_TARGET_FILE=surfdata_4x5_16pfts_Irrig_CMIP6_simyr2000_c170824.nc

# Save files
export REGRID_SAVE=regridder.nc
export OUTPUT=LUH2_historical_0850_2015_4x5.nc

# Combine strings
export STATES=${DATA_LOC}/${STATES_FILE}
export TRANSITIONS=${DATA_LOC}/${TRANSITIONS_FILE}
export MANAGE=${DATA_LOC}/${MANAGE_FILE}
export STATIC=${DATA_LOC}/${STATIC_FILE}
export REGRID_TARGET=${DATA_LOC}/${REGRID_TARGET_FILE}
export REGRIDDER=${DATA_LOC}/${REGRID_SAVE}

# Create copies of the luh2 data that have the time data modified
# Comment this out if the user already has the modified datasets available

# Update that filename to point to

# Regrid the luh2 data against a target surface data set and then remove the states_modified file
echo "starting storage"
du -h ${DATA_LOC}
echo "Correcting LUH2 time data for ${STATES}"
python -c "from luh2mod import AttributeUpdateLUH2; AttributeUpdateLUH2('${STATES}')"
export STATES=${DATA_LOC}/states_modified.nc
python luh2.py  -l ${STATES} -s ${STATIC} -r ${REGRID_TARGET} -w ${REGRIDDER} -o ${DATA_LOC}/states_regrid.nc
echo -e"storage status:\n"
du -h ${DATA_LOC}
rm ${STATES}

# Regrid the luh2 transitions data using the saved regridder weights file and merge into previous regrid output
echo "Correcting LUH2 time data for ${TRANSITIONS}"
python -c "from luh2mod import AttributeUpdateLUH2; AttributeUpdateLUH2('${TRANSITIONS}')"
export TRANSITIONS=${DATA_LOC}/transitions_modified.nc
python luh2.py  -l ${TRANSITIONS} -s ${STATIC} -r ${REGRID_TARGET} -w ${REGRIDDER} \
                -m ${DATA_LOC}/states_regrid.nc -o ${DATA_LOC}/states_trans_regrid.nc
echo -e"storage status:\n"
du -h ${DATA_LOC}
rm ${TRANSITIONS}
rm ${DATA_LOC}/states_regrid.nc

# Regrid the luh2 management data using the saved regridder file and merge into previous regrid output
echo "Correcting LUH2 time data for ${MANAGE}"
python -c "from luh2mod import AttributeUpdateLUH2; AttributeUpdateLUH2('${MANAGE}')"
export MANAGE=${DATA_LOC}/management_modified.nc
python luh2.py  -l ${MANAGE} -s ${STATIC} -r ${REGRID_TARGET} -w ${REGRIDDER} \
                -m ${DATA_LOC}/states_trans_regrid.nc -o ${OUTPUT}
echo -e"storage status:\n"
du -h ${DATA_LOC}
rm ${MANAGE}
rm ${DATA_LOC}/states_trans_regrid.nc
rm ${REGRIDDER}
