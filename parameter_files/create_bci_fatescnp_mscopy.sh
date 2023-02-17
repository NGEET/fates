#!/bin/sh
# =======================================================================================


SITE_BASE_DIR=/home/rgknox/Models/E3SM/cime/scripts
DIN_LOC_ROOT=/home/rgknox/Models/InputDatasets/e3sm_input_datasets
CASEROOT=/home/rgknox/LandRuns/E3SM
MACH=eddi
PROJECT=ac_ngeet
COMPILER=gnu

# SETTINGS SPECIFIC TO THIS SCRIPT (DON'T CHANGE)
CIME_MODEL=e3sm
SITE_NAME=bci_0.1x0.1_bothmet_v7i.c210414
CLM_USRDAT_DOMAIN=domain_bci_elm.c070319_c190703.nc
CLM_USRDAT_SURDAT=surfdata_allnatveg_bci_elm.c211027.nc
RES=ELM_USRDAT


# SETTINGS AND PATHS THAT ARE DEPENDENT (DONT CHANGE)
CLM_SURFDAT_DIR=${SITE_BASE_DIR}/${SITE_NAME}
CLM_DOMAIN_DIR=${SITE_BASE_DIR}/${SITE_NAME}
DIN_LOC_ROOT_FORCE=${SITE_BASE_DIR}
CLM_HASH=`git log -n 1 --pretty=%h`
FATES_HASH=`(cd ../../components/elm/src/external_models/fates;git log -n 1 --pretty=%h)`
GIT_HASH=C${CLM_HASH}-F${FATES_HASH}


# COMPSET and SPINUP/NORMAL PHASE SETTINGS
# ==============================================================================


PHASE='SPINUP'
#PHASE='TRANSIENT'

if [ $PHASE = 'SPINUP' ]
then
	 COMPSET=1850_DATM%QIA_ELM%BGC-FATES_SICE_SOCN_SROF_SGLC_SWAV
	 NAMETAG='bcnp-v72-retest-build'
fi
if [ $PHASE = 'TRANSIENT' ]
then
	 COMPSET=2000_DATM%QIA_ELM%BGC-FATES_SICE_SOCN_SROF_SGLC_SWAV
	 NAMETAG='bcnp-v72-so1-nrl1-sf1xnpp-vc55-rl3-v1.75e7-api24-p2'
fi	



CASE_NAME=${CASEROOT}/${NAMETAG}.${MACH}.${COMPILER}.${GIT_HASH}.`date +"%Y-%m-%d"`


# REMOVE EXISTING CASE DIRECTORY IF PRESENT 
# (I like using the next command, good for re-making cases
# that were originally created wrongly, but it could be dangerous
# if you forget it is there and can accidentaly remove your stuff)
rm -r ${CASE_NAME}

echo $CASE_NAME

# CREATE THE CASE
./create_newcase --case=${CASE_NAME} --res=${RES} --compset=${COMPSET} --mach=${MACH} --compiler=${COMPILER} --project=${PROJECT} 


cd ${CASE_NAME}


# SET PATHS TO SCRATCH ROOT, DOMAIN AND MET DATA (USERS WILL PROB NOT CHANGE THESE)
# =================================================================================
./xmlchange ATM_DOMAIN_FILE=${CLM_USRDAT_DOMAIN}
./xmlchange ATM_DOMAIN_PATH=${CLM_DOMAIN_DIR}
./xmlchange LND_DOMAIN_FILE=${CLM_USRDAT_DOMAIN}
./xmlchange LND_DOMAIN_PATH=${CLM_DOMAIN_DIR}
./xmlchange DATM_MODE=CLM1PT
./xmlchange ELM_USRDAT_NAME=${SITE_NAME}
./xmlchange DIN_LOC_ROOT_CLMFORC=${DIN_LOC_ROOT_FORCE}
./xmlchange CIME_OUTPUT_ROOT=${CASE_NAME}
./xmlchange NTASKS=1


# TURN ON ECA, century and nitrif/denitrif
# ==============================================================================

#./xmlchange --append ELM_BLDNML_OPTS="-nutrient cnp -nutrient_comp_pathway rd -soil_decomp century"
./xmlchange --append ELM_BLDNML_OPTS="-nutrient cnp -nutrient_comp_pathway eca -soil_decomp century"

# Set simulation start date and specify if this is spin-up or normal
# Note about spinup:
#      Make sure the spinup ALWAYS ALWAYS ALWAYS starts on year 0001 !!!
#      The normal Phase simulation must use a spinup restart file as it's "finitdat"... but,
#      it doesn't have to be from the same date as your normal phase starting date. The
#      normal phase starting date can be any time you want.

if [ $PHASE = 'SPINUP' ]
then
    ./xmlchange --append ELM_BLDNML_OPTS="-bgc_spinup on"
    ./xmlchange STOP_N=500
    ./xmlchange REST_N=1
    ./xmlchange RUN_STARTDATE='0001-01-01'
fi

if [ $PHASE = 'TRANSIENT' ]
then

    # Stop a 20tr run at 2010, and then continue
    # with SSP2-4.5 (or whatever you want)
    # to get a continuous CO2 timeseries from
    # past to future
    ./xmlchange STOP_N=510   #214
    ./xmlchange REST_N=10
    ./xmlchange RUN_STARTDATE='1501-01-01'
    ./xmlchange DATM_CO2_TSERIES=20tr
    #./xmlchange DATM_CO2_TSERIES=SSP2-4.5
    ./xmlchange CCSM_BGC=CO2A
    ./xmlchange ELM_CO2_TYPE=diagnostic
fi


# SPECIFY RUN TYPE PREFERENCES (USERS MIGHT CHANGE THESE)
# =================================================================================

./xmlchange MPILIB=mpi-serial
./xmlchange DEBUG=True
./xmlchange STOP_OPTION=nyears
./xmlchange DATM_CLMNCEP_YR_START=2003
./xmlchange DATM_CLMNCEP_YR_END=2016
./xmlchange DIN_LOC_ROOT=$DIN_LOC_ROOT



# Comment this out if you get PIO errors
# =================================================================================
./xmlchange PIO_VERSION=2


# MACHINE SPECIFIC, AND/OR USER PREFERENCE CHANGES (USERS MIGHT CHANGE THESE, JUST PREFERENCE)
# =================================================================================

./xmlchange DOUT_S_SAVE_INTERIM_RESTART_FILES=TRUE
./xmlchange DOUT_S=TRUE
./xmlchange DOUT_S_ROOT='$CASEROOT/run'
./xmlchange RUNDIR=${CASE_NAME}/run
./xmlchange EXEROOT=${CASE_NAME}/bld

# Modify the paramter file

cp /home/rgknox/Models/E3SM/components/elm/src/external_models/fates/parameter_files/fates_params_opt224_040822_v71_api25.nc fates_params_opt224_040822_local.nc

alias modparm="python /home/rgknox/Models/E3SM/components/elm/src/external_models/fates/tools/modify_fates_paramfile.py --fin fates_params_opt224_040822_local.nc --fout fates_params_opt224_040822_local.nc --O"

# ADD A FIXER TO SECOND PFT
#modparm --var fates_cnp_nfix1 --val 0.2 --pft 2

modparm --var fates_leaf_vcmax25top --val 55.7 --allpfts  # base: 30.94711
modparm --var fates_turnover_fnrt --val 3.0 --allpfts # base 1.0

# 0.304347826 * 2.5e-9 = 0.76

# def fates_eca_vmax_nh4 = 5e-09

modparm --var fates_cnp_vmax_nh4 --val 1.75e-7 --allpfts        # 2.5
modparm --var fates_cnp_vmax_no3 --val 1.75e-7 --allpfts

modparm --var fates_cnp_vmax_p --val 5e-6 --allpfts  # This is super high- only valid for a non-P, N-only run!!!


#modparm --var fates_cnp_pid_ki --pft 1 --val 5e-7   # This defaults to zero
modparm --var fates_cnp_pid_kp --allpfts --val 5e-4
modparm --var fates_cnp_pid_kd --allpfts --val 0.1
modparm --var fates_cnp_pid_ki --allpfts --val 0

modparm --var fates_cnp_nitr_store_ratio --val 1.0 --allpfts        # def 1
modparm --var fates_cnp_phos_store_ratio --val 1.0 --allpfts

modparm --var fates_cnp_store_ovrflw_frac --val 1.0 --allpfts       # def 1

modparm --var fates_allom_l2fr --val 0.4863088 --allpfts  # Starting L2FR value



# MODIFY THE CLM NAMELIST (USERS MODIFY AS NEEDED)

cat >> user_nl_elm <<EOF
fsurdat = '${SITE_BASE_DIR}/${SITE_NAME}/${CLM_USRDAT_SURDAT}'
fates_parteh_mode = 2
use_lch4 = .true.
use_fates = .true.
use_fates_nocomp = .false.
paramfile = '${SITE_BASE_DIR}/${SITE_NAME}/elm_params.cbgc.base_c201202.nc'
fates_paramfile = '${CASE_NAME}/fates_params_opt224_040822_local.nc'
hist_fincl1='FATES_STOREN_TF_CANOPY_SZPF','FATES_STOREN_TF_USTORY_SZPF',
'FATES_STOREP_TF_CANOPY_SZPF','FATES_STOREP_TF_USTORY_SZPF',
'FATES_NH4UPTAKE_SZPF','FATES_NO3UPTAKE_SZPF','FATES_PUPTAKE_SZPF',
'FATES_NEFFLUX_SZPF','FATES_VEGN_SZPF','FATES_DDBH_CANOPY_SZPF','FATES_PEFFLUX_SZPF',
'FATES_DDBH_USTORY_SZPF',
'FATES_NPP_SZPF','FATES_STOREC_TF_CANOPY_SZPF','FATES_STOREC_TF_USTORY_SZPF','FATES_FROOTCTURN_USTORY_SZ','FATES_FROOTCTURN_CANOPY_SZ'
EOF

#'FATES_L2FR_CLSZPF',
#fates_paramfile = '/home/rgknox/E3SM/components/elm/src/external_models/fates/parameter_files/fates_params_evrgn_opt224_vmaxlow_fnrt486_nfixdeva0.nc'
#fates_paramfile = '/home/rgknox/E3SM/components/elm/src/external_models/fates/parameter_files/fates_params_evrgn_opt224_vmaxhigh_fnrt486_nfixdeva.nc
#use_fates_nocomp = .false.
#fates_params_evrgn_opt224_vmaxlow_v1_011722_hydro.nc
# OLD VARS:
#hist_fincl1='STOREN_TFRAC_CANOPY_SCPF','STOREN_TFRAC_UNDERSTORY_SCPF','STOREP_TFRAC_CANOPY_SCPF','STOREP_TFRAC_UNDERSTORY_SCPF',
#'NH4UPTAKE_SCPF','NO3UPTAKE_SCPF','PUPTAKE_SCPF','NEFFLUX_SCPF','TOTVEGN_SCPF','DDBH_UNDERSTORY_SCPF','DDBH_CANOPY_SCPF',
#'LEAF2FNRT_CANOPY_SCPF','LEAF2FNRT_UNDERSTORY_SCPF','LEAF2FNRT_SCPF'

#/home/rgknox/E3SM/components/elm/src/external_models/fates/parameter_files/fates_params_evrgn_opt224_vmaxlow_l2min.5_v1_011722_hydro.nc

if [ $PHASE = 'FOCUSED' ]
then
	 echo 'hist_mfilt  = 480' >> user_nl_elm
	 echo 'hist_nhtfrq = -1' >> user_nl_elm
else
	 echo 'hist_mfilt = 12' >> user_nl_elm
	 echo 'hist_nhtfrq = 0' >> usrr_nl_elm
fi

if [ $PHASE = 'SPINUP' ]
then
    echo 'suplphos = '\'ALL\' >> user_nl_elm
    echo 'suplnitro = '\'NONE\' >> user_nl_elm
    echo 'finidat = '\'\' >> user_nl_elm
    echo 'nyears_ad_carbon_only = 30' >> user_nl_elm
else
    echo 'suplphos = '\'ALL\' >> user_nl_elm
    echo 'suplnitro = '\'NONE\' >> user_nl_elm
    echo 'finidat = '\'/home/rgknox/Models/E3SM/cime/scripts/bcnp-v72-so1-nrl1-sf1xnpp-vc55-v2.5e7-api24-p1.eddi.gnu.C98106c908-Fa3a5fde9.2022-12-27.elm.r.0501-01-01-00000.nc \' >> user_nl_elm
fi


# MODIFY THE DATM NAMELIST (DANGER ZONE - USERS BEWARE CHANGING)

cat >> user_nl_datm <<EOF
taxmode = "cycle", "cycle", "cycle","extend"
EOF

./case.setup

# HERE WE NEED TO MODIFY THE STREAM FILE (DANGER ZONE - USERS BEWARE CHANGING)
if [ $RES = ELM_USRDAT ]
then
	./preview_namelists
	cp run/datm.streams.txt.CLM1PT.ELM_USRDAT user_datm.streams.txt.CLM1PT.ELM_USRDAT
	`sed -i '/FLDS/d' user_datm.streams.txt.CLM1PT.ELM_USRDAT`
fi

./case.build
