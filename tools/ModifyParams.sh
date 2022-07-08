#!/bin/sh
# Put in ~/global/scratch/users/yanlanliu/E3SM/components/elm/src/external_models/fates/tools
cp fates_parameter_default.nc fates_params_seednew_bkp.nc

./modify_fates_paramfile.py --input fates_params_seednew_bkp.nc --output fates_params_seednew_bkp.nc --O --var fates_recruit_initd --value 0.0 --allPFTs

./modify_fates_paramfile.py --input fates_params_seednew_bkp.nc --output fates_params_seednew_bkp.nc --O --var fates_recruit_initd --value 0.4 --pft 7

./modify_fates_paramfile.py --input fates_params_seednew_bkp.nc --output fates_params_seednew_bkp.nc --O --var fates_recruit_initd --value 0.4 --pft 9

./modify_fates_paramfile.py --input fates_params_seednew_bkp.nc --output fates_params_seednew_bkp.nc --O --var fates_recruit_initd --value 0.4 --pft 10


./modify_fates_paramfile.py --input fates_params_seednew_bkp.nc --output fates_params_seednew_bkp.nc --O --var fates_seed_alloc --value 0.3 --pft 7

./modify_fates_paramfile.py --input fates_params_seednew_bkp.nc --output fates_params_seednew_bkp.nc --O --var fates_seed_alloc --value 0.3 --pft 9

./modify_fates_paramfile.py --input fates_params_seednew_bkp.nc --output fates_params_seednew_bkp.nc --O --var fates_seed_alloc --value 0.3 --pft 10

./modify_fates_paramfile.py --input fates_params_seednew_bkp.nc --output fates_params_seednew_bkp.nc --O --var fates_seed_alloc_mature --value 0.7 --pft 7

./modify_fates_paramfile.py --input fates_params_seednew_bkp.nc --output fates_params_seednew_bkp.nc --O --var fates_seed_alloc_mature --value 0.7 --pft 9

./modify_fates_paramfile.py --input fates_params_seednew_bkp.nc --output fates_params_seednew_bkp.nc --O --var fates_seed_alloc_mature --value 0.7 --pft 10

./modify_fates_paramfile.py --input fates_params_seednew_bkp.nc --output fates_params_seednew_bkp.nc --O --var fates_seed_dbh_repro_threshold --value 0.5 --pft 7

./modify_fates_paramfile.py --input fates_params_seednew_bkp.nc --output fates_params_seednew_bkp.nc --O --var fates_seed_dbh_repro_threshold --value 0.5 --pft 9

./modify_fates_paramfile.py --input fates_params_seednew_bkp.nc --output fates_params_seednew_bkp.nc --O --var fates_seed_dbh_repro_threshold --value 0.5 --pft 10

./modify_fates_paramfile.py --input fates_params_seednew_bkp.nc --output fates_params_seednew_bkp.nc --O --var fates_phen_chiltemp --value -2.0

./modify_fates_paramfile.py --input fates_params_seednew_bkp.nc --output fates_params_seednew_bkp.nc --O --var fates_phen_coldtemp --value -0.0

