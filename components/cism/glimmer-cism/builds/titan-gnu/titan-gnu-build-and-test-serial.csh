#!/bin/csh



# Master build script for mac laptops. Last updated 2/28/2013 by SFP.
# This is a hacked version of Kate's original script for use on Hopper.
# For now, only supports parallel build with Trilinos using gnu and cmake. 
# Only a subset of the small, standard tests are run, on both 1 and 4 procs.
 
# (1) execute from the builds/mac-gnu subdirectory of CISM

#add logic at the top to decide which versions to build 

# PARALLEL BUILD WITH CMAKE 

# setenv TEST_DIR "/USERS/$USER/work/modeling/cism/seacism-oceans11/tests/higher-order"

# 5/7/2014 DMR -- added performance tests:

## This will automatically submit dome60-500 ijobs. gis_1km and gis_4km will not be submitted
## automatically because you will have to build and run Felix/Albany on hopper first. Once you do that,
## you can go to lines #193-194, 197-198, 201-202, and uncomment them.
setenv PERF_TEST 0

@ run_perf_tests = (($1 == run-perf-tests) || ($2 == run-perf-tests) || ($3 == run-perf-tests) || ($4 == run-perf-tests))
if ($run_perf_tests) then
  setenv PERF_TEST 1
endif

#**!move this and source it to your .bashrc (wherever your higher-order directory is located)
#setenv TEST_DIR /global/scratch2/sd/$USER/cism2/higher-order

if (! -d $TEST_DIR) mkdir -p $TEST_DIR

setenv TEST_SUITE_DEFAULT_LOC  http://oceans11.lanl.gov/cism/livv
#setenv TEST_SUITE_DEFAULT_LOC /ccs/proj/cli062/test_suite

setenv build_problem 0

set COMPILER_NAME = gnu
set PLATFORM_NAME = titan

# set PLATFORM_NAME = $1
# set COMPILER_NAME = $2

set CMAKE_SCRIPT = $PLATFORM_NAME'-'$COMPILER_NAME'-serial-cmake'
set CMAKE_CONF_OUT = 'conf_'$COMPILER_NAME'.out'
set CMAKE_BUILD_OUT = 'cmake_'$COMPILER_NAME'_build.out'
#set CISM_RUN_SCRIPT = $PLATFORM_NAME'job' 
#set CISM_RUN_SCRIPT = 'hopjob'
set CISM_RUN_SCRIPT = 'ijob' 
#set CISM_VV_SCRIPT = $PLATFORM_NAME'_VV.bash'
set CISM_VV_SCRIPT = 'rhea_VV.bash'

echo
echo 'To use this script, type: csh '$PLATFORM_NAME'-'$COMPILER_NAME'-build-and-test.csh'
echo
#echo 'For a quick test (dome only), type: csh '$PLATFORM_NAME'-'$COMPILER_NAME'-build-and-test.csh quick-test'
echo
echo "Call with no-copy to prevent copying of the reg_test and livv defaults."
echo "Call with run-perf-tests to run the performance tests."
echo "Call with skip-tests to skip testing (builds executable and copies it to TEST_DIR)."


echo
echo 'See the LIVV documentation for instructions on setting up the test directory (TEST_DIR).'
echo


#echo 'The following environment variables must be set: TEST_DIR, GLIMMER_TRILINOS_DIR'
#echo 'Examples (place in .cshrc or .bashrc):'
#echo 'csh, tcsh:  setenv GLIMMER_TRILINOS_DIR "/Users/$USER/Trilinos/gcc-build/install"'
#echo 'bash:       export GLIMMER_TRILINOS_DIR="/Users/$USER/Trilinos/gcc-build/install"'
echo
echo 'Setting TEST_DIR to the location: '
echo 'TEST_DIR =' $TEST_DIR
echo 'TEST_DIR must also be set in your .bashrc file.'

# PARALLEL BUILD WITH CMAKE

echo
echo "Configuring and building in directory: " $PWD
echo 

echo 'Configuring '$COMPILER_NAME' cmake build...'
source ./$CMAKE_SCRIPT >& $CMAKE_CONF_OUT
echo 'Making serial '$COMPILER_NAME'...'
make -j 8 >& $CMAKE_BUILD_OUT

#if ( -e example-drivers/simple_glide/src/simple_glide ) then
# echo 'Copying '$COMPILER_NAME' simple_glide_serial to test directory'
# cp -f example-drivers/simple_glide/src/simple_glide $TEST_DIR/simple_glide_serial
#else
# echo "cmake '$COMPILER_NAME' build failed, no executable"
# @ build_problem = 1
#endif

if ( -e cism_driver/cism_driver ) then
 echo 'Copying '$COMPILER_NAME' cism_driver_serial to test directory'
 cp -f cism_driver/cism_driver $TEST_DIR/cism_driver_serial
else
 echo "cmake '$COMPILER_NAME' build failed, no executable"
 @ build_problem = 1
endif


if ($build_problem == 1 ) then
  echo "No job submitted -- cmake build failed."
else  # execute tests:
 


@ no_copy_set = (($1 == no-copy) || ($2 == no-copy) || ($3 == no-copy) || ($4 == no-copy))

 # Make copy of test suite in $TEST_DIR:
if (!($no_copy_set)) then
 echo "Copying default reg_test and LIVV to $TEST_DIR"
 pushd . > /dev/null
 cd $TEST_DIR
 if ( -e reg_test_default.tgz ) rm -f reg_test_default.tgz 
 wget $TEST_SUITE_DEFAULT_LOC/reg_test_default.tgz
 tar xfz reg_test_default.tgz
 popd > /dev/null

 if ($PERF_TEST) then
    echo "Copying default perf_test to $TEST_DIR"
   pushd . > /dev/null
   cd $TEST_DIR
   if ( -e perf_test_default.tgz ) rm -f perf_test_default.tgz 
   wget $TEST_SUITE_DEFAULT_LOC/perf_test_default.tgz
   tar xfz perf_test_default.tgz
   popd > /dev/null
 endif

 cp -rf ../../tests/higher-order/livv $TEST_DIR
endif

 if (($1 == "skip-tests") || ($2 == "skip-tests") || ($3 == "skip-tests") || ($4 == "skip-tests")) then
   echo "Skipping tests."
   exit
 endif

 echo 'Submitting test jobs to compute nodes.'

 setenv run_all_tests 1
 if (($1 == "quick-test") || ($2 == "quick-test") || ($3 == "quick-test") || ($4 == "quick-test")) then
   setenv run_all_tests 0 
 endif




 #diagnostic dome test case
 cd $TEST_DIR/reg_test/dome30/diagnostic
 qsub $CISM_RUN_SCRIPT


 if ($run_all_tests == 1) then

  #evolving dome test case
  cd $TEST_DIR/reg_test/dome30/evolving
  qsub $CISM_RUN_SCRIPT

  # confined shelf to periodic BC
  cd $TEST_DIR/reg_test/confined-shelf
  qsub $CISM_RUN_SCRIPT

  # circular shelf to periodic BC
  cd $TEST_DIR/reg_test/circular-shelf
  qsub $CISM_RUN_SCRIPT

  # ISMIP test case A, 80 km 
  cd $TEST_DIR/reg_test/ismip-hom-a/80km
  qsub $CISM_RUN_SCRIPT

  # ISMIP test case A, 20 km 
  cd $TEST_DIR/reg_test/ismip-hom-a/20km
  qsub $CISM_RUN_SCRIPT

  ## ISMIP test case C, 80 km - not operational for glide
  cd $TEST_DIR/reg_test/ismip-hom-c/80km
  qsub $CISM_RUN_SCRIPT
  
  ## ISMIP test case C, 20 km - not operational for glide
  cd $TEST_DIR/reg_test/ismip-hom-c/20km
  qsub $CISM_RUN_SCRIPT
 endif

  if ($PERF_TEST == 0 ) then
    echo "No performance suite jobs were submitted."
  else
    echo 'Submitting performance jobs to compute nodes.'
    echo 'Go to rhea to complete Visualization and Verification (LIVV)'

  #dome 60 test case
    cd $TEST_DIR/perf_test/dome60
    qsub $CISM_RUN_SCRIPT

  #dome 120 test case
    cd $TEST_DIR/perf_test/dome120
    qsub $CISM_RUN_SCRIPT

  #dome 240 test case
    cd $TEST_DIR/perf_test/dome240
    qsub $CISM_RUN_SCRIPT

  #dome 500 test case
    cd $TEST_DIR/perf_test/dome500
    qsub $CISM_RUN_SCRIPT

  #dome 1000 test case - not operational currently
  #  cd $TEST_DIR/perf_test/dome1000
  #  qsub $CISM_RUN_SCRIPT
  
  #gis 4km test case
  #  cd $TEST_DIR/perf_test/gis_4km
  #  qsub $CISM_RUN_SCRIPT
  
  #gis 2km test case
  #  cd $TEST_DIR/perf_test/gis_2km
  #  qsub $CISM_RUN_SCRIPT
  
  #gis 1km test case
  #  cd $TEST_DIR/perf_test/gis_1km
  #  qsub $CISM_RUN_SCRIPT
  endif
endif


 echo
 echo "Test Suite jobs started -- using qstat to monitor."
 echo 

 set still_running = 1
 set counter = 0
 set timeout_error = 0

 set run_list = "dome_30_test dome_30_evolve conf_shelf circ_shelf ishoma_80 ishoma_20 dome_60_test dome_120_test dome_240_test dome_500_test dome_1000_test"

 while ($still_running)
  set ls_out = `qstat | grep $USER`

  set found = 0 
  foreach cur ($run_list)
   foreach elem ($ls_out)
    if ("$cur" == "$elem") then
     if (($counter % 5) == 0) echo "Still running: $cur" 
     set found = 1
    endif
    # if ($found == 1) break 
   end 
  end
  if ($found == 0) then 
   echo "All jobs completed."
   set still_running = 0
  else 
   sleep 60
  endif
  @ counter = $counter + 1
  if ($counter == 120) then
   set still_running = 0
   set timeout_error = 1
   echo "Timeout error -- jobs are taking too long. Exiting script."
  endif
  if (($counter % 5) == 0) echo "Minutes: $counter"
 end

 if ($timeout_error == 0) then
  echo "Total minutes: $counter"
  echo

  echo "Call disabled to: $CISM_VV_SCRIPT, which is located in:" 
  echo "$TEST_DIR/livv"
  echo
  echo "Perform this step on rhea after the Test Suite jobs have completed."
  # cd $TEST_DIR/livv
  # bash $CISM_VV_SCRIPT from-script $1
 endif

 echo
 # echo "If there were errors finding ncl, add the ncl installation directory to your PATH in ~/.bashrc."
 echo

endif
