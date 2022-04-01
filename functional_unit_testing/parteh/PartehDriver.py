# =======================================================================================
#
## @package PARTEH (Plant Allocatoin and Reactive Transport Exensible Hypotheses
#
# For usage: $python PartehDriver.py --help
#
# This script is designed to run PARTEH offline (ie not coupled with an ecosystem model).
# It will interpret user input, and provide synthetic initial conditions and boundary
# conditions to the plant.
#
# Step 1) Read in User arguments
#      1a) Define simulation conditions (initial conditions,timing,parameters,etc)
#      1b) Define state variables
#      1c) Define fluxes terms (and their forms)
#      1d) Define source-sink (boundary conditions) terms
# Step 2) Cycle through flux terms, perform allocations and determine construction of Dx/Dt
# Step 3) Initialize Simulation
# Step 4) Time-step simulation
#      4a) calculate derivative
#      4b) Integrate (either internally or via numerical integration package)
#
# =======================================================================================

import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt
from datetime import datetime
import argparse
import platform
import numpy as np
import os
import sys
import getopt
import code  # For development: code.interact(local=dict(globals(), **locals()))
import time
import imp
import ctypes
from ctypes import *
from operator import add

PartehInterpretParameters = imp.load_source('PartehInterpretParameters', \
                                            'py_modules/PartehInterpretParameters.py')
PartehTypes = imp.load_source('PartehTypes', 'py_modules/PartehTypes.py')
SyntheticBoundaries = imp.load_source('SyntheticBoundaries','py_modules/SyntheticBoundaries.py')
CDLParse = imp.load_source('CDLParse','py_modules/CDLParse.py')
F90ParamParse = imp.load_source('F90ParamParse','py_modules/F90ParamParse.py')

from PartehInterpretParameters import load_xml
from CDLParse import CDLParseDims, CDLParseParam, cdl_param_type
from F90ParamParse import f90_param_type, GetParamsInFile, GetPFTParmFileSymbols, MakeListUnique

f90_fates_integrators_obj_name = 'bld/FatesIntegratorsMod.o'
f90_fates_parteh_params_obj_name = 'bld/PRTParametersMod.o'
f90_fates_partehwrap_obj_name  = 'bld/FatesPARTEHWrapMod.o'
f90_fates_lossfluxes_obj_name  = 'bld/PRTLossFluxesMod.o'
f90_fates_parteh_generic_obj_name = 'bld/PRTGenericMod.o'
f90_fates_unitwrap_obj_name = 'bld/UnitWrapMod.o'
f90_fates_parteh_callom_obj_name = 'bld/PRTAllometricCarbonMod.o'
f90_fates_parteh_cnpallom_obj_name = 'bld/PRTAllometricCNPMod.o'
f90_fates_cohortwrap_obj_name = 'bld/FatesCohortWrapMod.o'
f90_fates_allom_obj_name = 'bld/FatesAllometryMod.o'

# -----------------------------------------------------------------------------------
#
# We may be calling fortran, if so, we need to initialize the modules
# This includes building the library objects, calling those objects
# and possibly allocating memory in those objects.  The fortran libraries
# and functions are held inside globally defined objects fates_f90_obj
#
# -----------------------------------------------------------------------------------

# Define the F90 objects
# These must be loaded according to the module dependency order
# Note that these calls instantiate the modules

f90_fates_integrators_obj    = ctypes.CDLL(f90_fates_integrators_obj_name,mode=ctypes.RTLD_GLOBAL)
f90_fates_parteh_params_obj  = ctypes.CDLL(f90_fates_parteh_params_obj_name,mode=ctypes.RTLD_GLOBAL)
f90_fates_unitwrap_obj        = ctypes.CDLL(f90_fates_unitwrap_obj_name,mode=ctypes.RTLD_GLOBAL)
f90_fates_parteh_generic_obj = ctypes.CDLL(f90_fates_parteh_generic_obj_name,mode=ctypes.RTLD_GLOBAL)
f90_fates_allom_obj          = ctypes.CDLL(f90_fates_allom_obj_name,mode=ctypes.RTLD_GLOBAL)
f90_fates_parteh_callom_obj  = ctypes.CDLL(f90_fates_parteh_callom_obj_name,mode=ctypes.RTLD_GLOBAL)
f90_fates_lossfluxes_obj     = ctypes.CDLL(f90_fates_lossfluxes_obj_name,mode=ctypes.RTLD_GLOBAL)
f90_fates_parteh_cnpallom_obj = ctypes.CDLL(f90_fates_parteh_cnpallom_obj_name,mode=ctypes.RTLD_GLOBAL)
f90_fates_partehwrap_obj = ctypes.CDLL(f90_fates_partehwrap_obj_name,mode=ctypes.RTLD_GLOBAL)
f90_fates_cohortwrap_obj = ctypes.CDLL(f90_fates_cohortwrap_obj_name,mode=ctypes.RTLD_GLOBAL)





# ========================================================================================
# ========================================================================================
#                                        Main
# ========================================================================================
# ========================================================================================


def main():

    # First check to make sure python 2.7 is being used
    version = platform.python_version()
    verlist = version.split('.')

    if( not ((verlist[0] == '2') & (verlist[1] ==  '7') & (int(verlist[2])>=15) )  ):
        print("The PARTEH driver must be run with python 2.7")
        print(" with tertiary version >=15.")
        print(" your version is {}".format(version))
        print(" exiting...")
        sys.exit(2)

    parser = argparse.ArgumentParser(description='Parse command line arguments to this script.')
    parser.add_argument('--xml-file', dest='xml_file', type=str, \
                        help="The path to the XML file controling this simulation.",required=True)
    args = parser.parse_args()

    xml_file = args.xml_file

    # This loads the dictionaries of, and lists of objects that
    # define the variables, parameters and forms that govern the
    # system of equations and solution
    [time_control, fates_cdl_file, driver_params, boundary_method,use_pfts] = load_xml(xml_file)

    num_plants = len(use_pfts)
    

    # -------------------------------------------------------------------------------------
    # Check through the fortran Code we are coupling with, determine the list of parameters
    # that we need.
    # -------------------------------------------------------------------------------------

    var_list = GetParamsInFile('../../parteh/PRTParametersMod.F90')


    # Now look through EDPftvarcon.F90 to determine the variable name in file
    # that is associated with the variable pointer

    var_list = GetPFTParmFileSymbols(var_list,'../../parteh/PRTParamsFATESMod.F90')
    

    # This variable is not added to the list we send to fortran, this
    # is only for the initial condition on the python side
    var_list.append(f90_param_type('hgt_min','fates_recruit_hgt_min',False))
    


    # -------------------------------------------------------------
    # We can now cross reference our list of parameters against
    # the parameter file. This will create a new list of parameters
    # however in the form of a dictionary. This dictionary of
    # entries is accessible by its symbol name, and will also
    # read in and store the actual parameter values from the file.
    # -------------------------------------------------------------

    dims = CDLParseDims(fates_cdl_file)

    fates_params = {}
    for elem in var_list:
        fates_params[elem.var_sym] = CDLParseParam(fates_cdl_file,cdl_param_type(elem.var_name,elem.in_f90),dims)
    print('Finished loading PFT parameters')

    
    
    num_pfts   = dims['fates_pft']
    num_organs = dims['fates_prt_organs']

    # Initialize the PARTEH instance
    iret=f90_fates_partehwrap_obj.__fatespartehwrapmod_MOD_spmappyset()

    # Allocate the PFT and ORGAN arrays  (leaf+root+sap+store+structure+repro = 6)

    WrapPFTAllocArbitrary([val for key,val in dims.iteritems()])

    
    # Set the phenology type
    phen_type = []
    for iplnt in range(num_plants):

        ipft = use_pfts[iplnt]
        evergreen        = np.int(fates_params['evergreen'].data[ipft])
        cold_deciduous   = np.int(fates_params['season_decid'].data[ipft])
        stress_deciduous = np.int(fates_params['stress_decid'].data[ipft])
        if(evergreen==1):
            if(cold_deciduous==1):
                print("Poorly defined phenology mode 0")
                exit(2)
            if(stress_deciduous==1):
                print("Poorly defined phenology mode 1")
                exit(2)
            phen_type.append(1)
        elif(cold_deciduous==1):
            if(evergreen==1):
                print("Poorly defined phenology mode 2")
                exit(2)
            if(stress_deciduous==1):
                print("Poorly defined phenology mode 3")
                exit(2)
            phen_type.append(2)
        elif(stress_deciduous==1):
            if(evergreen==1):
                print("Poorly defined phenology mode 4")
                exit(2)
            if(cold_deciduous==1):
                print("Poorly defined phenology mode 5")
                exit(2)
            phen_type.append(3)
        else:
            print("Unknown phenology mode ? {} {} {}".format(evergreen,cold_deciduous,stress_deciduous))
            exit(2)


    # -------------------------------------------------------------------------
    # Loop through all parameters in the "fates_params"
    # dictionary, send their data to the FORTRAN code
    # ------------------------------------------------------------------------

    # Loop through parameters
    for parm_key, parm_obj in fates_params.iteritems():

        # Loop through their dimensions
        # 2D case
        if(parm_obj.in_f90):
            if(parm_obj.ndims>1):
                for idx0 in range(parm_obj.dim_sizelist[0]):
                    for idx1 in range(parm_obj.dim_sizelist[1]):
                        iret = f90_fates_unitwrap_obj.__prtparamsgeneric_MOD_prtparamspyset(byref(c_double(parm_obj.data[idx0,idx1])), \
                                                                                            byref(c_int(0)), \
                                                                                            byref(c_int(idx1+1)), \
                                                                                            byref(c_int(idx0+1)), \
                                                                                            c_char_p(parm_obj.symbol), \
                                                                                            c_long(len(parm_obj.symbol )))

            else:
                idx1=0
                for idx0 in range(parm_obj.dim_sizelist[0]):
                    iret = f90_fates_unitwrap_obj.__prtparamsgeneric_MOD_prtparamspyset(byref(c_double(parm_obj.data[idx0])), \
                                                                                        byref(c_int(0)), \
                                                                                        byref(c_int(idx0+1)), \
                                                                                        byref(c_int(idx1+1)), \
                                                                                        c_char_p(parm_obj.symbol), \
                                                                                        c_long(len(parm_obj.symbol )))



    # Allocate the cohort array (We create on cohort per PFT)
    iret=f90_fates_cohortwrap_obj.__fatescohortwrapmod_MOD_cohortinitalloc(byref(c_int(num_plants)))

   
    
    for iplnt in range(num_plants):
        ipft = use_pfts[iplnt]
        hgt_min = np.float(fates_params['hgt_min'].data[ipft])
        init_canopy_trim = 1.0
        iret=f90_fates_cohortwrap_obj.__fatescohortwrapmod_MOD_cohortmodeset(byref(c_int(ipft)), \
                                                                             byref(c_int(int(driver_params['parteh_model'].param_vals[ipft]))))
        iret=f90_fates_cohortwrap_obj.__fatescohortwrapmod_MOD_cohortpyset(byref(c_int(ipft)), \
                                                                           byref(c_double(hgt_min)), \
                                                                           byref(c_double(init_canopy_trim)))


    # Initialize diagnostics
    diagnostics = []
    for iplnt in range(num_plants):
        ipft = use_pfts[iplnt]
        diagnostics.append(PartehTypes.diagnostics_type())


    # --------------------------------------------------------------------------------
    # Time Initialization
    # --------------------------------------------------------------------------------
    time_control.ResetTime()

    # --------------------------------------------------------------------------------
    # Time integration (outer) loop
    # --------------------------------------------------------------------------------
    while (time_control.sim_complete != True):

        print('Simulating Date: {}'.format(time_control.datetime.item()))

        # Start the integration substep loop
        endtime = time_control.datetime+np.timedelta64(int(time_control.dt_fullstep),'s')

        for iplnt in range(num_plants):

            ipft = use_pfts[iplnt]
            # Generate the boundary condition for the current time-step
            # ---------------------------------------------------------------------------

            # First lets query this pft-cohort and return a smattering of indices

            leaf_area  = c_double(0.0)
            agb        = c_double(0.0)
            crown_area = c_double(0.0)
            dbh        = c_double(0.0)
            target_leaf_c = c_double(-9.9)
            leaf_c     = c_double(0.0)
            fnrt_c     = c_double(0.0)
            sapw_c     = c_double(0.0)
            store_c    = c_double(0.0)
            struct_c   = c_double(0.0)
            repro_c    = c_double(0.0)
            root_c_exudate = c_double(0.0)
            growth_resp    = c_double(0.0)
            leaf_cturn     = c_double(0.0)
            fnrt_cturn     = c_double(0.0)
            sapw_cturn     = c_double(0.0)
            store_cturn    = c_double(0.0)
            struct_cturn   = c_double(0.0)

            leaf_n     = c_double(0.0)
            fnrt_n     = c_double(0.0)
            sapw_n     = c_double(0.0)
            store_n    = c_double(0.0)
            struct_n   = c_double(0.0)
            repro_n    = c_double(0.0)
            root_n_exudate = c_double(0.0)
            leaf_nturn     = c_double(0.0)
            fnrt_nturn     = c_double(0.0)
            sapw_nturn     = c_double(0.0)
            store_nturn    = c_double(0.0)
            struct_nturn   = c_double(0.0)

            leaf_p     = c_double(0.0)
            fnrt_p     = c_double(0.0)
            sapw_p     = c_double(0.0)
            store_p    = c_double(0.0)
            struct_p   = c_double(0.0)
            repro_p    = c_double(0.0)
            root_p_exudate = c_double(0.0)
            leaf_pturn     = c_double(0.0)
            fnrt_pturn     = c_double(0.0)
            sapw_pturn     = c_double(0.0)
            store_pturn    = c_double(0.0)
            struct_pturn   = c_double(0.0)

            iret=f90_fates_cohortwrap_obj.__fatescohortwrapmod_MOD_wrapqueryvars(byref(c_int(ipft)), \
                                                                                 byref(leaf_area), \
                                                                                 byref(crown_area), \
                                                                                 byref(agb), \
                                                                                 byref(store_c),\
                                                                                 byref(target_leaf_c))



            doy = time_control.datetime.astype(object).timetuple().tm_yday



            # Call phenology module, if no leaves... then npp should be zero...
            flush_c,drop_frac_c,leaf_status = SyntheticBoundaries.DeciduousPhenology(doy, \
                                                                                     target_leaf_c.value, \
                                                                                     store_c.value, phen_type[iplnt])



            if(boundary_method=="DailyCFromCArea"):

                presc_npp_p1     = driver_params['fates_prescribed_npp_p1'].param_vals[iplnt]

                net_daily_c = SyntheticBoundaries.DailyCFromCArea(presc_npp_p1, \
                                                                  crown_area.value, \
                                                                  phen_type[iplnt], \
                                                                  leaf_status)
                net_daily_n = 0.0
                net_daily_p = 0.0
                r_maint_demand = 0.0


            elif(boundary_method=="DailyCNPFromCArea"):

                presc_npp_p1   = driver_params['fates_prescribed_npp_p1'].param_vals[iplnt]
                presc_nflux_p1 = driver_params['fates_prescribed_nflux_p1'].param_vals[iplnt]
                presc_pflux_p1 = driver_params['fates_prescribed_pflux_p1'].param_vals[iplnt]

                net_daily_c, net_daily_n, net_daily_p = SyntheticBoundaries.DailyCNPFromCArea(presc_npp_p1, \
                                                                                              presc_nflux_p1, \
                                                                                              presc_pflux_p1, \
                                                                                              crown_area.value, \
                                                                                              phen_type[iplnt], \
                                                                                              leaf_status)
                r_maint_demand = 0.0


            elif(boundary_method=="DailyCNPFromStorageSinWaveNoMaint"):

                presc_npp_amp  = driver_params['fates_prescribed_npp_amp'].param_vals[iplnt]
                presc_npp_p1   = driver_params['fates_prescribed_npp_p1'].param_vals[iplnt]
                presc_nflux_p1 = driver_params['fates_prescribed_nflux_p1'].param_vals[iplnt]
                presc_pflux_p1 = driver_params['fates_prescribed_pflux_p1'].param_vals[iplnt]

                net_daily_c, net_daily_n, net_daily_p = SyntheticBoundaries.DailyCNPFromStorageSinWave(doy,\
                                                                                 store_c.value,\
                                                                                 presc_npp_p1, \
                                                                                 presc_nflux_p1, \
                                                                                 presc_pflux_p1, \
                                                                                 crown_area.value, \
                                                                                 presc_npp_amp, \
                                                                                 phen_type[iplnt], \
                                                                                 leaf_status )
                r_maint_demand = 0.0

            else:
                print("An unknown boundary method was specified\n")
                print("type: {} ? ... quitting.".format(boundary_method))
                exit()






            # This function will pass in all boundary conditions, some will be dummy arguments
            init_canopy_trim = 1.0
            iret=f90_fates_cohortwrap_obj.__fatescohortwrapmod_MOD_wrapdailyprt(byref(c_int(ipft)), \
                                                                                byref(c_double(net_daily_c)), \
                                                                                byref(c_double(init_canopy_trim)), \
                                                                                byref(c_double(flush_c)), \
                                                                                byref(c_double(drop_frac_c)), \
                                                                                byref(c_int(leaf_status)), \
                                                                                byref(c_double(net_daily_n)), \
                                                                                byref(c_double(net_daily_p)), \
                                                                                byref(c_double(r_maint_demand)))



            # This function will retrieve diagnostics
            iret=f90_fates_cohortwrap_obj.__fatescohortwrapmod_MOD_wrapquerydiagnostics(byref(c_int(ipft)),  \
                                                                                        byref(dbh),     \
                                                                                        byref(leaf_c),  \
                                                                                        byref(fnrt_c),  \
                                                                                        byref(sapw_c),  \
                                                                                        byref(store_c), \
                                                                                        byref(struct_c), \
                                                                                        byref(repro_c), \
                                                                                        byref(leaf_cturn),  \
                                                                                        byref(fnrt_cturn),  \
                                                                                        byref(sapw_cturn),  \
                                                                                        byref(store_cturn), \
                                                                                        byref(struct_cturn), \
                                                                                        byref(leaf_n),  \
                                                                                        byref(fnrt_n),  \
                                                                                        byref(sapw_n),  \
                                                                                        byref(store_n), \
                                                                                        byref(struct_n), \
                                                                                        byref(repro_n), \
                                                                                        byref(leaf_nturn),  \
                                                                                        byref(fnrt_nturn),  \
                                                                                        byref(sapw_nturn),  \
                                                                                        byref(store_nturn), \
                                                                                        byref(struct_nturn), \
                                                                                        byref(leaf_p),  \
                                                                                        byref(fnrt_p),  \
                                                                                        byref(sapw_p),  \
                                                                                        byref(store_p), \
                                                                                        byref(struct_p), \
                                                                                        byref(repro_p), \
                                                                                        byref(leaf_pturn),  \
                                                                                        byref(fnrt_pturn),  \
                                                                                        byref(sapw_pturn),  \
                                                                                        byref(store_pturn), \
                                                                                        byref(struct_pturn), \
                                                                                        byref(crown_area), \
                                                                                        byref(root_c_exudate), \
                                                                                        byref(root_n_exudate), \
                                                                                        byref(root_p_exudate), \
                                                                                        byref(growth_resp))


            diagnostics[iplnt].dates.append(time_control.datetime.astype(datetime))
            diagnostics[iplnt].dbh.append(dbh.value)
            diagnostics[iplnt].leaf_c.append(leaf_c.value)
            diagnostics[iplnt].fnrt_c.append(fnrt_c.value)
            diagnostics[iplnt].sapw_c.append(sapw_c.value)
            diagnostics[iplnt].store_c.append(store_c.value)
            diagnostics[iplnt].struct_c.append(struct_c.value)
            diagnostics[iplnt].repro_c.append(repro_c.value)
            diagnostics[iplnt].leaf_cturn.append(leaf_cturn.value)
            diagnostics[iplnt].fnrt_cturn.append(fnrt_cturn.value)
            diagnostics[iplnt].sapw_cturn.append(sapw_cturn.value)
            diagnostics[iplnt].store_cturn.append(store_cturn.value)
            diagnostics[iplnt].struct_cturn.append(struct_cturn.value)
            diagnostics[iplnt].dailyc.append(net_daily_c)
            diagnostics[iplnt].crown_area.append(crown_area.value)

            diagnostics[iplnt].growth_resp.append(growth_resp.value)

            diagnostics[iplnt].leaf_n.append(leaf_n.value)
            diagnostics[iplnt].fnrt_n.append(fnrt_n.value)
            diagnostics[iplnt].sapw_n.append(sapw_n.value)
            diagnostics[iplnt].store_n.append(store_n.value)
            diagnostics[iplnt].struct_n.append(struct_n.value)
            diagnostics[iplnt].repro_n.append(repro_n.value)
            diagnostics[iplnt].leaf_nturn.append(leaf_nturn.value)
            diagnostics[iplnt].fnrt_nturn.append(fnrt_nturn.value)
            diagnostics[iplnt].sapw_nturn.append(sapw_nturn.value)
            diagnostics[iplnt].store_nturn.append(store_nturn.value)
            diagnostics[iplnt].struct_nturn.append(struct_nturn.value)

            diagnostics[iplnt].leaf_p.append(leaf_p.value)
            diagnostics[iplnt].fnrt_p.append(fnrt_p.value)
            diagnostics[iplnt].sapw_p.append(sapw_p.value)
            diagnostics[iplnt].store_p.append(store_p.value)
            diagnostics[iplnt].struct_p.append(struct_p.value)
            diagnostics[iplnt].repro_p.append(repro_p.value)
            diagnostics[iplnt].leaf_pturn.append(leaf_pturn.value)
            diagnostics[iplnt].fnrt_pturn.append(fnrt_pturn.value)
            diagnostics[iplnt].sapw_pturn.append(sapw_pturn.value)
            diagnostics[iplnt].store_pturn.append(store_pturn.value)
            diagnostics[iplnt].struct_pturn.append(struct_pturn.value)

            diagnostics[iplnt].root_c_exudate.append(root_c_exudate.value)
            diagnostics[iplnt].root_n_exudate.append(root_n_exudate.value)
            diagnostics[iplnt].root_p_exudate.append(root_p_exudate.value)


        # We don't have a fancy time integrator so we simply update with
        # a full step

        time_control.UpdateTime()

        # ---------------------------------------------------------------------------
        # Timestep complete, check the time
        # ---------------------------------------------------------------------------
        #        time_control.CheckFullStepTime(endtime)


#    fig0, ax = plt.subplots()
#    for ipft in range(parameters.num_pfts):
#        ax.plot_date(diagnostics[0].dates,diagnostics[0].dbh)
#        ax.set_xlim(diagnostics[0].dates[0],diagnostics[0].dates[-1])

#    plt.show()
#    code.interact(local=locals())


    linestyles  = ['-','-.','--','-',':','-.','--',':','-','-.','--',':' ]




    fig1, ((ax1, ax2, ax3, ax4), (ax5, ax6, ax7, ax8)) = plt.subplots(2, 4 , sharex='col') #, sharey='row')
    fig1.set_size_inches(12, 6)
    for iplnt in range(num_plants):
        ipft = use_pfts[iplnt]
        ax1.plot_date(diagnostics[iplnt].dates,diagnostics[iplnt].struct_c,linestyles[iplnt],label='{}'.format(iplnt))
    ax1.set_title('Structural\n Carbon')
    ax1.legend(loc='upper left')
    ax1.set_ylabel('[kg C]')
    ax1.grid(True)

    for iplnt in range(num_plants):
        ax2.plot_date(diagnostics[iplnt].dates,diagnostics[iplnt].leaf_c,linestyles[iplnt])
    ax2.set_title('Leaf\n Carbon')
    ax2.grid(True)

    for iplnt in range(num_plants):
        ax3.plot_date(diagnostics[iplnt].dates,diagnostics[iplnt].fnrt_c,linestyles[iplnt])
    ax3.set_title('Fineroot\n Carbon')
    ax3.grid(True)

    for iplnt in range(num_plants):
        ax4.plot_date(diagnostics[iplnt].dates,diagnostics[iplnt].sapw_c,linestyles[iplnt])
    ax4.set_title('Sapwood\n Carbon')
    ax4.set_ylabel('[kg C]')
    ax4.grid(True)

    for iplnt in range(num_plants):
        ax5.plot_date(diagnostics[iplnt].dates,diagnostics[iplnt].store_c,linestyles[iplnt])
    ax5.set_title('Storage\n Carbon')
    ax5.set_xlabel('Year')
    ax5.grid(True)

    for iplnt in range(num_plants):
        ax6.plot_date(diagnostics[iplnt].dates,diagnostics[iplnt].repro_c,linestyles[iplnt])
    ax6.set_title('Integrated\n Reproductive\n Carbon')
    ax6.set_xlabel('Year')
    ax6.grid(True)

    for iplnt in range(num_plants):
        ax7.plot_date(diagnostics[iplnt].dates,np.cumsum(diagnostics[iplnt].root_c_exudate),linestyles[iplnt])
    ax7.set_title('Integrated\n Exudated\n Carbon')
    ax7.set_xlabel('Year')
    ax7.grid(True)

    for iplnt in range(num_plants):
        ax8.plot_date(diagnostics[iplnt].dates,np.cumsum(diagnostics[iplnt].growth_resp),linestyles[iplnt])
    ax8.set_title('Integrated\n Growth\n Respiration')
    ax8.set_xlabel('Year')
    ax8.grid(True)





    plt.tight_layout()

    #  Plant proportions
    #  ---------------------------------------------------------------------------------
    fig2, ( (ax1,ax2),(ax3,ax4) ) = plt.subplots(2,2)
    fig2.set_size_inches(7, 6)
    for iplnt in range(num_plants):
        ipft = use_pfts[iplnt]
        ax1.plot_date(diagnostics[iplnt].dates,diagnostics[iplnt].dbh,linestyles[iplnt],label='{}'.format(iplnt))
    ax1.set_xlabel('Date')
    ax1.set_title('DBH [cm]')
    ax1.legend(loc='upper left')
    ax1.grid(True)

    for iplnt in range(num_plants):
        ax2.plot_date(diagnostics[iplnt].dates,diagnostics[iplnt].crown_area,linestyles[iplnt])
    ax2.set_xlabel('Date')
    ax2.set_title('Crown Area [m2]')
    ax2.grid(True)

    for iplnt in range(num_plants):
        ax3.plot(diagnostics[iplnt].dbh,1000.0*np.array(diagnostics[iplnt].dailyc))

    ax3.set_xlabel('DBH [cm]')
    ax3.set_title('Daily Carbon Gain [g]')
    ax3.grid(True)

    for iplnt in range(num_plants):
        ax4.plot(diagnostics[iplnt].dbh,diagnostics[iplnt].crown_area)
    ax4.set_xlabel('DBH [cm]')
    ax4.set_title('Crown Area [m2]')
    ax4.grid(True)







    plt.tight_layout()


    # Error (bias)
    # ---------------------------------------------------------------------------------

    fig4 = plt.figure()
    for iplnt in range(num_plants):
        ipft = use_pfts[iplnt]
        total_plant_carbon0 = np.array(diagnostics[iplnt].struct_c[0]) + \
                              np.array(diagnostics[iplnt].leaf_c[0])   + \
                              np.array(diagnostics[iplnt].fnrt_c[0])   + \
                              np.array(diagnostics[iplnt].sapw_c[0])   + \
                              np.array(diagnostics[iplnt].store_c[0])  + \
                              np.array(diagnostics[iplnt].repro_c[0])

        total_plant_carbon = np.array(diagnostics[iplnt].struct_c) + \
                             np.array(diagnostics[iplnt].leaf_c)   + \
                             np.array(diagnostics[iplnt].fnrt_c)   + \
                             np.array(diagnostics[iplnt].sapw_c)   + \
                             np.array(diagnostics[iplnt].store_c)  + \
                             np.array(diagnostics[iplnt].repro_c)

        integrated_plant_turnover = np.cumsum(diagnostics[iplnt].struct_cturn) + \
                                    np.cumsum(diagnostics[iplnt].leaf_cturn) +  \
                                    np.cumsum(diagnostics[iplnt].fnrt_cturn) +  \
                                    np.cumsum(diagnostics[iplnt].sapw_cturn) +  \
                                    np.cumsum(diagnostics[iplnt].store_cturn)


        plt.plot(np.cumsum(diagnostics[iplnt].dailyc), \
                 (np.cumsum(diagnostics[iplnt].dailyc) - \
                            (total_plant_carbon + \
                             integrated_plant_turnover - \
                             total_plant_carbon0 ) ) / total_plant_carbon )

    plt.xlabel('Integrated Daily Carbon Gain [kg]')
    plt.ylabel('Integrated Bias [kg]')
    plt.grid(True)

    # Plot out the input fluxes

    fig5= plt.figure()
    for iplnt in range(num_plants):
        ipft = use_pfts[iplnt]
        plt.plot_date(diagnostics[iplnt].dates,diagnostics[iplnt].dailyc,linestyles[iplnt],label='{}'.format(iplnt))

    plt.xlabel('Date')
    plt.ylabel('Daily Carbon Flux')
    plt.grid(True)
    plt.legend(loc='upper left')


    # Special Focus plots for a PFT of interest

    figs = {}
    for iplnt in range(num_plants):
        ipft = use_pfts[iplnt]
        figs[iplnt], (ax1, ax2, ax3) = plt.subplots(1, 3)

        figs[iplnt].set_size_inches(8, 4)
        ax1.stackplot(np.cumsum(diagnostics[iplnt].dailyc), \
                      np.array(diagnostics[iplnt].struct_c)+np.cumsum(diagnostics[iplnt].struct_cturn), \
                      np.array(diagnostics[iplnt].leaf_c)+np.cumsum(diagnostics[iplnt].leaf_cturn), \
                      np.array(diagnostics[iplnt].fnrt_c)+np.cumsum(diagnostics[iplnt].fnrt_cturn), \
                      np.array(diagnostics[iplnt].sapw_c)+np.cumsum(diagnostics[iplnt].sapw_cturn), \
                      np.array(diagnostics[iplnt].store_c)+np.cumsum(diagnostics[iplnt].store_cturn), \
                      np.array(diagnostics[iplnt].repro_c), \
                      labels = ["Struct","Leaf","FRoot","Sapw","Storage","Repro"])
        ax1.set_title('Allocated Mass \nby Pool [kg]')
        ax1.grid(True)

        ax2.stackplot(np.cumsum(diagnostics[iplnt].dailyc), \
                      np.cumsum(diagnostics[iplnt].struct_cturn), \
                      np.cumsum(diagnostics[iplnt].leaf_cturn), \
                      np.cumsum(diagnostics[iplnt].fnrt_cturn),  \
                      np.cumsum(diagnostics[iplnt].sapw_cturn), \
                      np.cumsum(diagnostics[iplnt].store_cturn), \
                      np.array(diagnostics[iplnt].repro_c), \
                      labels = ["Struct","Leaf","FRoot","Sapw","Storage","Repro"] )
        ax2.legend(loc=2)
        ax2.grid(True)
        ax2.set_xlabel('Integrated Daily\n Carbon Gain [kg]')
        ax2.set_title('Integrated Turnover\n by Pool [kg]')


        #code.interact(local=locals())
        npp_leaf = np.array(diagnostics[iplnt].leaf_c[1:]) - \
                   np.array(diagnostics[iplnt].leaf_c[0:-1]) + \
                   np.array(diagnostics[iplnt].leaf_cturn[1:])
        npp_fnrt = np.array(diagnostics[iplnt].fnrt_c[1:]) - \
                   np.array(diagnostics[iplnt].fnrt_c[0:-1]) + \
                   np.array(diagnostics[iplnt].fnrt_cturn[1:])
        npp_sapw = np.array(diagnostics[iplnt].sapw_c[1:]) - \
                   np.array(diagnostics[iplnt].sapw_c[0:-1]) + \
                   np.array(diagnostics[iplnt].sapw_cturn[1:])
        npp_store = np.array(diagnostics[iplnt].store_c[1:]) - \
                    np.array(diagnostics[iplnt].store_c[0:-1]) + \
                    np.array(diagnostics[iplnt].store_cturn[1:])
        npp_struct = np.array(diagnostics[iplnt].struct_c[1:]) - \
                     np.array(diagnostics[iplnt].struct_c[0:-1]) + \
                     np.array(diagnostics[iplnt].struct_cturn[1:])
        npp_repro = np.array(diagnostics[iplnt].repro_c[1:]) - \
                    np.array(diagnostics[iplnt].repro_c[0:-1])

        ax3.stackplot(np.cumsum(diagnostics[iplnt].dailyc[1:]), \
                      npp_struct, npp_leaf, npp_fnrt, npp_sapw, npp_store,  npp_repro)

        ax3.grid(True)
        ax3.set_title('Daily NPP \nby Pool [kg]')

        plt.figtext(0.1,0.05,"Plant: {}".format(iplnt),bbox={'facecolor':'red', 'alpha':0.5, 'pad':10}, fontsize=15)


        plt.tight_layout()

    plt.show()

    print('\nSimulation Complete \nThank You Come Again')
    #exit(0)



def WrapPFTAllocArbitrary(*args):

    nargs = len(args[0])

    if(nargs==1):
        iret=f90_fates_unitwrap_obj.__prtparamsgeneric_MOD_prtparamsalloc(byref(c_int(args[0][0])))
    elif(nargs==2):
        iret=f90_fates_unitwrap_obj.__prtparamsgeneric_MOD_prtparamsalloc(byref(c_int(args[0][0])), byref(c_int(args[0][1])))
    elif(nargs==3):
        iret=f90_fates_unitwrap_obj.__prtparamsgeneric_MOD_prtparamsalloc(byref(c_int(args[0][0])), byref(c_int(args[0][1])), byref(c_int(args[0][2])))
    elif(nargs==4):
        iret=f90_fates_unitwrap_obj.__prtparamsgeneric_MOD_prtparamsalloc(byref(c_int(args[0][0])), byref(c_int(args[0][1])), byref(c_int(args[0][2])), \
                                                                       byref(c_int(args[0][3])))
    elif(nargs==5):
        iret=f90_fates_unitwrap_obj.__prtparamsgeneric_MOD_prtparamsalloc(byref(c_int(args[0][0])), byref(c_int(args[0][1])), byref(c_int(args[0][2])), \
                                                                       byref(c_int(args[0][3])), byref(c_int(args[0][4])))
    elif(nargs==6):
        iret=f90_fates_unitwrap_obj.__prtparamsgeneric_MOD_prtparamsalloc(byref(c_int(args[0][0])), byref(c_int(args[0][1])), byref(c_int(args[0][2])), \
                                                                       byref(c_int(args[0][3])), byref(c_int(args[0][4])), byref(c_int(args[0][5])))
    elif(nargs==7):
        iret=f90_fates_unitwrap_obj.__prtparamsgeneric_MOD_prtparamsalloc(byref(c_int(args[0][0])), byref(c_int(args[0][1])), byref(c_int(args[0][2])), \
                                                                       byref(c_int(args[0][3])), byref(c_int(args[0][4])), byref(c_int(args[0][5])), \
                                                                       byref(c_int(args[0][6])))
    elif(nargs==8):
        iret=f90_fates_unitwrap_obj.__prtparamsgeneric_MOD_prtparamsalloc(byref(c_int(args[0][0])), byref(c_int(args[0][1])), byref(c_int(args[0][2])), \
                                                                       byref(c_int(args[0][3])), byref(c_int(args[0][4])), byref(c_int(args[0][5])), \
                                                                       byref(c_int(args[0][6])), byref(c_int(args[0][7])))

    elif(nargs==9):
        iret=f90_fates_unitwrap_obj.__prtparamsgeneric_MOD_prtparamsalloc(byref(c_int(args[0][0])), byref(c_int(args[0][1])), byref(c_int(args[0][2])), \
                                                                       byref(c_int(args[0][3])), byref(c_int(args[0][4])), byref(c_int(args[0][5])), \
                                                                       byref(c_int(args[0][6])), byref(c_int(args[0][7])), byref(c_int(args[0][8])))
    elif(nargs==10):
        iret=f90_fates_unitwrap_obj.__prtparamsgeneric_MOD_prtparamsalloc(byref(c_int(args[0][0])), byref(c_int(args[0][1])), byref(c_int(args[0][2])), \
                                                                       byref(c_int(args[0][3])), byref(c_int(args[0][4])), byref(c_int(args[0][5])), \
                                                                       byref(c_int(args[0][6])), byref(c_int(args[0][7])), byref(c_int(args[0][8])), \
                                                                       byref(c_int(args[0][9])))
    elif(nargs==11):
        iret=f90_fates_unitwrap_obj.__prtparamsgeneric_MOD_prtparamsalloc(byref(c_int(args[0][0])), byref(c_int(args[0][1])), byref(c_int(args[0][2])), \
                                                                       byref(c_int(args[0][3])), byref(c_int(args[0][4])), byref(c_int(args[0][5])), \
                                                                       byref(c_int(args[0][6])), byref(c_int(args[0][7])), byref(c_int(args[0][8])), \
                                                                       byref(c_int(args[0][9])), byref(c_int(args[0][10])))


    else:
        print('So many dimensions...')
        print('add more clauses')
        exit(2)

# =======================================================================================



def usage():
     print('')
     print('=======================================================================')
     print('')
     print(' python PartehDriver.py --help --cdlfile=<path-to-file>')
     print('')
     print(' This is a driver script for PARTEH')
     print(' (Plant Allocation and Reactive Transport Extensible Hypotheses)')
     print(' Only 1 option is currently relevent, and that is a path to the ')
     print(' XML file that controls this simulation. ')
     print('')
     print(' Arguments:')
     print('')
     print(' -h --help ')
     print(' print this help message')
     print('')
     print(' --xmlfile = <path-to-file>')
     print(' the relative or full file path to the xml file that controls')
     print(' this simulation.')
     print('')



     
# =======================================================================================
# This is the actual call to main

if __name__ == "__main__":
    main()
