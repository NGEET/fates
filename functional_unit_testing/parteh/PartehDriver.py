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
#from matplotlib.backends.backend_pdf import PdfPages
import platform
import numpy as np
import os
import sys
import getopt
import code  # For development: code.interact(local=locals())
import time
import imp
import ctypes
from ctypes import *
from operator import add

PartehInterpretParameters = imp.load_source('PartehInterpretParameters', \
                                            'py_modules/PartehInterpretParameters.py')
PartehTypes = imp.load_source('PartehTypes', 'py_modules/PartehTypes.py')
SyntheticBoundaries = imp.load_source('SyntheticBoundaries','py_modules/SyntheticBoundaries.py')

from PartehInterpretParameters import load_xml

f90_fates_wrap_obj_name = 'bld/FatesWrapMod.o'
f90_fates_integrators_obj_name = 'bld/FatesIntegratorsMod.o'
f90_fates_partehwrap_obj_name  = 'bld/FatesPARTEHWrapMod.o'
f90_fates_lossfluxes_obj_name  = 'bld/PRTLossFluxesMod.o'
f90_fates_parteh_generic_obj_name = 'bld/PRTGenericMod.o'
f90_fates_pftwrap_obj_name = 'bld/FatesPFTWrapMod.o'
f90_fates_parteh_callom_obj_name = 'bld/PRTAllometricCarbonMod.o'
f90_fates_parteh_cnpallom_obj_name = 'bld/PRTAllometricCNPMod.o'
f90_fates_cohortwrap_obj_name = 'bld/FatesCohortWrapMod.o'
f90_fates_allom_obj_name = 'bld/FatesAllometryMod.o'

# =======================================================================================
# Some Global Parmaeters

## The name of the xml file containing site data (should not change)
xml_file = ''


# ========================================================================================
# ========================================================================================
#                                        Main
# ========================================================================================
# ========================================================================================


def main(argv):

    # First check to make sure python 2.7 is being used
    version = platform.python_version()
    verlist = version.split('.')

    if( not ((verlist[0] == '2') & (verlist[1] ==  '7') & (int(verlist[2])>=15) )  ):
        print("The PARTEH driver mus be run with python 2.7")
        print(" with tertiary version >=15.")
        print(" your version is {}".format(version))
        print(" exiting...")
        sys.exit(2)


    # Retrieve the name and path to the xml control file
    # from the input arguments
    xml_file = interp_args(argv)

    # Initialize the time structure
    time_control = PartehTypes.timetype()

    # Initialize the parameter structure
    parameters = PartehTypes.param_type()

    # This loads the dictionaries of, and lists of objects that
    # define the variables, parameters and forms that govern the
    # system of equations and solution
    load_xml(xml_file,time_control,parameters)

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
    f90_fates_wrap_obj           = ctypes.CDLL(f90_fates_wrap_obj_name,mode=ctypes.RTLD_GLOBAL)
    f90_fates_integrators_obj    = ctypes.CDLL(f90_fates_integrators_obj_name,mode=ctypes.RTLD_GLOBAL)
    f90_fates_pftwrap_obj        = ctypes.CDLL(f90_fates_pftwrap_obj_name,mode=ctypes.RTLD_GLOBAL)
    f90_fates_parteh_generic_obj = ctypes.CDLL(f90_fates_parteh_generic_obj_name,mode=ctypes.RTLD_GLOBAL)
    f90_fates_allom_obj          = ctypes.CDLL(f90_fates_allom_obj_name,mode=ctypes.RTLD_GLOBAL)
    f90_fates_parteh_callom_obj  = ctypes.CDLL(f90_fates_parteh_callom_obj_name,mode=ctypes.RTLD_GLOBAL)
    f90_fates_lossfluxes_obj     = ctypes.CDLL(f90_fates_lossfluxes_obj_name,mode=ctypes.RTLD_GLOBAL)
    f90_fates_parteh_cnpallom_obj = ctypes.CDLL(f90_fates_parteh_cnpallom_obj_name,mode=ctypes.RTLD_GLOBAL)
    f90_fates_partehwrap_obj = ctypes.CDLL(f90_fates_partehwrap_obj_name,mode=ctypes.RTLD_GLOBAL)
    f90_fates_cohortwrap_obj = ctypes.CDLL(f90_fates_cohortwrap_obj_name,mode=ctypes.RTLD_GLOBAL)
    
    # Initialize the PARTEH instance
    iret=f90_fates_partehwrap_obj.__fatespartehwrapmod_MOD_spmappyset() #byref(c_int(parameters.prt_model)))
    
    # Allocate the PFT and ORGAN arrays  (leaf+root+sap+store+structure+repro = 6)
    max_num_organs = 6
    iret=f90_fates_pftwrap_obj.__edpftvarcon_MOD_edpftvarconalloc(byref(c_int(parameters.num_pfts)), \
                                                                  byref(c_int(max_num_organs)))
    
    # Set the phenology type
    phen_type = []
    for pft_idx,pft_obj in enumerate(parameters.parteh_pfts):

        evergreen        = np.int(parameters.parteh_pfts[pft_idx].param_dic['fates_phen_evergreen'][0])
        cold_deciduous   = np.int(parameters.parteh_pfts[pft_idx].param_dic['fates_phen_season_decid'][0])
        stress_deciduous = np.int(parameters.parteh_pfts[pft_idx].param_dic['fates_phen_stress_decid'][0])
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



    # Loop through each pft and pft's parameters and pass them to the fortran object
    # Also, some parameters may be arrays (like organ number)
    for pft_idx,pft_obj in enumerate(parameters.parteh_pfts):

        for par_idx, par_key in enumerate(pft_obj.param_dic.iterkeys()):
            pval = pft_obj.param_dic[par_key]
            print("{} {} {}".format(par_idx,par_key,pval))

            # The dictionary of parameters is populated with lists of floats, even
            # scalars are single entry lists
            
            if( len(pval)==1 ):
                iret = f90_fates_pftwrap_obj.__edpftvarcon_MOD_edpftvarconpyset(byref(c_int(pft_idx+1)), \
                                                                                byref(c_int(0)), \
                                                                                byref(c_double(pval[0])), \
                                                                                c_char_p(par_key.strip()), \
                                                                                c_long(len(par_key.strip())))
            else:
                for i2d in range(len(pval)):
                    iret = f90_fates_pftwrap_obj.__edpftvarcon_MOD_edpftvarconpyset(byref(c_int(pft_idx+1)), \
                                                                                    byref(c_int(i2d+1)), \
                                                                                    byref(c_double(pval[i2d])), \
                                                                                    c_char_p(par_key.strip()), \
                                                                                    c_long(len(par_key.strip())))
            
    # Allocate the cohort array (We create on cohort per PFT)
    iret=f90_fates_cohortwrap_obj.__fatescohortwrapmod_MOD_cohortinitalloc(byref(c_int(parameters.num_pfts)))

    for pft_idx, pft_obj in enumerate(parameters.parteh_pfts):
        hgt_min = pft_obj.param_dic['fates_recruit_hgt_min']
        init_canopy_trim = 1.0
        iret=f90_fates_cohortwrap_obj.__fatescohortwrapmod_MOD_cohortpyset(byref(c_int(pft_idx+1)), \
                                                                           byref(c_double(hgt_min[0])), \
                                                                           byref(c_double(init_canopy_trim))) #, \
            #                                                                           byref(c_int(parameters.prt_model)))
    
    # Initialize diagnostics
    diagnostics = [] 
    for pft_idx, pft_obj in enumerate(parameters.parteh_pfts):
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

        for pft_idx, pft_obj in enumerate(parameters.parteh_pfts):

            
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

            iret=f90_fates_cohortwrap_obj.__fatescohortwrapmod_MOD_wrapqueryvars(byref(c_int(pft_idx+1)), \
                                                                                 byref(leaf_area), \
                                                                                 byref(crown_area), \
                                                                                 byref(agb), \
                                                                                 byref(store_c),\
                                                                                 byref(target_leaf_c))

           

            doy = time_control.datetime.astype(object).timetuple().tm_yday



            # Call phenology module, if no leaves... then npp should be zero...
            flush_c,drop_frac_c,leaf_status = SyntheticBoundaries.DeciduousPhenology(doy, \
                                                                                     target_leaf_c.value, \
                                                                                     store_c.value, phen_type[pft_idx])

            if(parameters.boundary_method=="DailyCFromCArea"):
                
                presc_npp_p1     = parameters.boundary_pfts[pft_idx].param_dic['fates_prescribed_npp_p1']

                net_daily_c = SyntheticBoundaries.DailyCFromCArea(presc_npp_p1, \
                                                                  crown_area.value, \
                                                                  phen_type[pft_idx], \
                                                                  leaf_status)
                net_daily_n = 0.0
                net_daily_p = 0.0
                r_maint_demand = 0.0
                

            elif(parameters.boundary_method=="DailyCNPFromCArea"):

                presc_npp_p1   = parameters.boundary_pfts[pft_idx].param_dic['fates_prescribed_npp_p1']
                presc_nflux_p1 = parameters.boundary_pfts[pft_idx].param_dic['fates_prescribed_nflux_p1']
                presc_pflux_p1 = parameters.boundary_pfts[pft_idx].param_dic['fates_prescribed_pflux_p1']

                net_daily_c, net_daily_n, net_daily_p = SyntheticBoundaries.DailyCNPFromCArea(presc_npp_p1, \
                                                                                              presc_nflux_p1, \
                                                                                              presc_pflux_p1, \
                                                                                              crown_area.value, \
                                                                                              phen_type[pft_idx], \
                                                                                              leaf_status)
                r_maint_demand = 0.0


            elif(parameters.boundary_method=="DailyCNPFromStorageSinWaveNoMaint"):

                presc_npp_amp  = parameters.boundary_pfts[pft_idx].param_dic['fates_prescribed_npp_amp']
                presc_npp_p1   = parameters.boundary_pfts[pft_idx].param_dic['fates_prescribed_npp_p1']
                presc_nflux_p1 = parameters.boundary_pfts[pft_idx].param_dic['fates_prescribed_nflux_p1']
                presc_pflux_p1 = parameters.boundary_pfts[pft_idx].param_dic['fates_prescribed_pflux_p1']
                
                
               

                net_daily_c, net_daily_n, net_daily_p = SyntheticBoundaries.DailyCNPFromStorageSinWave(doy,\
                                                                                 store_c.value,\
                                                                                 presc_npp_p1, \
                                                                                 presc_nflux_p1, \
                                                                                 presc_pflux_p1, \
                                                                                 crown_area.value, \
                                                                                 presc_npp_amp, \
                                                                                 phen_type[pft_idx], \
                                                                                 leaf_status )
                r_maint_demand = 0.0

            else:
                print("An unknown boundary method was specified\n")
                print("type: {} ? ... quitting.".format(parameters.boundary_method))
                exit()


            
            


            # This function will pass in all boundary conditions, some will be dummy arguments
            init_canopy_trim = 1.0
            iret=f90_fates_cohortwrap_obj.__fatescohortwrapmod_MOD_wrapdailyprt(byref(c_int(pft_idx+1)), \
                                                                                byref(c_double(net_daily_c)), \
                                                                                byref(c_double(init_canopy_trim)), \
                                                                                byref(c_double(flush_c)), \
                                                                                byref(c_double(drop_frac_c)), \
                                                                                byref(c_int(leaf_status)), \
                                                                                byref(c_double(net_daily_n)), \
                                                                                byref(c_double(net_daily_p)), \
                                                                                byref(c_double(r_maint_demand)))
                                                                                

            
            # This function will retrieve diagnostics
            iret=f90_fates_cohortwrap_obj.__fatescohortwrapmod_MOD_wrapquerydiagnostics(byref(c_int(pft_idx+1)),  \
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


            diagnostics[pft_idx].dates.append(time_control.datetime.astype(datetime))
            diagnostics[pft_idx].dbh.append(dbh.value)
            diagnostics[pft_idx].leaf_c.append(leaf_c.value)
            diagnostics[pft_idx].fnrt_c.append(fnrt_c.value)
            diagnostics[pft_idx].sapw_c.append(sapw_c.value)
            diagnostics[pft_idx].store_c.append(store_c.value)
            diagnostics[pft_idx].struct_c.append(struct_c.value) 
            diagnostics[pft_idx].repro_c.append(repro_c.value)
            diagnostics[pft_idx].leaf_cturn.append(leaf_cturn.value)
            diagnostics[pft_idx].fnrt_cturn.append(fnrt_cturn.value)
            diagnostics[pft_idx].sapw_cturn.append(sapw_cturn.value)
            diagnostics[pft_idx].store_cturn.append(store_cturn.value)
            diagnostics[pft_idx].struct_cturn.append(struct_cturn.value)
            diagnostics[pft_idx].dailyc.append(net_daily_c)
            diagnostics[pft_idx].crown_area.append(crown_area.value)
            
            diagnostics[pft_idx].growth_resp.append(growth_resp.value)

            diagnostics[pft_idx].leaf_n.append(leaf_n.value)
            diagnostics[pft_idx].fnrt_n.append(fnrt_n.value)
            diagnostics[pft_idx].sapw_n.append(sapw_n.value)
            diagnostics[pft_idx].store_n.append(store_n.value)
            diagnostics[pft_idx].struct_n.append(struct_n.value) 
            diagnostics[pft_idx].repro_n.append(repro_n.value)
            diagnostics[pft_idx].leaf_nturn.append(leaf_nturn.value)
            diagnostics[pft_idx].fnrt_nturn.append(fnrt_nturn.value)
            diagnostics[pft_idx].sapw_nturn.append(sapw_nturn.value)
            diagnostics[pft_idx].store_nturn.append(store_nturn.value)
            diagnostics[pft_idx].struct_nturn.append(struct_nturn.value)

            diagnostics[pft_idx].leaf_p.append(leaf_p.value)
            diagnostics[pft_idx].fnrt_p.append(fnrt_p.value)
            diagnostics[pft_idx].sapw_p.append(sapw_p.value)
            diagnostics[pft_idx].store_p.append(store_p.value)
            diagnostics[pft_idx].struct_p.append(struct_p.value) 
            diagnostics[pft_idx].repro_p.append(repro_p.value)
            diagnostics[pft_idx].leaf_pturn.append(leaf_pturn.value)
            diagnostics[pft_idx].fnrt_pturn.append(fnrt_pturn.value)
            diagnostics[pft_idx].sapw_pturn.append(sapw_pturn.value)
            diagnostics[pft_idx].store_pturn.append(store_pturn.value)
            diagnostics[pft_idx].struct_pturn.append(struct_pturn.value)

            diagnostics[pft_idx].root_c_exudate.append(root_c_exudate.value)
            diagnostics[pft_idx].root_n_exudate.append(root_n_exudate.value)
            diagnostics[pft_idx].root_p_exudate.append(root_p_exudate.value)


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
    for ipft in range(parameters.num_pfts):
        ax1.plot_date(diagnostics[ipft].dates,diagnostics[ipft].struct_c,linestyles[ipft],label=parameters.parteh_pfts[ipft].name)
    ax1.set_title('Structural\n Carbon')
    ax1.legend(loc='upper left')
    ax1.set_ylabel('[kg C]')
    ax1.grid(True)

    for ipft in range(parameters.num_pfts):
        ax2.plot_date(diagnostics[ipft].dates,diagnostics[ipft].leaf_c,linestyles[ipft])
    ax2.set_title('Leaf\n Carbon')
    ax2.grid(True)

    for ipft in range(parameters.num_pfts):
        ax3.plot_date(diagnostics[ipft].dates,diagnostics[ipft].fnrt_c,linestyles[ipft])
    ax3.set_title('Fineroot\n Carbon')
    ax3.grid(True)

    for ipft in range(parameters.num_pfts):
        ax4.plot_date(diagnostics[ipft].dates,diagnostics[ipft].sapw_c,linestyles[ipft])
    ax4.set_title('Sapwood\n Carbon')
    ax4.set_ylabel('[kg C]')
    ax4.grid(True)

    for ipft in range(parameters.num_pfts):
        ax5.plot_date(diagnostics[ipft].dates,diagnostics[ipft].store_c,linestyles[ipft])
    ax5.set_title('Storage\n Carbon')
    ax5.set_xlabel('Year')
    ax5.grid(True)

    for ipft in range(parameters.num_pfts):
        ax6.plot_date(diagnostics[ipft].dates,diagnostics[ipft].repro_c,linestyles[ipft])
    ax6.set_title('Integrated\n Reproductive\n Carbon')
    ax6.set_xlabel('Year')
    ax6.grid(True)

    for ipft in range(parameters.num_pfts):
        ax7.plot_date(diagnostics[ipft].dates,np.cumsum(diagnostics[ipft].root_c_exudate),linestyles[ipft])
    ax7.set_title('Integrated\n Exudated\n Carbon')
    ax7.set_xlabel('Year')
    ax7.grid(True)

    for ipft in range(parameters.num_pfts):
        ax8.plot_date(diagnostics[ipft].dates,np.cumsum(diagnostics[ipft].growth_resp),linestyles[ipft])
    ax8.set_title('Integrated\n Growth\n Respiration')
    ax8.set_xlabel('Year')
    ax8.grid(True)


    


    plt.tight_layout()

    #  Plant proportions
    #  ---------------------------------------------------------------------------------
    fig2, ( (ax1,ax2),(ax3,ax4) ) = plt.subplots(2,2)
    fig2.set_size_inches(7, 6)
    for ipft in range(parameters.num_pfts):
        ax1.plot_date(diagnostics[ipft].dates,diagnostics[ipft].dbh,linestyles[ipft],label=parameters.parteh_pfts[ipft].name)
    ax1.set_xlabel('Date')
    ax1.set_title('DBH [cm]')
    ax1.legend(loc='upper left')
    ax1.grid(True)

    for ipft in range(parameters.num_pfts):
        ax2.plot_date(diagnostics[ipft].dates,diagnostics[ipft].crown_area,linestyles[ipft])
    ax2.set_xlabel('Date')
    ax2.set_title('Crown Area [m2]')
    ax2.grid(True)

    for ipft in range(parameters.num_pfts):
        ax3.plot(diagnostics[ipft].dbh,1000.0*np.array(diagnostics[ipft].dailyc))

    ax3.set_xlabel('DBH [cm]')
    ax3.set_title('Daily Carbon Gain [g]')
    ax3.grid(True)

    for ipft in range(parameters.num_pfts):
        ax4.plot(diagnostics[ipft].dbh,diagnostics[ipft].crown_area)
    ax4.set_xlabel('DBH [cm]')
    ax4.set_title('Crown Area [m2]')
    ax4.grid(True)

    


   
    

    plt.tight_layout()


    # Error (bias)
    # ---------------------------------------------------------------------------------
    
    fig4 = plt.figure()
    for ipft in range(parameters.num_pfts):
        
        total_plant_carbon0 = np.array(diagnostics[ipft].struct_c[0]) + \
                              np.array(diagnostics[ipft].leaf_c[0])   + \
                              np.array(diagnostics[ipft].fnrt_c[0])   + \
                              np.array(diagnostics[ipft].sapw_c[0])   + \
                              np.array(diagnostics[ipft].store_c[0])  + \
                              np.array(diagnostics[ipft].repro_c[0])
        
        total_plant_carbon = np.array(diagnostics[ipft].struct_c) + \
                             np.array(diagnostics[ipft].leaf_c)   + \
                             np.array(diagnostics[ipft].fnrt_c)   + \
                             np.array(diagnostics[ipft].sapw_c)   + \
                             np.array(diagnostics[ipft].store_c)  + \
                             np.array(diagnostics[ipft].repro_c)
        
        integrated_plant_turnover = np.cumsum(diagnostics[ipft].struct_cturn) + \
                                    np.cumsum(diagnostics[ipft].leaf_cturn) +  \
                                    np.cumsum(diagnostics[ipft].fnrt_cturn) +  \
                                    np.cumsum(diagnostics[ipft].sapw_cturn) +  \
                                    np.cumsum(diagnostics[ipft].store_cturn)


        plt.plot(np.cumsum(diagnostics[ipft].dailyc), \
                 (np.cumsum(diagnostics[ipft].dailyc) - \
                            (total_plant_carbon + \
                             integrated_plant_turnover - \
                             total_plant_carbon0 ) ) / total_plant_carbon )
                 
    plt.xlabel('Integrated Daily Carbon Gain [kg]')
    plt.ylabel('Integrated Bias [kg]')
    plt.grid(True)

    # Plot out the input fluxes

    fig5= plt.figure()
    for ipft in range(parameters.num_pfts):
        plt.plot_date(diagnostics[ipft].dates,diagnostics[ipft].dailyc,linestyles[ipft],label=parameters.parteh_pfts[ipft].name)

    plt.xlabel('Date')
    plt.ylabel('Daily Carbon Flux')
    plt.grid(True)
    plt.legend(loc='upper left')
    
        
    # Special Focus plots for a PFT of interest

    figs = {}
    for ipft in range(parameters.num_pfts):
        figs[ipft], (ax1, ax2, ax3) = plt.subplots(1, 3)
        
        figs[ipft].set_size_inches(8, 4)
        ax1.stackplot(np.cumsum(diagnostics[ipft].dailyc), \
                      np.array(diagnostics[ipft].struct_c)+np.cumsum(diagnostics[ipft].struct_cturn), \
                      np.array(diagnostics[ipft].leaf_c)+np.cumsum(diagnostics[ipft].leaf_cturn), \
                      np.array(diagnostics[ipft].fnrt_c)+np.cumsum(diagnostics[ipft].fnrt_cturn), \
                      np.array(diagnostics[ipft].sapw_c)+np.cumsum(diagnostics[ipft].sapw_cturn), \
                      np.array(diagnostics[ipft].store_c)+np.cumsum(diagnostics[ipft].store_cturn), \
                      np.array(diagnostics[ipft].repro_c), \
                      labels = ["Struct","Leaf","FRoot","Sapw","Storage","Repro"])
        ax1.set_title('Allocated Mass \nby Pool [kg]')
        ax1.grid(True)

        ax2.stackplot(np.cumsum(diagnostics[ipft].dailyc), \
                      np.cumsum(diagnostics[ipft].struct_cturn), \
                      np.cumsum(diagnostics[ipft].leaf_cturn), \
                      np.cumsum(diagnostics[ipft].fnrt_cturn),  \
                      np.cumsum(diagnostics[ipft].sapw_cturn), \
                      np.cumsum(diagnostics[ipft].store_cturn), \
                      np.array(diagnostics[ipft].repro_c), \
                      labels = ["Struct","Leaf","FRoot","Sapw","Storage","Repro"] )
        ax2.legend(loc=2)
        ax2.grid(True)
        ax2.set_xlabel('Integrated Daily\n Carbon Gain [kg]')
        ax2.set_title('Integrated Turnover\n by Pool [kg]')

    
        #code.interact(local=locals())
        npp_leaf = np.array(diagnostics[ipft].leaf_c[1:]) - \
                   np.array(diagnostics[ipft].leaf_c[0:-1]) + \
                   np.array(diagnostics[ipft].leaf_cturn[1:])
        npp_fnrt = np.array(diagnostics[ipft].fnrt_c[1:]) - \
                   np.array(diagnostics[ipft].fnrt_c[0:-1]) + \
                   np.array(diagnostics[ipft].fnrt_cturn[1:])
        npp_sapw = np.array(diagnostics[ipft].sapw_c[1:]) - \
                   np.array(diagnostics[ipft].sapw_c[0:-1]) + \
                   np.array(diagnostics[ipft].sapw_cturn[1:])
        npp_store = np.array(diagnostics[ipft].store_c[1:]) - \
                    np.array(diagnostics[ipft].store_c[0:-1]) + \
                    np.array(diagnostics[ipft].store_cturn[1:])
        npp_struct = np.array(diagnostics[ipft].struct_c[1:]) - \
                     np.array(diagnostics[ipft].struct_c[0:-1]) + \
                     np.array(diagnostics[ipft].struct_cturn[1:])    
        npp_repro = np.array(diagnostics[ipft].repro_c[1:]) - \
                    np.array(diagnostics[ipft].repro_c[0:-1])
    
        ax3.stackplot(np.cumsum(diagnostics[ipft].dailyc[1:]), \
                      npp_struct, npp_leaf, npp_fnrt, npp_sapw, npp_store,  npp_repro)

        ax3.grid(True)
        ax3.set_title('Daily NPP \nby Pool [kg]')

        plt.figtext(0.1,0.05,"PFT: {}".format(ipft+1),bbox={'facecolor':'red', 'alpha':0.5, 'pad':10}, fontsize=15)


        plt.tight_layout()


    plt.show()

    print('\nSimulation Complete \nThank You Come Again')
    #exit(0)


            
# =======================================================================================



def usage():
     print('')
     print('=======================================================================')
     print('')
     print(' python PartehDriver.py --help --xmlfile=<path-to-file>')
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

def interp_args(argv):

    argv.pop(0)  # The script itself is the first argument, forget it

    ## File path to the xml control card
    xmlfile = ''

    try:
        opts, args = getopt.getopt(argv, 'h',["help","xmlfile="])

    except getopt.GetoptError as err:
        print('Argument error, see usage')
        usage()
        sys.exit(2)


    if(len(opts)==0):
        print('\n\n')
        print('No arguments were specified')
        print('Exiting, see Usage below')
        print('\n\n')
        usage()
        sys.exit(0)

    for o, a in opts:
        if o in ("-h", "--help"):
            usage()
            sys.exit(0)
        elif o in ("--xmlfile"):
            xmlfile = a.strip()
            if(not os.path.isfile(xmlfile)):
                print('\n\n')
                print('The XML control file could not be found')
                print(' via argument --xmlfile')
                print(' xmlfile = ',xmlfile)
                print('\n\n')
                usage()
                sys.exit(0)
        else:
            assert False, "unhandled option"

    return(xmlfile)
# =======================================================================================
# This is the actual call to main
   
if __name__ == "__main__":
    main(sys.argv)
