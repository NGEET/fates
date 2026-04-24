# =======================================================================================
#
# For usage: $python LeafBiophysDriver.py --help
#
# This script runs unit tests on the leaf biophysics functions
#
#
# =======================================================================================

import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt
from datetime import datetime
import argparse
from matplotlib.backends.backend_pdf import PdfPages
import platform
import xml.etree.ElementTree as et
import numpy as np
import matplotlib
import os
import sys
import getopt
import code  # For development: code.interact(local=dict(globals(), **locals()))
import time
import importlib
import csv
import subprocess
import re
import CtypesLeafBiophys
import ctypes
from ctypes import *
from operator import add
sys.path.append('../shared/py_src')
from PyF90Utils import c8, ci, cchar, c8_arr, ci_arr, ccharnb

font = {'family' : 'sans-serif',
        'weight' : 'normal',
        'size'   : 11}

matplotlib.rc('font', **font)


# Global constants to use in all Leaf Biophysics unit testing
# =======================================================================================


# For debugging
dump_parameters = False

# Should we evaluate vcmax, jmax and kp actual?
do_evalvjkbytemp = False

# Do an analysis on convergence
do_test_citol = False

# Freezing point of water in Kelvin (at standard atmosphere)
tfrz_1atm = 273.15

# 25 degrees C in Kelvin (used because T25 functions)
leaf_tempk25 = tfrz_1atm + 25.0

# Daylight limitations can be imposed on Vcmax, a value of
# 1 means daylight length is at its maximum
dayl_factor_full = 1.0

# If Kumerathunga respiration is used, it requires moving averages
# of leaf temperature
t_growth_kum = -999
t_home_kum = -999

# Simple conversion, number of micro-moles in a mole
umol_per_mol = 1.e6
mol_per_umol = 1.e-6

# 1 standard atmosphere in [Pa]
can_press_1atm = 101325.0

# Atmospheric CO2 partial pressure [Pa] at 400 ppm
co2_ppress_400ppm = 0.0004*can_press_1atm

# Atmospheric O2 partial pressure [Pa] %29.5 of atmosphere
o2_ppress_209kppm = 0.2095*can_press_1atm

# 70% of atmospheric CO2 is a reasonablish guess for
# intercellular CO2 concentration during primary production
# We can use this to test the gross assimilation routines
# directly, without having to solve for the equilibrium
# intercellular CO2
ci_ppress_static = 0.7*co2_ppress_400ppm

# When there is hydrualic limitation on photosynthesis
# (via Vcmax reductions), then the btran factor is 1
btran_nolimit = 1.0

# Respiration scaler at canopy top
rdark_scaler_top = 1.0

# Nitrogen scaler at canopy top
nscaler_top = 1.0


# Create aliases for the ctype Fortran objects
# =======================================================================================

exec(open("CtypesLeafBiophys.py").read())


# Subroutines
# =======================================================================================


def GetJmaxKp25Top(vcmax25_top):

    # Calculate Jmax and Kp at the canopy top at 25C
    # they scale off of vcmax
    #
    # jmax25_top:  Canopy top maximum electron transport
    #              rate at 25C (umol electrons/m**2/s)
    #
    # kp25top      Canopy top initial slope of CO2 response
    #              curve (C4 plants) at 25C
    
    jmax25_top = 1.67   * vcmax25_top
    kp25_top   = 20000.  * vcmax25_top
    
    # q10 response of product limited psn.
    # co2_rcurve_islope = co2_rcurve_islope25 * 2._r8**((veg_tempk-(tfrz+25._r8))/10._r8)
    
    return jmax25_top, kp25_top



def main(argv):

    parser = argparse.ArgumentParser(description='Parse command line arguments to this script.')
    parser.add_argument('--fin', dest='xmlfile', type=str, help="XML control file, required.", required=True)
    parser.add_argument('--smoketest', action='store_true')
    args = parser.parse_args()

    # Load the xml control file
    # -----------------------------------------------------------------------------------
    xmlroot = et.parse(args.xmlfile).getroot()
    
    # We will allocate 1 token pft to hold data, we will change the values as needed
    numpft = 1

    # Allocating parameters
    print('Allocating parameter space for {} pfts'.format(numpft))
    iret = f90_alloc_leaf_param_sub(ci(numpft))

    
    
    # Push scalar parameters
    print('Pushing parameters from the xml file to the f90 lb_params datastructure')
    scalar_root = xmlroot.find('f90_params').find('scalar_dim')
    for param in scalar_root.iter('param'):
        iret = f90_set_leaf_param_sub(c8(float(param.text.split(',')[0])),ci(0),*ccharnb(param.attrib['name'].strip()))

    # Push pft parameters to fortran instantiations
    pft_root = xmlroot.find('f90_params').find('pft_dim')
    leaf_c3psn = []
    leaf_stomatal_intercept = []
    for param in pft_root.iter('param'):
        for pft in range(numpft):
            iret = f90_set_leaf_param_sub(c8(float(param.text.split(',')[pft])),ci(pft+1),*ccharnb(param.attrib['name'].strip()))
            if(param.attrib['name'].strip() == 'fates_leaf_c3psn'):
                leaf_c3psn.append(int(param.text.split(',')[pft]))
            if(param.attrib['name'].strip() == 'fates_leaf_stomatal_intercept'):
                leaf_stomatal_intercept.append(int(param.text.split(',')[pft]))
                
    # Read in non-fortran parameters from the xml file
    fates_leaf_vcmax25top    = []
    fates_stoich_nitr = []
    fates_leaf_slatop = []
    
    print('Reading non-fortran pft parameters')
    py_pft_root = xmlroot.find('py_params').find('pft_dim')
    for param in py_pft_root.iter('param'):
        for pft in range(numpft):
            if (param.attrib['name']=='fates_leaf_vcmax25top'):
                fates_leaf_vcmax25top.append(np.float64(param.text.split(',')[pft]))
            if (param.attrib['name']=='fates_leaf_slatop'):
                fates_leaf_slatop.append(np.float64(param.text.split(',')[pft]))
            if (param.attrib['name']=='fates_stoich_nitr'):
                fates_stoich_nitr.append(np.float64(param.text.split(',')[pft]))
                
    print('Reading non-fortran scalar parameters')
    py_scalar_root = xmlroot.find('py_params').find('scalar_dim')
    for param in py_scalar_root.iter('param'):
        if (param.attrib['name']=='fates_maintresp_leaf_model'):
            fates_maintresp_leaf_model = int(param.text.split(',')[0])


    # Axes to test:
    #   Temperature
    #   PAR
    #   Humidity
    #   Boundary Layer conductance
    #   BTRAN
    #     vcmax (via btran)
    #     stomatal intercept (via btran)
    #     stomatal slope 1   (via btran)
    #     stomatal slope 2   (via btran/medlyn)
    # Switches to test
    #   Stomtal Model (medlyn and BB)
    #   C3/C4
    #   fates_leaf_stomatal_btran_model (0,1,2,3)
    #   fates_leaf_agross_btran_model
    
    # Axes implicitly tested
    # daylength factors (scale vcmax, but covered by btran scaling)
    # 

    leaf_stomatal_btran_models = [0,1,2,3,4,5]   #4,5 medlyn only
    leaf_agross_btran_models = [0,1,2]
    
        
    # Leaf temperature ranges [C]
    leaf_tempc_min = -50.0
    leaf_tempc_max = 60.0
    leaf_tempc_n = 15
    leaf_tempc_vec = np.linspace(leaf_tempc_min,leaf_tempc_max,num=leaf_tempc_n)

    # Relative Humidity Ranges
    rh_max = 1.00
    rh_min = 0.001
    rh_n   = 15
    rh_vec = np.linspace(rh_min,rh_max,num=rh_n)
    
    # Absorbed PAR ranges [W/m2]
    par_abs_min = 0.0
    par_abs_max = 800
    par_abs_n  = 15
    par_abs_vec = np.linspace(par_abs_min,par_abs_max,num=par_abs_n)

    # Boundary Conductance ranges [umol/m2/s]
    gb_min =   50000.0            # Lower limit imposed by CLM/ELM 0.5 mol/m2/s
    gb_max = 8500000.0            # 50% larger than  Roughly largestthe largest values seen at BCI (which are 2.5mol/m2/s)
    gb_n  = 10
    gb_vec = np.linspace(gb_min,gb_max,num=gb_n)

    # btran ranges

    btran_n   = 15
    btran_min = -5
    btran_max = 0
    btran_vec = np.logspace(btran_min,btran_max,num=btran_n)

    # vcmax25top ranges
    vcmax25t_n = 10
    vcmax25t_min = 1
    vcmax25t_max = 250
    vcmax25t_vec = np.linspace(vcmax25t_min,vcmax25t_max,num=vcmax25t_n)


    # Set convergence tolerance
    ci_tol = 0.1
    
    # Initialize fortran return values
    # these are all temps and doubles
    vcmax_f      = c_double(-9.0)
    jmax_f       = c_double(-9.0)
    kp_f         = c_double(-9.0)
    agross_f     = c_double(-9.0)
    gstoma_f     = c_double(-9.0)
    anet_f       = c_double(-9.0)
    lmr_f        = c_double(-9.0)
    c13_f        = c_double(-9.0)
    ac_f         = c_double(-9.0)
    aj_f         = c_double(-9.0)
    ap_f         = c_double(-9.0)
    co2_interc_f = c_double(-9.0)
    veg_qs_f     = c_double(-9.0)
    veg_es_f     = c_double(-9.0)
    mm_kco2_f    = c_double(-9.0)
    mm_ko2_f     = c_double(-9.0)
    co2_cpoint_f = c_double(-9.0)
    qsdt_dummy_f = c_double(-9.0)
    esdt_dummy_f = c_double(-9.0)
    solve_iter_f = c_int(-9)
    gs0_f        = c_double(-9.0)
    gs1_f        = c_double(-9.0)
    gs2_f        = c_double(-9.0)

    pfails = 0
    ptests = 0
    printfail = True
    
    # unit conversion of W/m2 to umol photons/m^2/s
    wm2_to_umolm2s = 4.6
  
    print('\nStarting smoke tests with the following nested combinations:\n')

    print(' 2 conductance models, Medlyn and Ball-Berry')
    print(' 2 pathway (C3/C4) models')
    print(' {} stomatal btran options from {} to {}'.format(len(leaf_stomatal_btran_models),leaf_stomatal_btran_models[0],leaf_stomatal_btran_models[-1]))
    print(' {} agross btran options from {} to {}'.format(len(leaf_agross_btran_models),leaf_agross_btran_models[0],leaf_agross_btran_models[-1]))
    print(' {} leaf temperature values [C] from {} to {}'.format(leaf_tempc_n,leaf_tempc_min,leaf_tempc_max))
    print(' {} RH values [fraction] from {} to {}'.format(rh_n,rh_min,rh_max))
    print(' {} PAR Abs [W/m2] values from {} to {}'.format(par_abs_n,par_abs_min,par_abs_max))
    print(' {} BL conductance (gb) [umol/m2/s] values from {} to {}'.format(gb_n,gb_min,gb_max))
    print(' {} BTRAN values [fraction] from {} to {}'.format(btran_n,np.exp(btran_min),np.exp(btran_max)))
    print(' {} Vcmax 25 top values [umol/m2/s] from {} to {}'.format(vcmax25t_n,vcmax25t_min,vcmax25t_max))
    
    #print(' {} values from {} to {}'.format())

    ntests = 2*len(leaf_agross_btran_models)*vcmax25t_n*leaf_tempc_n*btran_n*gb_n*par_abs_n*rh_n*(4 + 6)

    # Report every 5%
    ntestmod = int(ntests/20)
    
    
    print('\nRunning a total of {} tests: \n'.format(ntests))
    time0 = time.process_time()
    # Switch between Medlyn and BB
    for ism in [1,2]:

        # Push the conductance choice to the fortran code
        iret = f90_set_leaf_param_sub(c8(float( ism  )),ci(0),*ccharnb('fates_leaf_stomatal_model'))
        
        if ism==1:            leaf_stomatal_btran_models = [0,1,2,3]
        else:
            leaf_stomatal_btran_models = [0,1,2,3,4,5]   #4,5 medlyn only
    
        for isb in leaf_stomatal_btran_models:

            iret = f90_set_leaf_param_sub(c8(float( isb  )),ci(1),*ccharnb('fates_leaf_stomatal_btran_model'))
            
            for iab in leaf_agross_btran_models:

                iret = f90_set_leaf_param_sub(c8(float( iab  )),ci(1),*ccharnb('fates_leaf_agross_btran_model'))

                for ic3 in [0,1]:

                    iret = f90_set_leaf_param_sub(c8(float( ic3  )),ci(1),*ccharnb('fates_leaf_c3psn'))

                    for vcmax25_top in vcmax25t_vec:
                    
                        jmax25_top,kp25_top =  GetJmaxKp25Top(vcmax25_top)
                    
                        for leaf_tempc in leaf_tempc_vec:
                        
                            leaf_tempk = leaf_tempc + tfrz_1atm
            
                            iret = f90_qsat_sub(c8(leaf_tempk),c8(can_press_1atm), \
                                                byref(veg_qs_f),byref(veg_es_f), \
                                                byref(qsdt_dummy_f),byref(esdt_dummy_f))

                            iret = f90_cangas_sub(c8(can_press_1atm), \
                                                  c8(o2_ppress_209kppm), \
                                                  c8(leaf_tempk), \
                                                  byref(mm_kco2_f), \
                                                  byref(mm_ko2_f), \
                                                  byref(co2_cpoint_f))
        
                            # Leaf Nitrogen Concentration at the top
                            lnc_top  = fates_stoich_nitr[0]/fates_leaf_slatop[0]

                            # Leaf Maintenance Respiration (temp and pft dependent)
                            if(fates_maintresp_leaf_model==1):
                                iret = f90_lmr_ryan_sub(c8(lnc_top),c8(nscaler_top), ci(1), c8(leaf_tempk), byref(lmr_f))
                            elif(fates_maintresp_leaf_model==2):
                                iret = f90_lmr_atkin_sub(c8(lnc_top),c8(rdark_scaler_top),c8(leaf_tempk),c8(atkin_mean_leaf_tempk),byref(lmr_f) )
                            else:
                                print('unknown leaf respiration model')
                                exit(1)
                                
                            for btran in btran_vec:
                                iret = f90_biophysrate_sub(ci(1), \
                                                           c8(vcmax25_top), c8(jmax25_top), c8(kp25_top), \
                                                           c8(nscaler_top), c8(leaf_tempk), c8(dayl_factor_full), \
                                                           c8(t_growth_kum), c8(t_home_kum), c8(btran), \
                                                           byref(vcmax_f), byref(jmax_f), byref(kp_f), byref(gs0_f), byref(gs1_f), byref(gs2_f))
            
                                for gb in gb_vec:
                                    for par_abs in par_abs_vec:
                                        par_abs_umol = par_abs*wm2_to_umolm2s
                                        for rh in rh_vec:
                                            vpress = rh * veg_es_f.value
                                            ptests = ptests + 1
                                            try:
                                                iret = f90_leaflayerphoto_sub(c8(par_abs_umol),  \
                                                                              ci(1),   \
                                                                              c8(vcmax_f.value),   \
                                                                              c8(jmax_f.value),    \
                                                                              c8(kp_f.value),      \
                                                                              c8(gs0_f.value), \
                                                                              c8(gs1_f.value), \
                                                                              c8(gs2_f.value), \
                                                                              c8(leaf_tempk), \
                                                                              c8(can_press_1atm), \
                                                                              c8(co2_ppress_400ppm), \
                                                                              c8(o2_ppress_209kppm), \
                                                                              c8(veg_es_f.value), \
                                                                              c8(gb), \
                                                                              c8(vpress), \
                                                                              c8(mm_kco2_f.value), \
                                                                              c8(mm_ko2_f.value), \
                                                                              c8(co2_cpoint_f.value), \
                                                                              c8(lmr_f.value), \
                                                                              c8(ci_tol), \
                                                                              byref(agross_f), \
                                                                              byref(gstoma_f), \
                                                                              byref(anet_f), \
                                                                              byref(c13_f), \
                                                                              byref(co2_interc_f), \
                                                                              byref(solve_iter_f) )
                                            except:
                                                pfails = pfails+1
                                                printfail=True
                                                
                                            if (np.mod(ptests,ntestmod)==0):
                                                print('Completed {} tests -- {} percent complete'.format(ptests,100*float(ptests)/float(ntests)))

                                            if (pfails>0 and np.mod(pfails,100)==0 and printfail):
                                                printfail=False
                                                print('\n{} fails so far\n'.format(pfails))
                                                

                                                
    print("\nCompleted Photosynthesis Smoke Test\n")
    print("\nElapsed Time [s]: {}\n".format(time.process_time() - time0))
    print("{} Failures out of {} Encountered; {}% of Tests\n".format(pfails,ptests,float(pfails)/float(ptests)))

    print('Deallocating parameter space')
    iret = f90_dealloc_leaf_param_sub()
    
    print('Functional Unit Testing Complete')
    exit(0)

    
# =======================================================================================
# This is the actual call to main

if __name__ == "__main__":
    main(sys.argv)
