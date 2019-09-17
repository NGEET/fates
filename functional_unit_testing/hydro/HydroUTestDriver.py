# =======================================================================================
#
# For usage: $python HydroUTestDriver.py --help
#
# This script runs unit tests on the hydraulics functions.
#
#
# =======================================================================================

import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt
from datetime import datetime
import argparse
#from matplotlib.backends.backend_pdf import PdfPages
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


CDLParse = imp.load_source('CDLParse','../shared/py_src/CDLParse.py')
F90ParamParse = imp.load_source('F90ParamParse','../shared/py_src/F90ParamParse.py')
PyF90Utils = imp.load_source('PyF90Utils','../shared/py_src/PyF90Utils.py')

from CDLParse import CDLParseDims, CDLParseParam, cdl_param_type
from F90ParamParse import f90_param_type, GetSymbolUsage, GetPFTParmFileSymbols, MakeListUnique
from PyF90Utils import c8, ci, cchar

# Load the fortran objects via CTYPES

f90_edparams_obj = ctypes.CDLL('bld/EDParamsHydroMod.o',mode=ctypes.RTLD_GLOBAL)
f90_constants_obj = ctypes.CDLL('bld/FatesConstantsMod.o',mode=ctypes.RTLD_GLOBAL)
f90_unitwrap_obj = ctypes.CDLL('bld/UnitWrapMod.o',mode=ctypes.RTLD_GLOBAL)
f90_hydrofuncs_obj = ctypes.CDLL('bld/FatesHydroUnitFunctionsMod.o',mode=ctypes.RTLD_GLOBAL)

# Alias the F90 functions, specify the return type
# -----------------------------------------------------------------------------------
psi_from_th = f90_hydrofuncs_obj.__fateshydrounitfunctionsmod_MOD_psi_from_th
psi_from_th.restype = c_double

dpsidth_from_th= f90_hydrofuncs_obj.__fateshydrounitfunctionsmod_MOD_dpsidth_from_th
dpsidth_from_th.restype = c_double

flc_from_psi = f90_hydrofuncs_obj.__fateshydrounitfunctionsmod_MOD_flc_from_psi
flc_from_psi.restype = c_double

dflcdpsi_from_psi= f90_hydrofuncs_obj.__fateshydrounitfunctionsmod_MOD_dflcdpsi_from_psi
dflcdpsi_from_psi.restype = c_double

solutepsi = f90_hydrofuncs_obj.__fateshydrounitfunctionsmod_MOD_solutepsi
pressurepsi = f90_hydrofuncs_obj.__fateshydrounitfunctionsmod_MOD_pressurepsi


# Some constants
rwcft  = [1.0,0.958,0.958,0.958]
rwccap = [1.0,0.947,0.947,0.947]
pm_leaf = 1
pm_stem = 2
pm_aroot = 4
pm_troot = 3

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

    # Read in the arguments
    # =======================================================================================

    parser = argparse.ArgumentParser(description='Parse command line arguments to this script.')
    parser.add_argument('--cdl-file', dest='cdlfile', type=str, \
                        help="Input CDL filename.  Required.", required=True)

    args = parser.parse_args()



    # -------------------------------------------------------------------------------------
    # Check through the fortran Code we are coupling with, determine the list of parameters
    # that we need.
    # The procedure GetSymbolUsage() returns a list of strings (non-unique)
    # -------------------------------------------------------------------------------------

    check_str = 'pft_p%'
    var_list0 = GetSymbolUsage('../../biogeophys/FatesHydroUnitFunctionsMod.F90',check_str)

    # This is the unique list of PFT parameters found in the salient Fortran code

    var_list = MakeListUnique(var_list0)

    # Now look through EDPftvarcon.F90 to determine the variable name in file
    # that is associated with the variable pointer

    var_list = GetPFTParmFileSymbols(var_list,'../../main/EDPftvarcon.F90')

    # -------------------------------------------------------------
    # We can now cross reference our list of parameters against
    # the parameter file. This will create a new list of parameters
    # however in the form of a dictionary. This dictionary of
    # entries is accessible by its symbol name, and will also
    # read in and store the actual parameter values from the file.
    # -------------------------------------------------------------

    dims = CDLParseDims(args.cdlfile)
    pftparms = {}
    for elem in var_list:
        pftparms[elem.var_sym] = CDLParseParam(args.cdlfile,cdl_param_type(elem.var_name),dims)
    print('Finished loading PFT parameters')

    num_pfts   = dims['fates_pft']

    scalarparms = {}
    scalarparms['hydr_psi0'] = CDLParseParam(args.cdlfile,cdl_param_type('fates_hydr_psi0'),dims)
    scalarparms['hydr_psicap'] = CDLParseParam(args.cdlfile,cdl_param_type('fates_hydr_psicap'),dims)




    # Allocate PFT arrays in the fortran objects
    iret=f90_unitwrap_obj.__edpftvarcon_MOD_edpftvarconalloc(ci(dims['fates_string_length']), \
                                                             ci(dims['fates_history_size_bins']), \
                                                             ci(dims['fates_NCWD']), \
                                                             ci(dims['fates_prt_organs']), \
                                                             ci(dims['fates_litterclass']), \
                                                             ci(dims['fates_history_height_bins']), \
                                                             ci(dims['fates_history_age_bins']), \
                                                             ci(dims['fates_hydr_organs']), \
                                                             ci(dims['fates_pft']), \
                                                             ci(dims['fates_variants']), \
                                                             ci(dims['fates_leafage_class']))



    # Set the PFT arrays
    for pft_key,pft_obj in pftparms.iteritems():
        for idim in range(np.int(np.prod(pft_obj.dim_sizelist))):
            if(pft_obj.ndims==1):
                idim1 = idim
                idim2 = 0
                rdata = pft_obj.data[idim]
                idata = np.int(pft_obj.data[idim])
            else:
                idim2 = np.mod(idim,num_pfts)
                idim1 = np.int(idim/num_pfts)
                rdata = pft_obj.data[idim1,idim2]
                idata = np.int(pft_obj.data[idim1,idim2])
            iret = f90_unitwrap_obj.__edpftvarcon_MOD_edpftvarconpyset(c8(rdata), \
                                                                       ci(idata), \
                                                                       ci(idim2+1), \
                                                                       ci(idim1+1), \
                                                                       c_char_p(pft_obj.symbol.strip()), \
                                                                       c_long(len(pft_obj.symbol.strip())))

    # Push the scalar params data to fortran

    iret = f90_edparams_obj.__edparamsmod_MOD_edparamspyset(c8(scalarparms['hydr_psi0'].data[0]), \
                                                            c_char_p(scalarparms['hydr_psi0'].symbol.strip()), \
                                                            c_long(len(scalarparms['hydr_psi0'].symbol.strip())))

    iret = f90_edparams_obj.__edparamsmod_MOD_edparamspyset(c8(scalarparms['hydr_psicap'].data[0]), \
                                                            c_char_p(scalarparms['hydr_psicap'].symbol.strip()), \
                                                            c_long(len(scalarparms['hydr_psicap'].symbol.strip())))


    # Initialize local objects in the unit test
    iret = f90_hydrofuncs_obj.__fateshydrounitfunctionsmod_MOD_initallocateplantmedia(ci(4))
    iret = f90_hydrofuncs_obj.__fateshydrounitfunctionsmod_MOD_setplantmediaparam(ci(1),c8(rwcft[0]),c8(rwccap[0]))
    iret = f90_hydrofuncs_obj.__fateshydrounitfunctionsmod_MOD_setplantmediaparam(ci(2),c8(rwcft[1]),c8(rwccap[1]))
    iret = f90_hydrofuncs_obj.__fateshydrounitfunctionsmod_MOD_setplantmediaparam(ci(3),c8(rwcft[2]),c8(rwccap[2]))
    iret = f90_hydrofuncs_obj.__fateshydrounitfunctionsmod_MOD_setplantmediaparam(ci(4),c8(rwcft[3]),c8(rwccap[3]))




    # Test 1 For a set of thetas, calculate psi for each pm.
    # ===================================================================================

    pft1 = 1
    pft2 = 2

    npts = 1000

    min_theta = 0.01
    max_theta = 0.99
    min_leaf_theta = pftparms['hydr_resid_node'].data[pm_leaf-1,pft1-1]
    max_leaf_theta = pftparms['hydr_thetas_node'].data[pm_leaf-1,pft1-1]
    min_stem_theta = pftparms['hydr_resid_node'].data[pm_stem-1,pft1-1]
    max_stem_theta = pftparms['hydr_thetas_node'].data[pm_stem-1,pft1-1]
    min_troot_theta = pftparms['hydr_resid_node'].data[pm_troot-1,pft1-1]
    max_troot_theta = pftparms['hydr_thetas_node'].data[pm_troot-1,pft1-1]
    min_aroot_theta = pftparms['hydr_resid_node'].data[pm_aroot-1,pft1-1]
    max_aroot_theta = pftparms['hydr_thetas_node'].data[pm_aroot-1,pft1-1]

    min_leaf_theta2 = pftparms['hydr_resid_node'].data[pm_leaf-1,pft2-1]
    min_stem_theta2 = pftparms['hydr_resid_node'].data[pm_stem-1,pft2-1]
    min_stem_theta2 = pftparms['hydr_resid_node'].data[pm_troot-1,pft2-1]
    min_stem_theta2 = pftparms['hydr_resid_node'].data[pm_aroot-1,pft2-1]



    theta       = np.linspace(min_theta, max_theta, num=npts-1)
    leaf_theta  = np.linspace(min_leaf_theta,max_leaf_theta, num=npts)
    stem_theta  = np.linspace(min_stem_theta,max_stem_theta, num=npts)
    troot_theta = np.linspace(min_troot_theta,max_troot_theta, num=npts)
    aroot_theta = np.linspace(min_aroot_theta,max_aroot_theta, num=npts)

    leaf_theta2  = np.linspace(min_leaf_theta2,max_leaf_theta, num=npts)

    # "full" range psi's
    leaf_fpsi  = np.zeros(shape=np.shape(theta),dtype=np.float64)
    stem_fpsi  = np.zeros(shape=np.shape(theta),dtype=np.float64)
    troot_fpsi = np.zeros(shape=np.shape(theta),dtype=np.float64)
    aroot_fpsi = np.zeros(shape=np.shape(theta),dtype=np.float64)

    # "constrained" range psi's
    leaf_cpsi  = np.full(shape=np.shape(leaf_theta),dtype=np.float64,fill_value=np.nan)
    leaf_cpsi2 = np.full(shape=np.shape(leaf_theta2),dtype=np.float64,fill_value=np.nan)
    stem_cpsi  = np.zeros(shape=np.shape(stem_theta),dtype=np.float64)
    troot_cpsi = np.zeros(shape=np.shape(troot_theta),dtype=np.float64)
    aroot_cpsi = np.zeros(shape=np.shape(aroot_theta),dtype=np.float64)


    mpl.rcParams.update({'font.size': 15})


    # Initialize the return variable
    cpsi = c_double(0)


    # Find PSI for each theta
    for i,th in enumerate(theta):
        leaf_fpsi[i] = psi_from_th(ci(pft1), ci(pm_leaf), c8(th))
        stem_fpsi[i] = psi_from_th(ci(pft1), ci(pm_stem), c8(th))
        troot_fpsi[i] = psi_from_th(ci(pft1), ci(pm_troot), c8(th))
        aroot_fpsi[i] = psi_from_th(ci(pft1), ci(pm_aroot), c8(th))

    for i,th in enumerate(leaf_theta):
        leaf_cpsi[i] = psi_from_th(ci(pft1), ci(pm_leaf), c8(th))

    for i,th in enumerate(leaf_theta2):
        leaf_cpsi2[i] = psi_from_th(ci(pft2),ci(pm_leaf),c8(th))

    for i,th in enumerate(stem_theta):
        stem_cpsi[i] = psi_from_th(ci(pft1), ci(pm_stem), c8(th))

    for i,th in enumerate(troot_theta):
        troot_cpsi[i] = psi_from_th(ci(pft1), ci(pm_troot), c8(th))

    for i,th in enumerate(aroot_theta):
        aroot_cpsi[i] = psi_from_th(ci(pft1), ci(pm_aroot), c8(th))



    fig0, (ax1,ax2) = plt.subplots(2)
    ax1.plot(theta,leaf_fpsi,label='Leaf')
    ax1.plot(theta,stem_fpsi,label='Stem')
    ax1.plot(theta,troot_fpsi,label='Troot')
    ax1.plot(theta,aroot_fpsi,label='Aroot')
    ax1.grid(True)
    ax1.set_ylabel('Psi')
    ax1.set_ylim((-10,10))
    ax1.set_xlim((0,1))
    ax1.set_title('Unconstrained Range')

    ax2.plot(leaf_theta,leaf_cpsi,label=r'Leaf')
    ax2.plot(stem_theta,stem_cpsi,label='Stem')
    ax2.plot(troot_theta,troot_cpsi,label='Troot')
    ax2.plot(aroot_theta,aroot_cpsi,label='Aroot')
    ax2.grid(True)
    ax2.set_xlabel('Theta')
    ax2.set_ylabel('Psi')
    ax2.set_ylim((-10,0))
    ax2.set_xlim((0,1))
    ax2.set_title('Constrained Residual and Sat')
    ax2.legend(loc='lower right')
    plt.tight_layout()


    # Compare how resid and sat effect the curve
    fig1, ax1 = plt.subplots(1)
    ax1.plot(leaf_theta,leaf_cpsi,label='pft1')
    ax1.plot(leaf_theta2,leaf_cpsi2,label='pft2')
    ax1.grid(True)
    ax1.set_ylabel('Psi')
    ax1.set_xlabel('theta')
#    ax1.set_ylim((-30,0))
    ax1.set_xlim((0,1))
    ax1.set_title('Comparing PFT 1 and PFT 2')


    # Lets look at the solute curve, is that the one that is breaking?

    resid = pftparms['hydr_resid_node'].data[pm_leaf-1,pft2-1]
    thetas = pftparms['hydr_thetas_node'].data[pm_leaf-1,pft2-1]
    pinot = pftparms['hydr_pinot_node'].data[pm_leaf-1,pft2-1]
    epsil =  pftparms['hydr_epsil_node'].data[pm_leaf-1,pft2-1]
    sol_psi = np.full(shape=np.shape(leaf_theta2),dtype=np.float64,fill_value=np.nan)
    ela_psi = np.full(shape=np.shape(leaf_theta2),dtype=np.float64,fill_value=np.nan)
    sol_psi2 = np.full(shape=np.shape(leaf_theta2),dtype=np.float64,fill_value=np.nan)
    ela_psi2 = np.full(shape=np.shape(leaf_theta2),dtype=np.float64,fill_value=np.nan)

    print(resid,thetas,pinot,epsil)


    for i,th in enumerate(leaf_theta2):
        sol_psi[i] = pinot*thetas*(rwcft[pm_leaf-1] - resid) / (th - thetas*resid)
        ela_psi[i] = epsil * (th - thetas*rwcft[pm_leaf-1]) / (thetas*(rwcft[pm_leaf-1]-resid)) - pinot

        iret = solutepsi(ci(pft2),ci(pm_leaf),c8(th),byref(cpsi))
        sol_psi2[i] = cpsi.value
        iret = pressurepsi(ci(pft2),ci(pm_leaf),c8(th),byref(cpsi))
        ela_psi2[i] = cpsi.value


    fig11,ax1 = plt.subplots(1)
    ax1.plot(leaf_theta2,sol_psi,label='solute psi')
    ax1.plot(leaf_theta2,ela_psi,label='press psi')
    ax1.plot(leaf_theta2,sol_psi2,label='f(sol)')
    ax1.plot(leaf_theta2,ela_psi2,label='f(ela)')
    ax1.legend(loc='lower right')
    ax1.grid(True)



    # Derivative Check on PSI
    # -----------------------------------------------------------------------------------

    leaf_dpsidth  = np.full(shape=np.shape(leaf_theta),dtype=np.float64,fill_value=np.nan)
    leaf_dpsidthc = np.full(shape=np.shape(leaf_theta),dtype=np.float64,fill_value=np.nan)
    for i,th in enumerate(leaf_theta[1:-2]):
        leaf_dpsidth[i] = dpsidth_from_th(ci(pft1), ci(pm_leaf), c8(th))
        leaf_dpsidthc[i] = (leaf_cpsi[i+1]-leaf_cpsi[i-1])/(leaf_theta[i+1]-leaf_theta[i-1])



    # Derivative Check
    stem_dpsidth  = np.full(shape=np.shape(stem_theta),dtype=np.float64,fill_value=np.nan)
    stem_dpsidthc = np.full(shape=np.shape(stem_theta),dtype=np.float64,fill_value=np.nan)
    for i,th in enumerate(stem_theta[1:-2]):
        stem_dpsidth[i] = dpsidth_from_th(ci(pft1), ci(pm_stem), c8(th))
        stem_dpsidthc[i] = (stem_cpsi[i+1]-stem_cpsi[i-1])/(stem_theta[i+1]-stem_theta[i-1])


    # Derivative Check
    troot_dpsidth  = np.full(shape=np.shape(troot_theta),dtype=np.float64,fill_value=np.nan)
    troot_dpsidthc = np.full(shape=np.shape(troot_theta),dtype=np.float64,fill_value=np.nan)
    for i,th in enumerate(troot_theta[1:-2]):
        troot_dpsidth[i] = dpsidth_from_th(ci(pft1), ci(pm_troot), c8(th))
        troot_dpsidthc[i] = (troot_cpsi[i+1]-troot_cpsi[i-1])/(troot_theta[i+1]-troot_theta[i-1])

    # Derivative Check
    aroot_dpsidth  = np.full(shape=np.shape(aroot_theta),dtype=np.float64,fill_value=np.nan)
    aroot_dpsidthc = np.full(shape=np.shape(aroot_theta),dtype=np.float64,fill_value=np.nan)
    for i,th in enumerate(aroot_theta[1:-2]):
        aroot_dpsidth[i] = dpsidth_from_th(ci(pft1), ci(pm_aroot), c8(th))
        aroot_dpsidthc[i] = (aroot_cpsi[i+1]-aroot_cpsi[i-1])/(aroot_theta[i+1]-aroot_theta[i-1])

    fig2, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(12,10))
    ax1.plot(leaf_theta,leaf_dpsidth,label='dpsidth')
    ax1.plot(leaf_theta,leaf_dpsidthc,label='DPsi/Dtheta')
    ax1.set_ylim((-10,100))
    ax1.legend(loc='lower left')
    ax1.set_title('Leaf')
    ax1.set_ylabel('dpsi/dth')

    ax2.plot(stem_theta,stem_dpsidth)
    ax2.plot(stem_theta,stem_dpsidthc)
    ax2.set_ylim((-10,100))
    ax2.set_title('Stem')
    ax2.legend(loc='lower left')

    ax3.plot(troot_theta,troot_dpsidth,label='dpsidth')
    ax3.plot(troot_theta,troot_dpsidthc,label='DPsi/Dtheta')
    ax3.set_ylim((-10,100))
    ax3.set_title('TRoot')
    ax3.set_ylabel('dpsi/dth')
    ax3.set_xlabel('theta')

    ax4.plot(aroot_theta,aroot_dpsidth,label='dpsidth')
    ax4.plot(aroot_theta,aroot_dpsidthc,label='DPsi/Dtheta')
    ax4.set_ylim((-10,100))
    ax4.set_title('ARoot')
    ax4.set_xlabel('theta')


    # Plot out FTC/PSI
    # Find PSI for each theta
    leaf_cflc      = np.full(shape=np.shape(aroot_theta),dtype=np.float64,fill_value=np.nan)
    leaf_dflcdpsi  = np.full(shape=np.shape(aroot_theta),dtype=np.float64,fill_value=np.nan)
    leaf_dflcdpsic = np.full(shape=np.shape(aroot_theta),dtype=np.float64,fill_value=np.nan)
    stem_cflc      = np.full(shape=np.shape(aroot_theta),dtype=np.float64,fill_value=np.nan)
    stem_dflcdpsi  = np.full(shape=np.shape(aroot_theta),dtype=np.float64,fill_value=np.nan)
    stem_dflcdpsic = np.full(shape=np.shape(aroot_theta),dtype=np.float64,fill_value=np.nan)
    troot_cflc      = np.full(shape=np.shape(aroot_theta),dtype=np.float64,fill_value=np.nan)
    troot_dflcdpsi  = np.full(shape=np.shape(aroot_theta),dtype=np.float64,fill_value=np.nan)
    troot_dflcdpsic = np.full(shape=np.shape(aroot_theta),dtype=np.float64,fill_value=np.nan)
    aroot_cflc      = np.full(shape=np.shape(aroot_theta),dtype=np.float64,fill_value=np.nan)
    aroot_dflcdpsi  = np.full(shape=np.shape(aroot_theta),dtype=np.float64,fill_value=np.nan)
    aroot_dflcdpsic = np.full(shape=np.shape(aroot_theta),dtype=np.float64,fill_value=np.nan)


    for i,psi in enumerate(leaf_cpsi):
        leaf_cflc[i]      = flc_from_psi(ci(pft1), ci(pm_leaf), c8(psi))
        leaf_dflcdpsi[i]  = dflcdpsi_from_psi(ci(pft1),ci(pm_leaf),c8(psi))

    for i,psi in enumerate(stem_cpsi):
        stem_cflc[i]     = flc_from_psi(ci(pft1), ci(pm_stem), c8(psi))
        stem_dflcdpsi[i] = dflcdpsi_from_psi(ci(pft1),ci(pm_stem),c8(psi))

    for i,psi in enumerate(troot_cpsi):
        troot_cflc[i] = flc_from_psi(ci(pft1), ci(pm_troot), c8(psi))
        troot_dflcdpsi[i] = dflcdpsi_from_psi(ci(pft1),ci(pm_troot),c8(psi))

    for i,psi in enumerate(aroot_cpsi):
        aroot_cflc[i] = flc_from_psi(ci(pft1), ci(pm_aroot), c8(psi))
        aroot_dflcdpsi[i] = dflcdpsi_from_psi(ci(pft1),ci(pm_aroot),c8(psi))

    # back-calculate the derivative
    for i,psi in enumerate(leaf_cpsi[1:-2]):
        leaf_dflcdpsic[i] = (leaf_cflc[i+1]-leaf_cflc[i-1]) / \
                            (leaf_cpsi[i+1]-leaf_cpsi[i-1])

    for i,psi in enumerate(stem_cpsi[1:-2]):
        stem_dflcdpsic[i] = (stem_cflc[i+1]-stem_cflc[i-1]) / \
                            (stem_cpsi[i+1]-stem_cpsi[i-1])

    for i,psi in enumerate(troot_cpsi[1:-2]):
        troot_dflcdpsic[i] = (troot_cflc[i+1]-troot_cflc[i-1]) / \
                             (troot_cpsi[i+1]-troot_cpsi[i-1])

    for i,psi in enumerate(aroot_cpsi[1:-2]):
        aroot_dflcdpsic[i] = (aroot_cflc[i+1]-aroot_cflc[i-1]) / \
                             (aroot_cpsi[i+1]-aroot_cpsi[i-1])

    fig3, ax1 = plt.subplots(1)
    ax1.plot(leaf_cpsi,leaf_cflc,label='Leaf')
    ax1.plot(stem_cpsi,stem_cflc,label='Stem')
    ax1.plot(troot_cpsi,troot_cflc,label='Troot')
    ax1.plot(aroot_cpsi,aroot_cflc,label='Aroot')
    ax1.grid(True)
    ax1.set_ylabel('FTC')
    ax1.set_xlabel('Psi')
    ax2.legend(loc='upper left')
    plt.tight_layout()


    fig4, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(12,10))
    ax1.plot(leaf_cpsi,leaf_dflcdpsi,label='function')
    ax1.plot(leaf_cpsi,leaf_dflcdpsic,label='discrete')
    ax1.legend(loc='upper left')
    ax1.set_ylabel('dFLC/dPSI')
    ax1.set_title('Leaf')

    ax2.plot(stem_cpsi,stem_dflcdpsi)
    ax2.plot(stem_cpsi,stem_dflcdpsic)
    ax2.set_title('Stem')

    ax3.plot(leaf_cpsi,leaf_dflcdpsi)
    ax3.plot(leaf_cpsi,leaf_dflcdpsic)
    ax3.set_title('TRoot')
    ax3.set_xlabel('Psi')
    ax3.set_ylabel('dFLC/dPSI')

    ax4.plot(leaf_cpsi,leaf_dflcdpsi,label='dpsidth')
    ax4.plot(leaf_cpsi,leaf_dflcdpsic,label='DPsi/Dtheta')
    ax4.set_title('ARoot')
    ax4.set_xlabel('Psi')


    plt.show()

#    code.interact(local=dict(globals(), **locals()))

# =======================================================================================
# This is the actual call to main

if __name__ == "__main__":
    main(sys.argv)
