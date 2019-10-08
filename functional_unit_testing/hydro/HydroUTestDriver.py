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
from PyF90Utils import c8, ci, cchar, c8_arr, ci_arr

# Load the fortran objects via CTYPES

f90_edparams_obj = ctypes.CDLL('bld/EDParamsHydroMod.o',mode=ctypes.RTLD_GLOBAL)
f90_constants_obj = ctypes.CDLL('bld/FatesConstantsMod.o',mode=ctypes.RTLD_GLOBAL)
f90_unitwrap_obj = ctypes.CDLL('bld/UnitWrapMod.o',mode=ctypes.RTLD_GLOBAL)
f90_wftfuncs_obj = ctypes.CDLL('bld/FatesHydroWTFMod.o',mode=ctypes.RTLD_GLOBAL)
f90_hydrounitwrap_obj = ctypes.CDLL('bld/HydroUnitWrapMod.o',mode=ctypes.RTLD_GLOBAL)

# Alias the F90 functions, specify the return type
# -----------------------------------------------------------------------------------

initalloc_wtfs = f90_hydrounitwrap_obj.__hydrounitwrapmod_MOD_initallocwtfs
setwrf = f90_hydrounitwrap_obj.__hydrounitwrapmod_MOD_setwrf
setwkf = f90_hydrounitwrap_obj.__hydrounitwrapmod_MOD_setwkf
th_from_psi = f90_hydrounitwrap_obj.__hydrounitwrapmod_MOD_wrapthfrompsi
th_from_psi.restype = c_double
psi_from_th = f90_hydrounitwrap_obj.__hydrounitwrapmod_MOD_wrappsifromth
psi_from_th.restype = c_double
dpsidth_from_th = f90_hydrounitwrap_obj.__hydrounitwrapmod_MOD_wrapdpsidth
dpsidth_from_th.restype = c_double
ftc_from_psi = f90_hydrounitwrap_obj.__hydrounitwrapmod_MOD_wrapftcfrompsi
ftc_from_psi.restype = c_double
dftcdpsi_from_psi = f90_hydrounitwrap_obj.__hydrounitwrapmod_MOD_wrapftcfrompsi
dftcdpsi_from_psi.restype = c_double


# Some constants
rwcft  = [1.0,0.958,0.958,0.958]
rwccap = [1.0,0.947,0.947,0.947]
pm_leaf = 1
pm_stem = 2
pm_troot = 3
pm_aroot = 4
pm_rhiz = 5

isoil1 = 0  # Top soil layer parameters (@BCI)
isoil2 = 1  # Bottom soil layer parameters

# Constants for rhizosphere
watsat = [0.567, 0.444]
sucsat = [159.659, 256.094]
bsw    = [6.408, 9.27]

unconstrained = True


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

    # Set number of analysis points
    npts = 1000


    #    min_theta = np.full(shape=(2),dtype=np.float64,fill_value=np.nan)

    th_ress = [0.10, 0.20, 0.30, 0.40]
    th_sats = [0.65, 0.65, 0.65, 0.65]
    alphas  = [1.0, 1.0, 1.0, 1.0]
    psds    = [2.7, 2.7, 2.7, 2.7]
    tort    = [0.5, 0.5, 0.5, 0.5]

    nwrf = 4
    nwkf = 2


    # Allocate memory to our objective classes
    iret = initalloc_wtfs(ci(nwrf),ci(nwkf))
    print('Allocated')

    # Push parameters to those classes
    # -------------------------------------------------------------------------
    # Generic VGs
    init_wrf_args  = (4 * c_double)(alphas[0],psds[0],th_sats[0],th_ress[0]) # alpha, psd, th_sat, th_res
    iret = setwrf(ci(1),ci(1),ci(len(init_wrf_args)),byref(init_wrf_args))

    init_wrf_args  = (4 * c_double)(alphas[1],psds[1],th_sats[1],th_ress[1]) # alpha, psd, th_sat, th_res
    iret = setwrf(ci(2),ci(1),ci(len(init_wrf_args)),byref(init_wrf_args))

    init_wrf_args  = (4 * c_double)(alphas[2],psds[2],th_sats[2],th_ress[2]) # alpha, psd, th_sat, th_res
    iret = setwrf(ci(3),ci(1),ci(len(init_wrf_args)),byref(init_wrf_args))

    init_wrf_args  = (4 * c_double)(alphas[3],psds[3],th_sats[3],th_ress[3]) # alpha, psd, th_sat, th_res
    iret = setwrf(ci(4),ci(1),ci(len(init_wrf_args)),byref(init_wrf_args))




    print('initialized WRF')

    theta = np.linspace(0.10, 0.7, num=npts)

    psi   = np.full(shape=(nwrf,len(theta)),dtype=np.float64,fill_value=np.nan)
    dpsidth = np.full(shape=(nwrf,len(theta)),dtype=np.float64,fill_value=np.nan)
    cdpsidth = np.full(shape=(nwrf,len(theta)),dtype=np.float64,fill_value=np.nan)

    for ic in range(nwrf):

        for i,th in enumerate(theta):
            psi[ic,i] = psi_from_th(ci(ic+1),c8(th))


    fig0, ax1 = plt.subplots(1,1,figsize=(9,6))
    for ic in range(nwrf):
        ax1.plot(theta,psi[ic,:])

    ax1.set_ylim((-10,5))
    ax1.set_ylabel('Matric Potential [MPa]')
    ax1.set_xlabel('VWC [m3/m3]')


    for ic in range(nwrf):
        for i in range(1,len(theta)-1):
            dpsidth[ic,i]  = dpsidth_from_th(ci(ic+1),c8(theta[i]))
            cdpsidth[ic,i] = (psi[ic,i+1]-psi[ic,i-1])/(theta[i+1]-theta[i-1])


    fig1, ax1 = plt.subplots(1,1,figsize=(9,6))
    for ic in range(nwrf):
        ax1.plot(theta,dpsidth[0,:],label='func')
        ax1.plot(theta,cdpsidth[0,:],label='check')
    ax1.set_ylim((0,1000))

    ax1.set_ylabel('dPSI/dTh [MPa m3 m-3]')
    ax1.set_xlabel('VWC [m3/m3]')
    ax1.legend(loc='upper right')
    plt.show()



    # Push parameters to WKF classes
    # -------------------------------------------------------------------------
    # Generic VGs
    init_wkf_args  = c8_arr(alphas[0],psds[0],th_sats[0],th_ress[0],tort[0]) # alpha, psd, th_sat, th_res
    iret = setwkf(ci(1),ci(1),ci(len(init_wkf_args)),byref(init_wkf_args))



    code.interact(local=dict(globals(), **locals()))

    exit(0)


    # Test 1 For a set of thetas, calculate psi for each pm.
    # ===================================================================================

    pft1 = 1
    pft2 = 2



    if(unconstrained):
        min_leaf_theta = 0.01
        max_leaf_theta = 0.99
        min_stem_theta = 0.01
        max_stem_theta = 0.99
        min_troot_theta = 0.01
        max_troot_theta = 0.99
        min_aroot_theta = 0.01
        max_aroot_theta = 0.99
        min_leaf_theta2 = 0.01
        min_stem_theta2 = 0.01
        min_troot_theta2 = 0.01
        min_aroot_theta2 = 0.01
        min_rhiz_theta = 0.01
        max_rhiz_theta = 0.99
    else:
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
        min_troot_theta2 = pftparms['hydr_resid_node'].data[pm_troot-1,pft2-1]
        min_aroot_theta2 = pftparms['hydr_resid_node'].data[pm_aroot-1,pft2-1]
        min_rhiz_theta   = 0.01
        max_rhiz_theta   = watsat[isoil1]

    # Rhizosphere
    # -----------------------------------------------------------------------------------

    rhiz_theta  = np.linspace(min_rhiz_theta,max_rhiz_theta, num=npts)
    rhiz_psi = np.full(shape=np.shape(rhiz_theta),dtype=np.float64,fill_value=np.nan)
    rhiz_psi2 = np.full(shape=np.shape(rhiz_theta),dtype=np.float64,fill_value=np.nan)
    rhiz_dpsidth  = np.full(shape=np.shape(rhiz_theta),dtype=np.float64,fill_value=np.nan)
    rhiz_dpsidthc = np.full(shape=np.shape(rhiz_theta),dtype=np.float64,fill_value=np.nan)
    rhiz_flc = np.full(shape=np.shape(rhiz_theta),dtype=np.float64,fill_value=np.nan)
    rhiz_dflcdpsi = np.full(shape=np.shape(rhiz_theta),dtype=np.float64,fill_value=np.nan)
    rhiz_dflcdpsic = np.full(shape=np.shape(rhiz_theta),dtype=np.float64,fill_value=np.nan)



    # Initialize Theta
    leaf_theta  = np.linspace(min_leaf_theta,max_leaf_theta, num=npts)
    stem_theta  = np.linspace(min_stem_theta,max_stem_theta, num=npts)
    troot_theta = np.linspace(min_troot_theta,max_troot_theta, num=npts)
    aroot_theta = np.linspace(min_aroot_theta,max_aroot_theta, num=npts)
    leaf_theta2  = np.linspace(min_leaf_theta2,max_leaf_theta, num=npts)

    # Initialize PSI
    leaf_psi  = np.full(shape=np.shape(leaf_theta),dtype=np.float64,fill_value=np.nan)
    leaf_psi2 = np.full(shape=np.shape(leaf_theta2),dtype=np.float64,fill_value=np.nan)
    stem_psi  = np.full(shape=np.shape(stem_theta),dtype=np.float64,fill_value=np.nan)
    troot_psi = np.full(shape=np.shape(troot_theta),dtype=np.float64,fill_value=np.nan)
    aroot_psi = np.full(shape=np.shape(aroot_theta),dtype=np.float64,fill_value=np.nan)

    # Initialize dPSI/dtheta derivative and discrete check
    leaf_dpsidth  = np.full(shape=np.shape(leaf_theta),dtype=np.float64,fill_value=np.nan)
    leaf_dpsidthc = np.full(shape=np.shape(leaf_theta),dtype=np.float64,fill_value=np.nan)
    stem_dpsidth  = np.full(shape=np.shape(stem_theta),dtype=np.float64,fill_value=np.nan)
    stem_dpsidthc = np.full(shape=np.shape(stem_theta),dtype=np.float64,fill_value=np.nan)
    troot_dpsidth  = np.full(shape=np.shape(troot_theta),dtype=np.float64,fill_value=np.nan)
    troot_dpsidthc = np.full(shape=np.shape(troot_theta),dtype=np.float64,fill_value=np.nan)
    aroot_dpsidth  = np.full(shape=np.shape(aroot_theta),dtype=np.float64,fill_value=np.nan)
    aroot_dpsidthc = np.full(shape=np.shape(aroot_theta),dtype=np.float64,fill_value=np.nan)


    # Initialize the FLC and its derivative
    leaf_flc      = np.full(shape=np.shape(leaf_theta),dtype=np.float64,fill_value=np.nan)
    leaf_dflcdpsi  = np.full(shape=np.shape(leaf_theta),dtype=np.float64,fill_value=np.nan)
    leaf_dflcdpsic = np.full(shape=np.shape(leaf_theta),dtype=np.float64,fill_value=np.nan)

    stem_flc      = np.full(shape=np.shape(stem_theta),dtype=np.float64,fill_value=np.nan)
    stem_dflcdpsi  = np.full(shape=np.shape(stem_theta),dtype=np.float64,fill_value=np.nan)
    stem_dflcdpsic = np.full(shape=np.shape(stem_theta),dtype=np.float64,fill_value=np.nan)

    troot_flc      = np.full(shape=np.shape(troot_theta),dtype=np.float64,fill_value=np.nan)
    troot_dflcdpsi  = np.full(shape=np.shape(troot_theta),dtype=np.float64,fill_value=np.nan)
    troot_dflcdpsic = np.full(shape=np.shape(troot_theta),dtype=np.float64,fill_value=np.nan)

    aroot_flc      = np.full(shape=np.shape(aroot_theta),dtype=np.float64,fill_value=np.nan)
    aroot_dflcdpsi  = np.full(shape=np.shape(aroot_theta),dtype=np.float64,fill_value=np.nan)
    aroot_dflcdpsic = np.full(shape=np.shape(aroot_theta),dtype=np.float64,fill_value=np.nan)



    mpl.rcParams.update({'font.size': 15})


    # Initialize the return variable
    cpsi = c_double(0)


    # Find PSI for each theta
    for i,th in enumerate(leaf_theta):
        leaf_psi[i] = psi_from_th(ci(pft1), ci(pm_leaf), c8(th))

    for i,th in enumerate(leaf_theta2):
        leaf_psi2[i] = psi_from_th(ci(pft2),ci(pm_leaf),c8(th))

    for i,th in enumerate(stem_theta):
        stem_psi[i] = psi_from_th(ci(pft1), ci(pm_stem), c8(th))

    for i,th in enumerate(troot_theta):
        troot_psi[i] = psi_from_th(ci(pft1), ci(pm_troot), c8(th))

    for i,th in enumerate(aroot_theta):
        aroot_psi[i] = psi_from_th(ci(pft1), ci(pm_aroot), c8(th))

    for i,th in enumerate(rhiz_theta):
        rhiz_psi[i] = psi_from_th(ci(pft1), ci(pm_rhiz), c8(th), \
                                  c8(watsat[0]), c8(sucsat[0]), c8(bsw[0]))
        rhiz_psi2[i] = psi_from_th(ci(pft1), ci(pm_rhiz), c8(th), \
                                   c8(watsat[1]), c8(sucsat[1]), c8(bsw[1]))


    fig0, (ax1,ax2) = plt.subplots(1,2,figsize=(9,6))
    ax1.plot(leaf_theta,leaf_psi,label='Leaf')
    ax1.plot(stem_theta,stem_psi,label='Stem')
    ax1.plot(troot_theta,troot_psi,label='Troot')
    ax1.plot(aroot_theta,aroot_psi,label='Aroot')
    ax1.plot(rhiz_theta,rhiz_psi,label='Rhiz')
    ax1.grid(True)
    ax1.set_ylabel('Psi')
    ax1.set_xlim((0,1))
    ax1.set_ylim((-20,10))
    ax1.set_xlabel('Theta')
    ax1.set_title('PFT: {}'.format(pft1))
    ax1.legend(loc='lower right')
    ax2.plot(leaf_theta,semilogneg(leaf_psi),label='Leaf')
    ax2.plot(stem_theta,semilogneg(stem_psi),label='Stem')
    ax2.plot(troot_theta,semilogneg(troot_psi),label='Troot')
    ax2.plot(aroot_theta,semilogneg(aroot_psi),label='Aroot')
    ax2.plot(rhiz_theta,semilogneg(rhiz_psi),label='Rhiz')
    ax2.grid(True)
    ax2.set_ylabel('log(Psi)')
    ax2.set_xlim((0,1))
    ax2.set_xlabel('Theta')


    plt.tight_layout()


    # Derivative Check on PSI
    # -----------------------------------------------------------------------------------
    for i in range(1,len(leaf_theta)-1):
        leaf_dpsidth[i] = dpsidth_from_th(ci(pft1), ci(pm_leaf), c8(leaf_theta[i]))
        leaf_dpsidthc[i] = (leaf_psi[i+1]-leaf_psi[i-1])/(leaf_theta[i+1]-leaf_theta[i-1])

    for i in range(1,len(stem_theta)-1):
        stem_dpsidth[i] = dpsidth_from_th(ci(pft1), ci(pm_stem), c8(stem_theta[i]))
        stem_dpsidthc[i] = (stem_psi[i+1]-stem_psi[i-1])/(stem_theta[i+1]-stem_theta[i-1])

    for i in range(1,len(troot_theta)-1):
        troot_dpsidth[i] = dpsidth_from_th(ci(pft1), ci(pm_troot), c8(troot_theta[i]))
        troot_dpsidthc[i] = (troot_psi[i+1]-troot_psi[i-1])/(troot_theta[i+1]-troot_theta[i-1])

    for i in range(1,len(aroot_theta)-1):
        aroot_dpsidth[i] = dpsidth_from_th(ci(pft1), ci(pm_aroot), c8(aroot_theta[i]))
        aroot_dpsidthc[i] = (aroot_psi[i+1]-aroot_psi[i-1])/(aroot_theta[i+1]-aroot_theta[i-1])

    for i in range(1,len(rhiz_theta)-1):
        rhiz_dpsidth[i]  = dpsidth_from_th(ci(pft1), ci(pm_rhiz), \
                                           c8(rhiz_theta[i]), c8(watsat[0]), \
                                           c8(sucsat[0]), c8(bsw[0]))
        rhiz_dpsidthc[i] = (rhiz_psi[i+1]-rhiz_psi[i-1])/ \
                           (rhiz_theta[i+1]-rhiz_theta[i-1])


    fig2, ((ax1,ax2),(ax3,ax4),(ax5,ax6)) = plt.subplots(3,2,figsize=(9,11))
    ax1.plot(leaf_theta,leaf_dpsidth,label='function')
    ax1.plot(leaf_theta,leaf_dpsidthc,label='discrete')
    ax1.set_xlim((0,1))
    ax1.legend(loc='upper right')
    ax1.set_title('Leaf')
    ax1.set_ylabel('dpsi/dth')
    ax1.grid(True)
    ax2.plot(stem_theta,stem_dpsidth)
    ax2.plot(stem_theta,stem_dpsidthc)
    ax2.set_xlim((0,1))
    ax2.set_title('Stem')
    ax2.grid(True)
    ax3.plot(troot_theta,troot_dpsidth)
    ax3.plot(troot_theta,troot_dpsidthc)
    ax3.set_xlim((0,1))
    ax3.set_title('TRoot')
    ax3.set_ylabel('dpsi/dth')
    ax3.grid(True)
    ax4.plot(aroot_theta,aroot_dpsidth)
    ax4.plot(aroot_theta,aroot_dpsidthc)
    ax4.set_xlim((0,1))
    ax4.set_title('ARoot')
    ax4.set_xlabel('theta')
    ax4.grid(True)
    ax5.plot(rhiz_theta,rhiz_dpsidth)
    ax5.plot(rhiz_theta,rhiz_dpsidthc)
    ax5.set_xlim((0,1))
    ax5.set_title('Rhiz')
    ax5.set_xlabel('theta')
    ax5.set_ylabel('dpsi/dth')
    ax5.grid(True)
    ax6.axis('off')

    plt.tight_layout()

    # Plot out FTC/PSI
    # Find PSI for each theta



    for i,psi in enumerate(leaf_psi):
        leaf_flc[i]      = flc_from_psi(ci(pft1), ci(pm_leaf), c8(leaf_theta[i]), c8(psi))

    for i,psi in enumerate(stem_psi):
        stem_flc[i]     = flc_from_psi(ci(pft1), ci(pm_stem), c8(stem_theta[i]), c8(psi))

    for i,psi in enumerate(troot_psi):
        troot_flc[i] = flc_from_psi(ci(pft1), ci(pm_troot), c8(troot_theta[i]), c8(psi))


    for i,psi in enumerate(aroot_psi):
        aroot_flc[i] = flc_from_psi(ci(pft1), ci(pm_aroot), c8(aroot_theta[i]), c8(psi))

    for i,psi in enumerate(rhiz_psi):
        rhiz_flc[i]  = flc_from_psi(ci(pft1), ci(pm_rhiz), c8(rhiz_theta[i]), \
                                    c8(psi), c8(sucsat[isoil1]), c8(bsw[isoil1]))


    # back-calculate the derivative
    for i in range(1,len(leaf_psi)-1):
        leaf_dflcdpsi[i]  = dflcdpsi_from_psi(ci(pft1),ci(pm_leaf), c8(leaf_theta[i]), c8(leaf_psi[i]))
        leaf_dflcdpsic[i] = (leaf_flc[i+1]-leaf_flc[i-1]) / \
                            (leaf_psi[i+1]-leaf_psi[i-1])

    for i in range(1,len(stem_psi)-1):
        stem_dflcdpsi[i] = dflcdpsi_from_psi(ci(pft1),ci(pm_stem), c8(stem_theta[i]), c8(stem_psi[i]))
        stem_dflcdpsic[i] = (stem_flc[i+1]-stem_flc[i-1]) / \
                            (stem_psi[i+1]-stem_psi[i-1])

    for i in range(1,len(troot_psi)-1):
        troot_dflcdpsi[i] = dflcdpsi_from_psi(ci(pft1),ci(pm_troot), c8(troot_theta[i]), c8(troot_psi[i]))
        troot_dflcdpsic[i] = (troot_flc[i+1]-troot_flc[i-1]) / \
                             (troot_psi[i+1]-troot_psi[i-1])

    for i in range(1,len(aroot_psi)-1):
        aroot_dflcdpsi[i] = dflcdpsi_from_psi(ci(pft1),ci(pm_aroot), c8(aroot_theta[i]), c8(aroot_psi[i]))
        aroot_dflcdpsic[i] = (aroot_flc[i+1]-aroot_flc[i-1]) / \
                             (aroot_psi[i+1]-aroot_psi[i-1])

    for i in range(1,len(rhiz_psi)-1):
        rhiz_dflcdpsi[i] = dflcdpsi_from_psi(ci(pft1), ci(pm_rhiz), c8(rhiz_theta[i]), \
                                             c8(rhiz_psi[i]), c8(sucsat[0]), c8(bsw[0]))
        rhiz_dflcdpsic[i] = (rhiz_flc[i+1] - rhiz_flc[i-1])/(rhiz_psi[i+1]-rhiz_psi[i-1])


    fig3, (ax1,ax2) = plt.subplots(1,2,figsize=(9,6))
    ax1.plot(leaf_psi,leaf_flc,label='Leaf')
    ax1.plot(stem_psi,stem_flc,label='Stem')
    ax1.plot(troot_psi,troot_flc,label='Troot')
    ax1.plot(aroot_psi,aroot_flc,label='Aroot')
    ax1.plot(rhiz_psi,rhiz_flc,label='Rhiz')
    ax1.grid(True)
    ax1.set_ylabel('FTC [-]')
    ax1.set_xlabel('Psi [MPa]')
    ax1.legend(loc='upper left')
    ax1.set_title('PFT: {}'.format(pft1))
    ax2.plot(leaf_theta,leaf_flc,label='leaf')
    ax2.plot(stem_theta,stem_flc,label='stem')
    ax2.plot(troot_theta,troot_flc,label='troot')
    ax2.plot(aroot_theta,aroot_flc,label='aroot')
    ax2.plot(rhiz_theta,rhiz_flc,label='rhiz')
    ax2.grid(True)
    ax2.set_ylabel('FTC [-]')
    ax2.set_xlabel('Theta [m3/m3]')
    plt.tight_layout()

    fig4, ((ax1,ax2),(ax3,ax4),(ax5,ax6)) = plt.subplots(3,2,figsize=(9,11))
    ax1.plot(leaf_psi,leaf_dflcdpsi,label='function')
    ax1.plot(leaf_psi,leaf_dflcdpsic,label='discrete')
    ax1.legend(loc='upper left')
    ax1.set_ylabel('dFLC/dPsi')
    ax1.set_title('Leaf')
    ax1.grid(True)
    ax2.plot(stem_psi,stem_dflcdpsi)
    ax2.plot(stem_psi,stem_dflcdpsic)
    ax2.set_title('Stem')
    ax2.grid(True)
    ax3.plot(leaf_psi,leaf_dflcdpsi)
    ax3.plot(leaf_psi,leaf_dflcdpsic)
    ax3.set_title('TRoot')
    ax3.set_ylabel('dFLC/dPsi')
    ax3.grid(True)
    ax4.plot(leaf_psi,leaf_dflcdpsi)
    ax4.plot(leaf_psi,leaf_dflcdpsic)
    ax4.set_title('ARoot')
    ax4.set_xlabel('Psi')
    ax4.grid(True)
    ax5.plot(rhiz_psi,semilogneg(rhiz_dflcdpsi))
    ax5.plot(rhiz_psi,semilogneg(rhiz_dflcdpsic))
    ax5.set_title('Rhiz')
    ax5.set_xlabel('Psi')
    ax5.set_ylabel('log(dFLC/dPsi)')
    ax5.grid(True)
    ax6.axis('off')
    plt.tight_layout()


    fig44, ax1 = plt.subplots(1,figsize=(7,7))
    ax1.plot(leaf_theta,leaf_dflcdpsi*leaf_dpsidth,label='leaf')
    ax1.plot(stem_theta,stem_dflcdpsi*stem_dpsidth,label='stem')
    ax1.plot(troot_theta,troot_dflcdpsi*troot_dpsidth,label='troot')
    ax1.plot(aroot_theta,aroot_dflcdpsi*aroot_dpsidth,label='aroot')
    ax1.plot(rhiz_theta,rhiz_dflcdpsi*rhiz_dpsidth,label='rhiz')
    ax1.legend(loc='upper left')
    ax1.set_ylabel('dFLC/dtheta')
    ax1.set_xlabel('theta')
    ax1.grid(True)
    plt.tight_layout()



    fig5, (ax1,ax2) = plt.subplots(2)
    ax1.plot(rhiz_theta,rhiz_psi,label='Sat={}, PSIsat={}, B={}'.format(watsat[0],-sucsat[0]*9.8*1.e-9*1000.0 ,bsw[0]))
    ax1.plot(rhiz_theta,rhiz_psi2,label='Sat={}, PSIsat={}, B={}'.format(watsat[1],-sucsat[1]*9.8*1.e-9*1000.0 ,bsw[1]))
    ax1.grid(True)
    ax1.set_ylabel('Psi [MPa]')
    ax1.set_xlim((0,1))
    ax1.set_xlabel('Theta [m3/m3]')
    ax1.set_title('Rhizosphere')
    ax1.legend(loc='lower right')
    plt.tight_layout()

    ax2.plot(leaf_theta,leaf_psi,label='Leaf')
    ax2.plot(stem_theta,stem_psi,label='Stem')
    ax2.plot(troot_theta,troot_psi,label='Troot')
    ax2.plot(aroot_theta,aroot_psi,label='Aroot')
    ax2.grid(True)
    ax2.set_ylabel('Psi [MPa]')
    ax2.set_xlim((0,1))
    ax2.set_xlabel('Theta [m3/m3]')
    ax2.set_title('PFT: {}'.format(pft1))
    ax2.legend(loc='lower right')
    plt.tight_layout()





    plt.show()

#    code.interact(local=dict(globals(), **locals()))

# Helper code to plot negative logs

def semilogneg(x):

    y = np.sign(x)*np.log(abs(x))
    return(y)

def semilog10net(x):

    y = np.sign(x)*np.log10(abs(x))
    return(y)


# =======================================================================================
# This is the actual call to main

if __name__ == "__main__":
    main(sys.argv)
