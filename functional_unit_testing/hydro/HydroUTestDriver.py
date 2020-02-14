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

f90_unitwrap_obj = ctypes.CDLL('bld/UnitWrapMod.o',mode=ctypes.RTLD_GLOBAL)
f90_constants_obj = ctypes.CDLL('bld/FatesConstantsMod.o',mode=ctypes.RTLD_GLOBAL)
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
dftcdpsi_from_psi = f90_hydrounitwrap_obj.__hydrounitwrapmod_MOD_wrapdftcdpsi
dftcdpsi_from_psi.restype = c_double


# Some constants
rwcft  = [1.0,0.958,0.958,0.958]
rwccap = [1.0,0.947,0.947,0.947]
pm_leaf = 1
pm_stem = 2
pm_troot = 3
pm_aroot = 4
pm_rhiz = 5

# These parameters are matched with the indices in FATES-HYDRO
vg_type = 1
cch_type = 2
tfs_type = 3

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


class vg_wrf:
    def __init__(self,index,alpha, psd, th_sat, th_res):
        self.alpha = alpha
        self.psd   = psd
        self.th_sat = th_sat
        self.th_res = th_res
        init_wrf_args = [self.alpha, self.psd, self.th_sat, self.th_res]
        iret = setwrf(ci(index),ci(vg_type),ci(len(init_wrf_args)),c8_arr(init_wrf_args))

class cch_wrf:
    def __init__(self,index,th_sat,psi_sat,beta):
        self.th_sat  = th_sat
        self.psi_sat = psi_sat
        self.beta    = beta
        init_wrf_args = [self.th_sat,self.psi_sat,self.beta]
        iret = setwrf(ci(index),ci(cch_type),ci(len(init_wrf_args)),c8_arr(init_wrf_args))

class vg_wkf:
    def __init__(self,index,alpha, psd, th_sat, th_res, tort):
        self.alpha  = alpha
        self.psd    = psd
        self.th_sat = th_sat
        self.th_res = th_res
        self.tort   = tort
        init_wkf_args = [self.alpha, self.psd,self.th_sat,self.th_res,self.tort]
        iret = setwkf(ci(index),ci(vg_type),ci(len(init_wkf_args)),c8_arr(init_wkf_args))

class cch_wkf:
    def __init__(self,index,th_sat,psi_sat,beta):
        self.th_sat  = th_sat
        self.psi_sat = psi_sat
        self.beta    = beta
        init_wkf_args = [self.th_sat,self.psi_sat,self.beta]
        iret = setwkf(ci(index),ci(cch_type),ci(len(init_wkf_args)),c8_arr(init_wkf_args))


class tfs_wrf:
    def __init__(self,index,th_sat,th_res,pinot,epsil,rwc_fd,cap_corr,cap_int,cap_slp,pmedia):
        self.th_sat = th_sat
        self.th_res = th_res
        self.pinot  = pinot
        self.epsil  = epsil
        self.rwc_fd = rwc_fd
        self.cap_corr = cap_corr
        self.cap_int  = cap_int
        self.cap_slp  = cap_slp
        self.pmedia   = pmedia
        init_wrf_args = [self.th_sat,self.th_res,self.pinot,self.epsil,self.rwc_fd,self.cap_corr,self.cap_int,self.cap_slp,self.pmedia]
        iret = setwrf(ci(index),ci(tfs_type),ci(len(init_wrf_args)),c8_arr(init_wrf_args))

class tfs_wkf:
    def __init__(self,index,p50,avuln):
        self.avuln = avuln
        self.p50   = p50
        init_wkf_args = [self.p50,self.avuln]
        iret = setwkf(ci(index),ci(tfs_type),ci(len(init_wkf_args)),c8_arr(init_wkf_args))


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

#    parser = argparse.ArgumentParser(description='Parse command line arguments to this script.')
#    parser.add_argument('--cdl-file', dest='cdlfile', type=str, \
#                        help="Input CDL filename.  Required.", required=True)

#    args = parser.parse_args()


    # Set number of analysis points
    npts = 1000


    #    min_theta = np.full(shape=(2),dtype=np.float64,fill_value=np.nan)

#    wrf_type = [vg_type, vg_type, cch_type, cch_type]
#    wkf_type = [vg_type, tfs_type, cch_type, tfs_type]

#    th_ress = [0.01, 0.10, -9, -9]
#    th_sats = [0.55, 0.55, 0.65, 0.65]
#    alphas  = [1.0, 1.0, 1.0, 1.0]
#    psds    = [2.7, 2.7, 2.7, 2.7]
#    tort    = [0.5, 0.5, 0.5, 0.5]
#    beta    = [-9, -9, 6, 9]
#    avuln   = [2.0, 2.0, 2.5, 2.5]
#    p50     = [-1.5, -1.5, -2.25, -2.25]

    ncomp= 3

    rwc_fd  = [1.0,0.958,0.958,0.958]
    rwccap  = [1.0,0.947,0.947,0.947]
    cap_slp = []
    cap_int = []
    cap_corr= []
    hydr_psi0 = 0.0
    hydr_psicap = -0.6

    for pm in range(4):
        if (pm == 0):
            cap_slp.append(0.0)
            cap_int.append(0.0)
            cap_corr.append(1.0)
        else:
            cap_slp.append((hydr_psi0 - hydr_psicap )/(1.0 - rwccap[pm]))
            cap_int.append(-cap_slp[pm] + hydr_psi0)
            cap_corr.append(-cap_int[pm]/cap_slp[pm])


    # Allocate memory to our objective classes
    iret = initalloc_wtfs(ci(ncomp),ci(ncomp))
    print('Allocated')


    # Define the funcions and their parameters
#    vg_wrf(1,alpha=1.0,psd=2.7,th_sat=0.55,th_res=0.1)
#    vg_wkf(1,alpha=1.0,psd=2.7,th_sat=0.55,th_res=0.1,tort=0.5)

    cch_wrf(1,th_sat=0.55, psi_sat=-1.56e-3, beta=6)
    cch_wkf(1,th_sat=0.55, psi_sat=-1.56e-3, beta=6)

#    cch_wrf(3,th_sat=0.55, psi_sat=-1.56e-3, beta=6)
#    tfs_wkf(3,p50=-2.25, avuln=2.0)

    names=['Soil','ARoot','Leaf']

    # Absorving root
    tfs_wrf(2,th_sat=0.75,th_res=0.15,pinot=-1.043478, \
            epsil=8,rwc_fd=rwc_fd[3],cap_corr=cap_corr[3], \
            cap_int=cap_int[3],cap_slp=cap_slp[3],pmedia=4)
    tfs_wkf(2,p50=-2.25, avuln=2.0)

    # Leaf
    tfs_wrf(3,th_sat=0.65,th_res=0.25,pinot=-1.47, \
            epsil=12,rwc_fd=rwc_fd[0],cap_corr=cap_corr[0], \
            cap_int=cap_int[0],cap_slp=cap_slp[0],pmedia=1)
    tfs_wkf(3,p50=-2.25, avuln=2.0)

    print('initialized WRF')

    theta = np.linspace(0.10, 0.7, num=npts)
    psi   = np.full(shape=(ncomp,len(theta)),dtype=np.float64,fill_value=np.nan)
    dpsidth = np.full(shape=(ncomp,len(theta)),dtype=np.float64,fill_value=np.nan)
    cdpsidth = np.full(shape=(ncomp,len(theta)),dtype=np.float64,fill_value=np.nan)

    for ic in range(ncomp):
        for i,th in enumerate(theta):
            psi[ic,i] = psi_from_th(ci(ic+1),c8(th))


    # Theta vs psi plots

    fig0, ax1 = plt.subplots(1,1,figsize=(9,6))
    for ic in range(ncomp):
        ax1.plot(theta,psi[ic,:],label='{}'.format(names[ic]))

    ax1.set_ylim((-30,5))
    ax1.set_ylabel('Matric Potential [MPa]')
    ax1.set_xlabel('VWC [m3/m3]')
    ax1.legend(loc='lower right')

    for ic in range(ncomp):
        for i in range(1,len(theta)-1):
            dpsidth[ic,i]  = dpsidth_from_th(ci(ic+1),c8(theta[i]))
            cdpsidth[ic,i] = (psi[ic,i+1]-psi[ic,i-1])/(theta[i+1]-theta[i-1])


    # Theta vs dpsi_dth (also checks deriv versus explicit)

    fig1, ax1 = plt.subplots(1,1,figsize=(9,6))
    for ic in range(ncomp):
        ax1.plot(theta,dpsidth[0,:],label='func')
        ax1.plot(theta,cdpsidth[0,:],label='check')
    ax1.set_ylim((0,1000))

    ax1.set_ylabel('dPSI/dTh [MPa m3 m-3]')
    ax1.set_xlabel('VWC [m3/m3]')
    ax1.legend(loc='upper right')

    # Push parameters to WKF classes
    # -------------------------------------------------------------------------
    # Generic VGs

    ftc   = np.full(shape=(ncomp,len(theta)),dtype=np.float64,fill_value=np.nan)
    dftcdpsi = np.full(shape=(ncomp,len(theta)),dtype=np.float64,fill_value=np.nan)
    cdftcdpsi = np.full(shape=(ncomp,len(theta)),dtype=np.float64,fill_value=np.nan)

    for ic in range(ncomp):
        for i in range(0,len(theta)):
            ftc[ic,i] = ftc_from_psi(ci(ic+1),c8(psi[ic,i]))

    for ic in range(ncomp):
        for i in range(1,len(theta)-1):
            dftcdpsi[ic,i]  = dftcdpsi_from_psi(ci(ic+1),c8(psi[ic,i]))
            cdftcdpsi[ic,i] = (ftc[ic,i+1]-ftc[ic,i-1])/(psi[ic,i+1]-psi[ic,i-1])


    # FTC versus Psi

    fig2, ax1 = plt.subplots(1,1,figsize=(9,6))
    for ic in range(ncomp):
        ax1.plot(psi[ic,:],ftc[ic,:],label='{}'.format(names[ic]))

    ax1.set_ylabel('FTC')
    ax1.set_xlabel('Psi [MPa]')
    ax1.set_xlim([-5,0])
    ax1.legend(loc='upper right')


    # FTC versus theta

    fig4, ax1 = plt.subplots(1,1,figsize=(9,6))
    for ic in range(ncomp):
        ax1.plot(theta,ftc[ic,:],label='{}'.format(names[ic]))

    ax1.set_ylabel('FTC')
    ax1.set_xlabel('Theta [m3/m3]')
    ax1.legend(loc='lower right')

    # dFTC/dPSI

    fig3,ax1 = plt.subplots(1,1,figsize=(9,6))
    for ic in range(ncomp):
#        ax1.plot(psi[ic,:],abs(dftcdpsi[ic,:]-cdftcdpsi[ic,:])/abs(cdftcdpsi[ic,:]),label='{}'.format(ic))
        ax1.plot(psi[ic,:],cdftcdpsi[ic,:],label='check')

    ax1.set_ylabel('dFTC/dPSI')
    ax1.set_xlabel('Psi [MPa]')
#    ax1.set_xlim([-30,3])
#    ax1.set_ylim([0,10])
    ax1.legend(loc='upper right')
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
