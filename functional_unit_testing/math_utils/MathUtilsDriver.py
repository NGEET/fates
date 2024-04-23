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


#CDLParse = imp.load_source('CDLParse','../shared/py_src/CDLParse.py')
#F90ParamParse = imp.load_source('F90ParamParse','../shared/py_src/F90ParamParse.py')
PyF90Utils = imp.load_source('PyF90Utils','../shared/py_src/PyF90Utils.py')


#from CDLParse import CDLParseDims, CDLParseParam, cdl_param_type
#from F90ParamParse import f90_param_type, GetSymbolUsage, GetPFTParmFileSymbols, MakeListUnique

from PyF90Utils import c8, ci, cchar, c8_arr, ci_arr

# Load the fortran objects via CTYPES

f90_unitwrap_obj = ctypes.CDLL('bld/UnitWrapMod.o',mode=ctypes.RTLD_GLOBAL)
f90_constants_obj = ctypes.CDLL('bld/FatesConstantsMod.o',mode=ctypes.RTLD_GLOBAL)
f90_fatesutils_obj = ctypes.CDLL('bld/FatesUtilsMod.o',mode=ctypes.RTLD_GLOBAL)

# Alias the F90 functions, specify the return type
# -----------------------------------------------------------------------------------

neighbor_dist = f90_fatesutils_obj.__fatesutilsmod_MOD_getneighbordistance
#quadratic_f = f90_fatesutils_obj.__fatesutilsmod_MOD_quadratic_f
quadratic_roots = f90_fatesutils_obj.__fatesutilsmod_MOD_quadraticroots
quadratic_sroots = f90_fatesutils_obj.__fatesutilsmod_MOD_quadraticrootssridharachary

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

def main(argv):

    # First check to make sure python 2.7 is being used
    version = platform.python_version()
    verlist = version.split('.')

    #if( not ((verlist[0] == '2') & (verlist[1] ==  '7') & (int(verlist[2])>=15) )  ):
    #    print("The PARTEH driver mus be run with python 2.7")
    #    print(" with tertiary version >=15.")
    #    print(" your version is {}".format(version))
    #    print(" exiting...")
    #    sys.exit(2)

    # Read in the arguments
    # =======================================================================================

    #    parser = argparse.ArgumentParser(description='Parse command line arguments to this script.')
    #    parser.add_argument('--cdl-file', dest='cdlfile', type=str, \
    #                        help="Input CDL filename.  Required.", required=True)
    #    args = parser.parse_args()

    # Set number of analysis points

    # y = ax2 + bx + c

    a = [1,1,5,1.5]
    b = [-2,7,10,3.2]
    c = [1,12,3,1.1]
    
    cd_r1 = c_double(-9.0)
    cd_r2 = c_double(-9.0)

    r1 = np.zeros([3,1])
    r2 = np.zeros([3,1])
    
    for ic in range(len(a)):

        #iret = quadratic_f(c8(a[ic]),c8(b[ic]),c8(c[ic]),byref(cd_r1),byref(cd_r2))
        #r1[0] = cd_r1.value
        #r2[0] = cd_r2.value
        
        iret = quadratic_roots(c8(a[ic]),c8(b[ic]),c8(c[ic]),byref(cd_r1),byref(cd_r2))
        r1[1] = cd_r1.value
        r2[1] = cd_r2.value

        iret = quadratic_sroots(c8(a[ic]),c8(b[ic]),c8(c[ic]),byref(cd_r1),byref(cd_r2))
        r1[2] = cd_r2.value
        r2[2] = cd_r1.value
        
        print(a[ic],b[ic],c[ic])
        print(r1)
        print(r2)
        
        #PlotQuadAndRoots(a[ic],b[ic],c[ic],r1,r2)
        
        
def PlotQuadAndRoots(a,b,c,d,r1,r2):
    
    fig, axs = plt.subplots(ncols=1,nrows=1,figsize=(8,8))
    ax1s = axs.reshape(-1)
    ic=0

    npts = 1000

    for i in range(npts):
        print(i)
    
    

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
