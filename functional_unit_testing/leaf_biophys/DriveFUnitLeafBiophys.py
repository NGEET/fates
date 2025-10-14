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


# Plot support routines
# ========================================================================

def EvalVJKByTemp(pft,fates_leaf_vcmax25top,leaf_c3psn,pdf):

    # Plot out vcmax, jmax and kp as a function of temperature
    # Assumes canopy top position, and is dependent on the
    # PFT base rate
        
    # Leaf temperature ranges [C]
    leaf_tempc_min = -20.0
    leaf_tempc_max = 50.0
    leaf_tempc_n = 100
    leaf_tempc_vec = np.linspace(leaf_tempc_min,leaf_tempc_max,num=leaf_tempc_n)

    vcmax_f = c_double(-9)
    jmax_f  = c_double(-9)
    kp_f    = c_double(-9)
    gs0_f   = c_double(-9)
    gs1_f   = c_double(-9)
    gs2_f   = c_double(-9)
    
    jmax25_top,kp25_top =  GetJmaxKp25Top(fates_leaf_vcmax25top)
    vcmax = np.zeros([leaf_tempc_n])
    jmax  = np.zeros([leaf_tempc_n])
    kp    = np.zeros([leaf_tempc_n])

    for it, leaf_tempc in enumerate(leaf_tempc_vec):

        leaf_tempk = leaf_tempc + tfrz_1atm

        iret = f90_biophysrate_sub(ci(pft+1), c8(fates_leaf_vcmax25top), \
                                   c8(jmax25_top), c8(kp25_top), \
                                   c8(nscaler_top), c8(leaf_tempk), c8(dayl_factor_full), \
                                   c8(t_growth_kum),c8(t_home_kum),c8(btran_nolimit), \
                                   byref(vcmax_f), byref(jmax_f), byref(kp_f),byref(gs0_f),byref(gs1_f),byref(gs2_f))

        vcmax[it] = vcmax_f.value
        jmax[it]  = jmax_f.value
        kp[it]    = kp_f.value

            
    if(leaf_c3psn == 0):
        fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(8.5,7.5))
    else:
        fig, (ax1,ax2) = plt.subplots(1,2,figsize=(8.5,5.5))
        
    ax1.plot(leaf_tempc_vec,vcmax)
    ax1.set_ylabel('Vcmax [umol/m2/s]')
    ax1.set_title('PFT: {}'.format(pft+1))
    ax1.grid(True)

    ax2.plot(leaf_tempc_vec,jmax)
    ax2.set_xlabel('Leaf Temperature [C]')
    ax2.set_ylabel('Jmax [umol/m2/s]')
    ax2.set_title('PFT: {}'.format(pft+1))
    ax2.grid(True)

    if(leaf_c3psn == 0):
        ax3.plot(leaf_tempc_vec,kp/umol_per_mol)
        ax3.set_xlabel('Leaf Temperature [C]')
        ax3.set_ylabel('Kp [mol/m2/s]')
        ax3.grid(True)
        ax3.set_xlabel('Leaf Temperature [C]')
        ax4.axis("off")
    else:
        ax1.set_xlabel('Leaf Temperature [C]')

    if(pdf):
        pdf.savefig(fig)
        plt.close(fig) 
        
            

# ========================================================================

    
def LinePlotY3dM1(ax,x1,x2,x3,y3d,str_x2,str_x3,add_labels):

    # This takes a 3d array and plots that array over
    # the first of its 3 dimensions.  It does this
    # by evaluating 3 percentiles in each of the other
    # two dimensions, such as the 0, 50% and 99%, and plotting
    # the nine different combinations. The plot axis
    # is passed in, and there is minimal stylizing of the
    # axis and plot formatting, which is left to be done
    # outside of this, and in the calling script
    
    # Find the indices on the second and third dimensions
    # that are the 10th,50th and 90th percentiles
    ix2_10 = int(len(x2)*0.05)
    ix2_90 = int(len(x2)*0.90)
    ix3_10 = int(len(x3)*0.05)
    ix3_90 = int(len(x3)*0.90)

    redish = [0.8,0.3,0.2]
    bluish = [0.3,0.3,0.8]

    ix2 = ix2_10;ix3=ix3_10;
    fmt_str = '%s = %5.2f %s = %5.2f'%(str_x2,x2[ix2],str_x3,x3[ix3])
    if(add_labels):
        ap1 = ax.plot(x1,y3d[:,ix2,ix3],linestyle='dashed',color=bluish,label=fmt_str)
    else:
        ap1 = ax.plot(x1,y3d[:,ix2,ix3],linestyle='dashed',color=bluish)
        
    ix2 = ix2_10;ix3=ix3_90;
    fmt_str = '%s = %5.2f %s = %5.2f'%(str_x2,x2[ix2],str_x3,x3[ix3])
    if(add_labels):
        ap1 = ax.plot(x1,y3d[:,ix2,ix3],linestyle='solid',color=bluish,label=fmt_str)
    else:
        ap1 = ax.plot(x1,y3d[:,ix2,ix3],linestyle='solid',color=bluish)
        
    ix2 = ix2_90;ix3=ix3_10;
    fmt_str = '%s = %5.2f %s = %5.2f'%(str_x2,x2[ix2],str_x3,x3[ix3])
    if(add_labels):
        ap1 = ax.plot(x1,y3d[:,ix2,ix3],linestyle='dashed',color=redish,label=fmt_str)
    else:
        ap1 = ax.plot(x1,y3d[:,ix2,ix3],linestyle='dashed',color=redish)
        
    ix2 = ix2_90;ix3=ix3_90;
    fmt_str = '%s = %5.2f %s = %5.2f'%(str_x2,x2[ix2],str_x3,x3[ix3])
    if(add_labels):
        ap1 = ax.plot(x1,y3d[:,ix2,ix3],linestyle='solid',color=redish,label=fmt_str)
    else:
        ap1 = ax.plot(x1,y3d[:,ix2,ix3],linestyle='solid',color=redish)
        
    ax.grid(True)
    
    

# =======================================================================================
def TestCiTol(fates_leaf_vcmax25top,leaf_c3psn,fates_stoich_nitr,fates_leaf_slatop,fates_maintresp_leaf_model,pdf):

   
    
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
    
    pft_n = len(fates_leaf_vcmax25top)
    
    gb_min = 0.5e6
    gb_max = 3.e6
    gb_n = 5
    gb_vec = np.linspace(gb_min,gb_max,gb_n)

    ci_tol_n = 5
    ci_tol_min = -3
    ci_tol_max = 1
    ci_tol_vec = np.logspace(ci_tol_min, ci_tol_max, num=ci_tol_n, base=10.0)

    # Leaf temperature ranges [C]
    leaf_tempc_min = -30.0
    leaf_tempc_max = 50.0
    leaf_tempc_n = 15
    leaf_tempc_vec = np.linspace(leaf_tempc_min,leaf_tempc_max,num=leaf_tempc_n)

    # Relative Humidity Ranges
    rh_max = 0.95
    rh_min = 0.10
    rh_n   = 15
    rh_vec = np.linspace(rh_min,rh_max,num=rh_n)
    
    # Absorbed PAR ranges [W/m2]
    par_abs_min = 0
    par_abs_max = 300
    par_abs_n  = 15
    par_abs_vec = np.linspace(par_abs_min,par_abs_max,num=par_abs_n)
    
    agross = np.zeros([leaf_tempc_n,rh_n,par_abs_n,gb_n,pft_n,ci_tol_n])
    gstoma = np.zeros([leaf_tempc_n,rh_n,par_abs_n,gb_n,pft_n,ci_tol_n])
    solves = np.zeros([leaf_tempc_n,rh_n,par_abs_n,gb_n,pft_n,ci_tol_n])
    

    for ipft in range(pft_n):

        jmax25_top,kp25_top =  GetJmaxKp25Top(fates_leaf_vcmax25top[ipft])
        # Leaf Nitrogen Concentration at the top
        lnc_top  = fates_stoich_nitr[ipft]/fates_leaf_slatop[ipft]
    
        for it, leaf_tempc in enumerate(leaf_tempc_vec):

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
        
            iret = f90_biophysrate_sub(ci(ipft+1), c8(fates_leaf_vcmax25top[ipft]), \
                                       c8(jmax25_top), c8(kp25_top), \
                                       c8(nscaler_top), c8(leaf_tempk), c8(dayl_factor_full), \
                                       c8(t_growth_kum),c8(t_home_kum),c8(btran_nolimit), \
                                       byref(vcmax_f), byref(jmax_f), byref(kp_f),byref(gs0_f),byref(gs1_f),byref(gs2_f))
        
            # Leaf Maintenance Respiration (temp and pft dependent)
            if(fates_maintresp_leaf_model==1):
                iret = f90_lmr_ryan_sub(c8(lnc_top),c8(nscaler_top), ci(ipft+1), c8(leaf_tempk), byref(lmr_f))
            elif(fates_maintresp_leaf_model==2):
                iret = f90_lmr_atkin_sub(c8(lnc_top),c8(rdark_scaler),c8(leaf_tempk),c8(atkin_mean_leaf_tempk),byref(lmr_f) )
            else:
                print('unknown leaf respiration model')
                exit(1)
            
            for ip, par_abs in enumerate(par_abs_vec):

                for ir, rh in enumerate(rh_vec):

                    vpress = rh * veg_es_f.value
                
                    for ig, gb_umol in enumerate(gb_vec):
                    
                        for ic, ci_tol in enumerate(ci_tol_vec):
                        
                            iret = f90_leaflayerphoto_sub(c8(par_abs),  \
                                                          ci(ipft+1),   \
                                                          c8(vcmax_f.value),   \
                                                          c8(jmax_f.value),    \
                                                          c8(kp_f.value),      \
                                                          c8(gs0_f.value),      \
                                                          c8(gs1_f.value),      \
                                                          c8(gs2_f.value),      \
                                                          c8(leaf_tempk), \
                                                          c8(can_press_1atm), \
                                                          c8(co2_ppress_400ppm), \
                                                          c8(o2_ppress_209kppm), \
                                                          c8(veg_es_f.value), \
                                                          c8(gb_umol), \
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
                        
                            agross[it,ir,ip,ig,ipft,ic] = agross_f.value
                            gstoma[it,ir,ip,ig,ipft,ic] = gstoma_f.value*1.e-6
                            solves[it,ir,ip,ig,ipft,ic] = float(solve_iter_f.value)

                            
    mpl.rcParams['lines.markersize'] = 2                
    fig10,axs = plt.subplots(4,3,figsize=(8.5,9.5))

    axs = axs.reshape(-1)

    agmax = np.max(agross)
    gsmax = np.max(gstoma)
    svmax = np.max(solves)

    dsmax = -10
    dsmin = 10

    for ic, ci_tol in enumerate(ci_tol_vec):
        dsmax = np.max([dsmax,np.max( solves[:,:,:,:,:,0].reshape(-1)-solves[:,:,:,:,:,ic].reshape(-1)) ])
        dsmin = np.min([dsmax,np.min( solves[:,:,:,:,:,0].reshape(-1)-solves[:,:,:,:,:,ic].reshape(-1)) ])

    
    hbins = np.linspace(dsmin-0.5,dsmax+0.5,int(dsmax-dsmin+2.0))
    hbinc = np.linspace(dsmin,dsmax,int(dsmax-dsmin+1))
    
    ax0 = axs[0]
    ax0.scatter( agross[:,:,:,:,:,0].reshape(-1),agross[:,:,:,:,:,1].reshape(-1),marker='.')
    ax0.set_xlim([0,agmax])
    ax0.set_ylim([0,agmax])
    ax0.set_title('Ag [umol/m2/s]')
    ax0.set_ylabel('Ci_tol = %7.3f Pa'%(ci_tol_vec[1]))
    ax0.set_xticklabels('')
    ax0.grid(True)

    ax1 = axs[1]
    ax1.scatter( gstoma[:,:,:,:,:,0].reshape(-1),gstoma[:,:,:,:,:,1].reshape(-1),marker='.')
    ax1.set_xlim([0,gsmax])
    ax1.set_ylim([0,gsmax])
    ax1.set_title('gs [mol/m2/s]')
    ax1.set_xticklabels('')
    ax1.grid(True)

    ax2 = axs[2]
    delhist = solves[:,:,:,:,:,0].reshape(-1)-solves[:,:,:,:,:,1].reshape(-1)
    ax2.hist(delhist,bins=hbins)
    ax2.set_xticks(hbinc)
    ax2.set_title('Solver Count Reduction')
    ax2.grid(True)
    ax2.set_xticklabels('')
    ax2.set_yticklabels('')
    
    ax3 = axs[3]
    ax3.scatter( agross[:,:,:,:,:,0].reshape(-1),agross[:,:,:,:,:,2].reshape(-1),marker='.')
    ax3.set_xlim([0,agmax])
    ax3.set_ylim([0,agmax])
    ax3.set_ylabel('Ci_tol = %7.3f Pa'%(ci_tol_vec[2]))
    ax3.grid(True)
    ax3.set_xticklabels('')
    
    ax4 = axs[4]
    ax4.scatter( gstoma[:,:,:,:,:,0].reshape(-1),gstoma[:,:,:,:,:,2].reshape(-1),marker='.')
    ax4.set_xlim([0,gsmax])
    ax4.set_ylim([0,gsmax])
    ax4.grid(True)
    ax4.set_xticklabels('')
    
    ax5 = axs[5]
    delhist = solves[:,:,:,:,:,0].reshape(-1)-solves[:,:,:,:,:,2].reshape(-1)
    ax5.hist(delhist,bins=hbins)
    ax5.grid(True)
    ax5.set_xticks(hbinc)
    ax5.set_xticklabels('')
    ax5.set_yticklabels('')
    
    ax6 = axs[6]
    ax6.scatter( agross[:,:,:,:,:,0].reshape(-1),agross[:,:,:,:,:,3].reshape(-1),marker='.')
    ax6.set_xlim([0,agmax])
    ax6.set_ylim([0,agmax])
    ax6.set_ylabel('Ci_tol = %7.3f Pa'%(ci_tol_vec[3]))
    ax6.grid(True)
    ax6.set_xticklabels('')
    
    ax7 = axs[7]
    ax7.scatter( gstoma[:,:,:,:,:,0].reshape(-1),gstoma[:,:,:,:,:,3].reshape(-1),marker='.')
    ax7.set_xlim([0,gsmax])
    ax7.set_ylim([0,gsmax])
    ax7.grid(True)
    ax7.set_xticklabels('')
    
    ax8 = axs[8]
    delhist = solves[:,:,:,:,:,0].reshape(-1)-solves[:,:,:,:,:,3].reshape(-1)
    ax8.hist(delhist,bins=hbins)
    ax8.grid(True)
    ax8.set_xticks(hbinc)
    ax8.set_xticklabels('')
    ax8.set_yticklabels('')
    
    ax9 = axs[9]
    ax9.scatter( agross[:,:,:,:,:,0].reshape(-1),agross[:,:,:,:,:,4].reshape(-1),marker='.')
    ax9.set_xlim([0,agmax])
    ax9.set_ylim([0,agmax])
    ax9.set_ylabel('Ci_tol = %7.3f Pa'%(ci_tol_vec[4]))
    ax9.set_xlabel('Ci_tol = %7.3f Pa'%(ci_tol_vec[0]))
    ax9.grid(True)

    ax10 = axs[10]
    ax10.scatter( gstoma[:,:,:,:,:,0].reshape(-1),gstoma[:,:,:,:,:,4].reshape(-1),marker='.')
    ax10.set_xlim([0,gsmax])
    ax10.set_ylim([0,gsmax])
    ax10.grid(True)
    ax10.set_xlabel('Ci_tol = %7.3f Pa'%(ci_tol_vec[0]))

    ax11 = axs[11]
    delhist = solves[:,:,:,:,:,0].reshape(-1)-solves[:,:,:,:,:,4].reshape(-1)
    ax11.hist(delhist,bins=hbins)
    ax11.grid(True)
    ax11.set_xticks(hbinc)
    ax11.set_yticklabels('')
    
    plt.subplots_adjust(wspace=0.25,hspace=0.0)

    if(pdf):
        pdf.savefig(fig10)
    else:
        plt.show()


def main(argv):

    parser = argparse.ArgumentParser(description='Parse command line arguments to this script.')
    parser.add_argument('--fin', dest='xmlfile', type=str, help="XML control file, required.", required=True)
    parser.add_argument('--fout', dest='pdffile', type=str, help="PDF output file, not required..", required=False)
    args = parser.parse_args()

    
    # Load the xml control file
    # -----------------------------------------------------------------------------------
    xmlroot = et.parse(args.xmlfile).getroot()

    # Determine if we are generating a pdf
    # -----------------------------------------------------------------------------------
    try:
        pdf_file = args.pdffile
        pdf = PdfPages(pdf_file)
    except:
        pdf = None
    
    numpft = int(xmlroot.find('numpft').text.strip())

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
                
    # Dump parameters
    if(dump_parameters):
        iret = f90_dump_param_sub()

    
    # Read in non-fortran parameters from the xml file
    leafn_vert_scaler_coeff1 = []
    leafn_vert_scaler_coeff2 = []
    fates_leaf_vcmax25top    = []
    fates_stoich_nitr = []
    fates_leaf_slatop = []
    
    print('Reading non-fortran pft parameters')
    py_pft_root = xmlroot.find('py_params').find('pft_dim')
    for param in py_pft_root.iter('param'):
        for pft in range(numpft):
            if (param.attrib['name']=='fates_leafn_vert_scaler_coeff1'):
                leafn_vert_scaler_coeff1.append(np.float64(param.text.split(',')[pft]))
            if (param.attrib['name']=='fates_leafn_vert_scaler_coeff2'):
                leafn_vert_scaler_coeff2.append(np.float64(param.text.split(',')[pft]))
            if (param.attrib['name']=='fates_leaf_vcmax25top'):
                fates_leaf_vcmax25top.append(np.float64(param.text.split(',')[pft]))
            if (param.attrib['name']=='fates_stoich_nitr'):
                fates_stoich_nitr.append(np.float64(param.text.split(',')[pft]))
            if (param.attrib['name']=='fates_leaf_slatop'):
                fates_leaf_slatop.append(np.float64(param.text.split(',')[pft]))
                
    print('Reading non-fortran scalar parameters')
    py_scalar_root = xmlroot.find('py_params').find('scalar_dim')
    for param in py_scalar_root.iter('param'):
        if (param.attrib['name']=='fates_maintresp_leaf_model'):
            fates_maintresp_leaf_model = int(param.text.split(',')[0])


    # Look at CI tolerances
    if(do_test_citol):
        TestCiTol(fates_leaf_vcmax25top,leaf_c3psn, \
                  fates_stoich_nitr,fates_leaf_slatop,fates_maintresp_leaf_model,pdf)

        
    # Leaf temperature ranges [C]
    leaf_tempc_min = -30.0
    leaf_tempc_max = 50.0
    leaf_tempc_n = 50
    leaf_tempc_vec = np.linspace(leaf_tempc_min,leaf_tempc_max,num=leaf_tempc_n)

    # Relative Humidity Ranges
    rh_max = 0.99
    rh_min = 0.01
    rh_n   = 10
    rh_vec = np.linspace(rh_min,rh_max,num=rh_n)

    
    
    # Absorbed PAR ranges [W/m2]
    par_abs_min = 1.0
    par_abs_max = 500
    par_abs_n  = 10
    par_abs_vec = np.linspace(par_abs_min,par_abs_max,num=par_abs_n)

    # Boundary Conductance ranges [umol/m2/s]
    gb_min =  500000.0            # Lower limit imposed by CLM/ELM
    gb_max = 2500000.0            # 50% larger than  Roughly largestthe largest values seen at BCI (which are 2.5mol/m2/s)
    gb_n  = 5
    gb_vec = np.linspace(gb_min,gb_max,num=gb_n)

    # btran ranges

    btran_n   = 50
    btran_min = -2
    btran_max = 0
    btran_vec = np.logspace(btran_min,btran_max,num=btran_n)
    
    # These variables are mostly dependent on leaf temperature,
    # and weakly dependent on atmospheric pressure
    veg_qsat_vec = np.zeros([leaf_tempc_n])
    veg_esat_vec = np.zeros([leaf_tempc_n])
    mm_kco2_vec = np.zeros([leaf_tempc_n])
    mm_ko2_vec = np.zeros([leaf_tempc_n])
    co2_cpoint_vec = np.zeros([leaf_tempc_n])

    
    # Generic boundary layer conductance from Kimura et al ~1.2 cm/s
    # reasonable ranges for green houses are 0-3 cm/s
    # perhaps storms would be >10cm/s?
    # Convert to molar form using 1 standard atm at 25C
    # Units:  umol/m2/s
    # 1 mol/m2/s was roughly the middle of the distribution
    # when generating conductances at BCI using a static canopy
    # and local driver data.
    
    bl_cond_const = 1.e6

    # Lets look at canopy top
    # kn = DecayCoeffVcmax(currentCohort%vcmax25top, &
    #                     EDPftvarcon_inst%maintresp_leaf_vert_scaler_coeff1(ft), &
    #                     EDPftvarcon_inst%maintresp_leaf_vert_scaler_coeff2(ft))
    # rdark_scaler = exp(-kn * cumulative_lai)

    rdark_scaler = 1.0

    # kn = DecayCoeffVcmax(currentCohort%vcmax25top, &
    #                      prt_params%leafn_vert_scaler_coeff1(ft), &
    #                      prt_params%leafn_vert_scaler_coeff2(ft))

    nscaler = 1.0

    # Set convergence tolerance
    #ci_tol = 2.*can_press_1atm
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
    
    print('Prepping Canopy Gas Parameters')
    
    for it, leaf_tempc in enumerate(leaf_tempc_vec):

        leaf_tempk = leaf_tempc + tfrz_1atm
            
        iret = f90_qsat_sub(c8(leaf_tempk),c8(can_press_1atm), \
                            byref(veg_qs_f),byref(veg_es_f), \
                            byref(qsdt_dummy_f),byref(esdt_dummy_f))

        
        veg_qsat_vec[it] = veg_qs_f.value
        veg_esat_vec[it] = veg_es_f.value
        
        iret = f90_cangas_sub(c8(can_press_1atm), \
                              c8(o2_ppress_209kppm), \
                              c8(leaf_tempk), \
                              byref(mm_kco2_f), \
                              byref(mm_ko2_f), \
                              byref(co2_cpoint_f))

        mm_kco2_vec[it] = mm_kco2_f.value
        mm_ko2_vec[it] = mm_ko2_f.value
        co2_cpoint_vec[it] = co2_cpoint_f.value

    
        

    print('\n')
    print('Experiment 1: Evaluating Photosynthesis Equations by pft/Tl/RH/PR')

    elapsed_time = 0.0
    for pft in range(numpft):

        if(do_evalvjkbytemp):
            print('\n')
            print('Experiment 1: Evaluating Vcmax,Jmax,Kp by Temperature')
            EvalVJKByTemp(pft,fates_leaf_vcmax25top[pft],leaf_c3psn[pft],pdf)

        
        print('Evaluating PFT {}'.format(pft+1))
        
        jmax25_top,kp25_top =  GetJmaxKp25Top(fates_leaf_vcmax25top[pft])

        vcmax = np.zeros([leaf_tempc_n])
        jmax  = np.zeros([leaf_tempc_n])
        kp    = np.zeros([leaf_tempc_n])
        lmr    = np.zeros([leaf_tempc_n])
        agross = np.zeros([leaf_tempc_n,rh_n,par_abs_n,gb_n,btran_n])
        gstoma = np.zeros([leaf_tempc_n,rh_n,par_abs_n,gb_n,btran_n])
        co2_interc = np.zeros([leaf_tempc_n,rh_n,par_abs_n,gb_n,btran_n])
        anet   = np.zeros([leaf_tempc_n,rh_n,par_abs_n,gb_n])
        ac     = np.zeros([leaf_tempc_n,rh_n,par_abs_n,gb_n])
        aj     = np.zeros([leaf_tempc_n,rh_n,par_abs_n,gb_n])
        ap     = np.zeros([leaf_tempc_n,rh_n,par_abs_n,gb_n])
        iters  = np.zeros([leaf_tempc_n,rh_n,par_abs_n,gb_n])
        
        # When calling component limitations exclusively
        # using an approximation of interstitial co2 as
        # 0.7*canopy_co2
        
        ac2    = np.zeros([leaf_tempc_n])
        aj2    = np.zeros([leaf_tempc_n,par_abs_n])
        ap2    = np.zeros([leaf_tempc_n])
        
        # Leaf Nitrogen Concentration at the top
        lnc_top  = fates_stoich_nitr[pft]/fates_leaf_slatop[pft]
        
        for it, leaf_tempc in enumerate(leaf_tempc_vec):

            leaf_tempk = leaf_tempc + tfrz_1atm
            
            iret = f90_biophysrate_sub(ci(pft+1), c8(fates_leaf_vcmax25top[pft]), \
                                       c8(jmax25_top), c8(kp25_top), \
                                       c8(nscaler), c8(leaf_tempk), c8(dayl_factor_full), \
                                       c8(t_growth_kum),c8(t_home_kum),c8(btran_nolimit), \
                                       byref(vcmax_f), byref(jmax_f), byref(kp_f), byref(gs0_f), byref(gs1_f), byref(gs2_f))
            
            vcmax[it] = vcmax_f.value
            jmax[it]  = jmax_f.value
            kp[it]    = kp_f.value

            if(leaf_c3psn[pft] == 1):
                ap2[it] = 0.0
            else:
                ap2[it] = f90_agross_pepc4(c8(co2_ppress_400ppm),c8(kp[it]),c8(can_press_1atm))
                
            # Leaf Maintenance Respiration (temp and pft dependent)
            if(fates_maintresp_leaf_model==1):
                iret = f90_lmr_ryan_sub(c8(lnc_top),c8(nscaler), ci(pft+1), c8(leaf_tempk), byref(lmr_f))
            elif(fates_maintresp_leaf_model==2):
                iret = f90_lmr_atkin_sub(c8(lnc_top),c8(rdark_scaler),c8(leaf_tempk),c8(atkin_mean_leaf_tempk),byref(lmr_f) )
            else:
                print('unknown leaf respiration model')
                exit(1)

            lmr[it] = lmr_f.value

            if(leaf_c3psn[pft] == 1):
                
                ac2[it] = f90_agross_rubiscoc3(c8(vcmax[it]),c8(ci_ppress_static),c8(o2_ppress_209kppm), \
                                             c8(co2_cpoint_vec[it]),c8(mm_kco2_vec[it]),c8(mm_ko2_vec[it]))
            else:
                ac2[it] = vcmax[it]

                 
            for ip, par_abs in enumerate(par_abs_vec):

                if(leaf_c3psn[pft] == 1):
                    aj2[it,ip] = f90_agross_rubpc3(c8(par_abs),c8(jmax[it]),c8(ci_ppress_static),c8(co2_cpoint_vec[it]))
                else:
                    aj2[it,ip] = f90_agross_rubpc4(c8(par_abs))
                    
                for ir, rh in enumerate(rh_vec):

                    vpress = rh * veg_esat_vec[it]

                    for ig, gb in enumerate(gb_vec):

                        for ib, btran in enumerate(btran_vec):

                            iret = f90_biophysrate_sub(ci(pft+1), c8(fates_leaf_vcmax25top[pft]), \
                                                       c8(jmax25_top), c8(kp25_top), \
                                                       c8(nscaler), c8(leaf_tempk), c8(dayl_factor_full), \
                                                       c8(t_growth_kum),c8(t_home_kum),c8(btran), \
                                                       byref(vcmax_f), byref(jmax_f), byref(kp_f), byref(gs0_f), byref(gs1_f), byref(gs2_f))
                        
                            iret = f90_leaflayerphoto_sub(c8(par_abs), \
                                                          ci(pft+1),   \
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
                                                          c8(veg_esat_vec[it]), \
                                                          c8(gb), \
                                                          c8(vpress), \
                                                          c8(mm_kco2_vec[it]), \
                                                          c8(mm_ko2_vec[it]), \
                                                          c8(co2_cpoint_vec[it]), \
                                                          c8(lmr[it]), \
                                                          c8(ci_tol), \
                                                          byref(agross_f), \
                                                          byref(gstoma_f), \
                                                          byref(anet_f), \
                                                          byref(c13_f), \
                                                          byref(co2_interc_f), \
                                                          byref(solve_iter_f) )

                            # Call the medlyn solve to get timing info
                            iret = f90_qsat_sub(c8(leaf_tempk),c8(can_press_1atm), \
                                                byref(veg_qs_f),byref(veg_es_f), \
                                                byref(qsdt_dummy_f),byref(esdt_dummy_f))
                            
                           
                            
                            agross[it,ir,ip,ig,ib] = agross_f.value
                            gstoma[it,ir,ip,ig,ib] = gstoma_f.value
                            anet[it,ir,ip,ig] = anet_f.value
                            ac[it,ir,ip,ig] = ac_f.value
                            aj[it,ir,ip,ig] = aj_f.value
                            ap[it,ir,ip,ig] = ap_f.value
                            co2_interc[it,ir,ip,ig] = co2_interc_f.value
                            iters[it,ir,ip,ig] = solve_iter_f.value
                            time0 = time.process_time()
                            
                            iret = f90_gs_medlyn(c8(anet_f.value),c8(veg_es_f.value),c8(vpress),c8(gs0_f.value),c8(gs1_f.value), \
                                                 c8(gs2_f.value),c8(co2_ppress_400ppm),c8(can_press_1atm),c8(gb),byref(gstoma_f))
                            elapsed_time = elapsed_time + (time.process_time() - time0)


        print("\nElapsed Time for Conductance: {}\n".format(elapsed_time))
        
        # Plot out component gross assimilation rates
        # by temperature with constant Ci, and by Ci with
        # constant temperature
            
        fig2,ax1 = plt.subplots(1,1,figsize=(6.5,5.5))

        ax1.plot(leaf_tempc_vec,ac2,label='Ac')
        ix2_10 = int(par_abs_n*0.1)
        ix2_50 = int(par_abs_n*0.5)
        ix2_90 = int(par_abs_n*0.9)
        ax1.plot(leaf_tempc_vec,aj2[:,ix2_10],color=[0.5,0.5,0.5],linestyle='dotted',label='Aj apar=%4.1f'%(par_abs_vec[ix2_10]))
        ax1.plot(leaf_tempc_vec,aj2[:,ix2_50],color=[0.5,0.5,0.5],linestyle='dashed',label='Aj apar=%4.1f'%(par_abs_vec[ix2_50]))
        ax1.plot(leaf_tempc_vec,aj2[:,ix2_90],color=[0.5,0.5,0.5],linestyle='solid',label='Aj apar=%4.1f'%(par_abs_vec[ix2_90]))
        ax1.plot(leaf_tempc_vec,ap2[:],color='orange',linestyle='solid',label='Ap')
        ax1.plot(leaf_tempc_vec, lmr, color=[0.7,0.5,0.3],linestyle='solid',label='Rdark')
        ax1.set_ylabel('[umol/m2/s]')
        ax1.set_xlabel('Leaf Temperature [C]')
        ax1.set_title('PFT: %3i, Vcmax25: %4.1f, Jmax25: %4.1f, Ci: %4.1f'%(pft+1,fates_leaf_vcmax25top[pft],jmax25_top,ci_ppress_static))
        ax1.grid(True)
        fig2.legend(loc='upper left')
        if(pdf):
            pdf.savefig(fig2)
            plt.close(fig2)
        else:
            plt.savefig('images/rates_pft%2.2i.png'%(pft))
            
            
        # Lets plot metrics by temperature, using the
        # 10th, and 90th percentiles of both RH and PAR
        # Agross, Anet, Gstoma, Ac, Aj, Ap
        
        fig5, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(8.,7.))
        gb_id = 4
        bt_id = -1
        LinePlotY3dM1(ax1,leaf_tempc_vec,rh_vec,par_abs_vec, \
                      agross[:,:,:,gb_id,bt_id].reshape([leaf_tempc_n,rh_n,par_abs_n]),'RH','APAR',True)
        ax1.set_ylabel('Agross \n [umol/m2/s]')
        ax1.set_xticklabels([])
        LinePlotY3dM1(ax2,leaf_tempc_vec,rh_vec,par_abs_vec,\
                      gstoma[:,:,:,gb_id,bt_id].reshape([leaf_tempc_n,rh_n,par_abs_n])*1.e-6,'RH','APAR',False)
        ax2.set_ylabel('Gs \n [mol/m2/s]')
        ax2.yaxis.set_label_position("right")
        ax2.yaxis.tick_right()

        LinePlotY3dM1(ax3,leaf_tempc_vec,rh_vec,par_abs_vec,\
                      co2_interc[:,:,:,gb_id,bt_id].reshape([leaf_tempc_n,rh_n,par_abs_n]),'RH','APAR',False)
        ax3.set_ylabel('Ci [Pa]')
        ax3.axhline(y=co2_ppress_400ppm,color=[0.3,0.3,0.3])
        ax3.plot(leaf_tempc_vec,co2_cpoint_vec)
        
        ax3.set_xlabel('Leaf Temperature [C]')
        ax1.set_xlabel('Leaf Temperature [C]')


        ax4.text(0.25,0.6,'PFT: %2i \ng_b = %4.2f [mol/m2/s]'%(pft+1,gb_vec[gb_id]*1.e-6))
        ax4.axis('off')
        
        plt.tight_layout()
        plt.subplots_adjust(wspace=0.02, hspace=0.03)
        #fig3.legend(loc='lower right',labelspacing = 0.2, fontsize=12)
        fig5.legend(loc='lower right', bbox_to_anchor=(0.87, 0.06),labelspacing = 0.2, fontsize=12, ncol=1) #, fancybox=True, shadow=True)
        if(pdf):
            pdf.savefig(fig5)
            plt.close(fig5)
        else:
            plt.savefig('images/AgGsCI_temp_v2_pft%2.2i.png'%(pft))
        
        
        # Metrics across the boundary layer conductivity gradient
        
        # Metrics across the humidity Gradient  (fixed conductance)
        fig4, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(8.,7.))
        gb_id = 3
        bt_id = -1
        LinePlotY3dM1(ax1,rh_vec,leaf_tempc_vec,par_abs_vec, \
                      np.transpose(agross,(1,0,2,3,4 ))[:,:,:,gb_id,bt_id].reshape([rh_n,leaf_tempc_n,par_abs_n]),'Tveg','APAR',True)
        
        ax1.set_ylabel('Agross \n [umol/m2/s]')
        ax1.set_xticklabels([])
        LinePlotY3dM1(ax2,rh_vec,leaf_tempc_vec,par_abs_vec,\
                      np.transpose(gstoma,(1,0,2,3,4 ))[:,:,:,gb_id,bt_id].reshape([rh_n,leaf_tempc_n,par_abs_n])*1.e-6,'Tveg','APAR',False)
        ax2.set_ylabel('Gs \n [mol/m2/s]')
        ax2.yaxis.set_label_position("right")
        ax2.yaxis.tick_right()

        LinePlotY3dM1(ax3,rh_vec,leaf_tempc_vec,par_abs_vec,\
                      np.transpose(co2_interc,(1,0,2,3,4 ))[:,:,:,gb_id,bt_id].reshape([rh_n,leaf_tempc_n,par_abs_n]),'Tveg','APAR',False)
        ax3.set_ylabel('Ci [Pa]')
        ax3.axhline(y=co2_ppress_400ppm,color=[0.3,0.3,0.3])

        ax3.set_xlabel('Relative Humidity [%]')
        ax1.set_xlabel('Relative Humidity [%]')


        ax4.text(0.25,0.6,'PFT: %2i \ng_b = %4.2f [mol/m2/s]'%(pft+1,gb_vec[gb_id]*1.e-6))
        ax4.axis('off')
        
        plt.tight_layout()
        plt.subplots_adjust(wspace=0.02, hspace=0.03)
        fig4.legend(loc='lower right', bbox_to_anchor=(0.87, 0.06),labelspacing = 0.2, fontsize=12, ncol=1) #, fancybox=True, shadow=True)
        if(pdf):
            pdf.savefig(fig4)
            plt.close(fig4)
        else:
            plt.savefig('images/AgGsCI_temp_v3_pft%2.2i.png'%(pft))



        # Look at the btran gradient

        fig44, ((ax1,ax2)) = plt.subplots(1,2,figsize=(8.,5.))
        gb_id = 3
        rh_id = np.argmin(np.abs(rh_vec-1.0))  # ID at 100% humidity
        par200_id = np.argmin(np.abs(par_abs_vec-200.))
        tv30_id = np.argmin(np.abs(leaf_tempc_vec-30.))
        par50_id = np.argmin(np.abs(par_abs_vec-50.))
        tv5_id = np.argmin(np.abs(leaf_tempc_vec-5.))

        ax1.plot(btran_vec,agross[tv30_id,rh_id,par200_id,gb_id,:],label='PAR=%5.1f W Tv=%3.1f C'%(par_abs_vec[par200_id],leaf_tempc_vec[tv30_id]),color = 'r',linestyle = '-')
        ax1.plot(btran_vec,agross[tv30_id,rh_id,par50_id,gb_id,:],label='PAR=%5.1f W Tv=%3.1f C'%(par_abs_vec[par50_id],leaf_tempc_vec[tv30_id]),color = 'r',linestyle = '--')
        ax1.plot(btran_vec,agross[tv5_id,rh_id,par200_id,gb_id,:],label='PAR=%5.1f W Tv=%3.1f C'%(par_abs_vec[par200_id],leaf_tempc_vec[tv5_id]),color = 'b',linestyle = '-')
        ax1.plot(btran_vec,agross[tv5_id,rh_id,par50_id,gb_id,:],label='PAR=%5.1f W Tv=%3.1f C'%(par_abs_vec[par50_id],leaf_tempc_vec[tv5_id]),color = 'b',linestyle = '--')
        ax1.set_title('Ag [umol/m2/s]')
        ax1.set_xlabel('BTRAN')
        ax1.grid('on')

        ax2.plot(btran_vec,gstoma[tv30_id,rh_id,par200_id,gb_id,:]*1.e-6,color = 'r',linestyle = '-')
        ax2.plot(btran_vec,gstoma[tv30_id,rh_id,par50_id,gb_id,:]*1.e-6,color = 'r',linestyle = '--')
        ax2.plot(btran_vec,gstoma[tv5_id,rh_id,par200_id,gb_id,:]*1.e-6,color = 'b',linestyle = '-')
        ax2.plot(btran_vec,gstoma[tv5_id,rh_id,par50_id,gb_id,:]*1.e-6,color = 'b',linestyle = '--')
        ax2.set_title('gs [mol/m2/s]')
        ax2.set_xlabel('BTRAN')
        ax2.grid('on')
        fig44.legend(loc='lower right', bbox_to_anchor=(0.87, 0.56),labelspacing = 0.2, fontsize=12, ncol=1)


            


        # Show for each pft
        if(not pdf):
            plt.show()

            
    if(pdf):
        pdf.close()
        
    print('Deallocating parameter space')
    iret = f90_dealloc_leaf_param_sub()
    
    print('Functional Unit Testing Complete')
    exit(0)

    
# =======================================================================================
# This is the actual call to main

if __name__ == "__main__":
    main(sys.argv)
