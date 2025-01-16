# =======================================================================================
#
# For usage: $python RadiationUTestDriver.py --help
#
# This script runs unit tests on the two-stream functions.
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
import xml.etree.ElementTree as ET
import numpy as np
import matplotlib
import os
import sys
import getopt
import code
import time
import importlib
import csv
import ctypes
from ctypes import *
from operator import add
sys.path.append('../shared/py_src')
from PyF90Utils import c8, ci, cchar, c8_arr, ci_arr, ccharnb

font = {'family' : 'sans-serif',
        'weight' : 'normal',
        'size'   : 11}

matplotlib.rc('font', **font)


# Instantiate the F90 modules
f90_shr_obj = ctypes.CDLL('bld/WrapShrMod.o',mode=ctypes.RTLD_GLOBAL)
f90_mem_obj = ctypes.CDLL('bld/FatesRadiationMemMod.o',mode=ctypes.RTLD_GLOBAL)
f90_twostr_obj = ctypes.CDLL('bld/TwoStreamMLPEMod.o',mode=ctypes.RTLD_GLOBAL)
f90_wrap_obj = ctypes.CDLL('bld/RadiationWrapMod.o',mode=ctypes.RTLD_GLOBAL)


# Create aliases for the calls and define arguments if it helps with clarity
alloc_twostream_call =  f90_wrap_obj.__radiationwrapmod_MOD_initallocate
dealloc_twostream_call = f90_wrap_obj.__radiationwrapmod_MOD_dealloc
alloc_radparams_call = f90_twostr_obj.__twostreammlpemod_MOD_allocateradparams
set_radparams_call   = f90_wrap_obj.__radiationwrapmod_MOD_setradparam
set_radparams_call.argtypes = [POINTER(c_double),POINTER(c_int),POINTER(c_int),c_char_p,c_long]

param_prep_call = f90_twostr_obj.__twostreammlpemod_MOD_paramprep

setup_canopy_call = f90_wrap_obj.__radiationwrapmod_MOD_setupcanopy
setup_canopy_call.argtypes = [POINTER(c_int),POINTER(c_int),POINTER(c_int), \
                              POINTER(c_double),POINTER(c_double),POINTER(c_double)]

grndsnow_albedo_call = f90_wrap_obj.__radiationwrapmod_MOD_setgroundsnow
grndsnow_albedo_call.argtypes = [POINTER(c_int),POINTER(c_double),c_char_p,c_long]

canopy_prep_call = f90_wrap_obj.__radiationwrapmod_MOD_wrapcanopyprep
zenith_prep_call = f90_wrap_obj.__radiationwrapmod_MOD_wrapzenithprep
solver_call = f90_wrap_obj.__radiationwrapmod_MOD_wrapsolve
setdown_call = f90_wrap_obj.__radiationwrapmod_MOD_wrapsetdownwelling

getintens_call = f90_wrap_obj.__radiationwrapmod_MOD_wrapgetintensity
getabsrad_call = f90_wrap_obj.__radiationwrapmod_MOD_wrapgetabsrad
getparams_call = f90_wrap_obj.__radiationwrapmod_MOD_wrapgetparams
forceparam_call = f90_wrap_obj.__radiationwrapmod_MOD_wrapforceparams
forceparam_call.argtypes = [POINTER(c_int),POINTER(c_int),POINTER(c_int),POINTER(c_double),c_char_p,c_long]

leaf_rhonir = [0.46, 0.41, 0.39, 0.46, 0.41, 0.41, 0.46, 0.41, 0.41, 0.28, 0.28, 0.28 ]
leaf_rhovis = [0.11, 0.09, 0.08, 0.11, 0.08, 0.08, 0.11, 0.08, 0.08, 0.05, 0.05, 0.05 ]
leaf_taunir = [0.33, 0.32, 0.42, 0.33, 0.43, 0.43, 0.33, 0.43, 0.43, 0.4,  0.4,  0.4 ]
leaf_tauvis = [0.06, 0.04, 0.06, 0.06, 0.06, 0.06, 0.06, 0.06, 0.06, 0.05, 0.05, 0.05]
leaf_xl     = [0.32, 0.01, 0.01, 0.32, 0.2, 0.59, 0.32, 0.59, 0.59, -0.23, -0.23, -0.23]
leaf_clumping_index = [0.85, 0.85, 0.8, 0.85, 0.85, 0.9, 0.85, 0.9, 0.9, 0.75, 0.75, 0.75]
stem_rhonir = [0.49, 0.36, 0.36, 0.49, 0.49, 0.49, 0.49, 0.49, 0.49, 0.53, 0.53, 0.53]
stem_rhovis = [0.21, 0.12, 0.12, 0.21, 0.21, 0.21, 0.21, 0.21, 0.21, 0.31, 0.31, 0.31]
stem_taunir = [0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.25, 0.25, 0.25]
stem_tauvis = [0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.12, 0.12, 0.12]


visb = 1
nirb = 2

normalized_boundary = 1
absolute_boundary = 2

class elem_type:
    def __init__(self,n_vai):

        self.area = -9.0
        self.lai  = -9.0
        self.sai  = -9.0

        self.n_vai = n_vai
        self.avai = np.zeros([n_vai])
        self.r_dn = np.zeros([n_vai])
        self.r_up = np.zeros([n_vai])
        self.r_b  = np.zeros([n_vai])
        self.r_abs = np.zeros([n_vai])
        #self.sunfrac = np.zeros([n_vai])

class patch_type:
    def __init__(self,ground_albedo_diff,ground_albedo_beam):
        self.ground_albedo_beam = ground_albedo_diff
        self.ground_albedo_beam = ground_albedo_beam
        self.cohorts = []

        # uses the form:
        # patch.cohorts.append(cohort_type(n_vai,lai,sai))
                 
class cohort_type:
    def __init__(self,n_vai,area_frac,lai,sai,pft):

        self.n_vai = n_vai
        #self.avai = np.zeros([n_vai])
        dvai = (lai+sai)/n_vai
        self.avai = np.linspace(dvai,lai+sai,num=n_vai)
        self.rd_abs_leaf = np.zeros([n_vai])
        self.rb_abs_leaf = np.zeros([n_vai])
        self.r_abs_stem = np.zeros([n_vai])
        self.sunfrac = np.zeros([n_vai])
        self.pft = pft
        
def main(argv):

    # All tests will use 2 bands 1=vis, 2=nir

    # Initialize radiation parameters
    n_bands = 2
    n_pft   = 12

    iret = alloc_radparams_call(ci(n_pft),ci(n_bands))

    for ft in range(n_pft):

        pft=ft+1
        # rho (vis+nir)
        iret = set_radparams_call(c_double(leaf_rhovis[ft]),c_int(pft),c_int(visb),*ccharnb("rhol"))
        iret = set_radparams_call(c_double(leaf_rhonir[ft]),c_int(pft),c_int(nirb),*ccharnb("rhol"))
        iret = set_radparams_call(c_double(stem_rhovis[ft]),c_int(pft),c_int(visb),*ccharnb("rhos"))
        iret = set_radparams_call(c_double(stem_rhonir[ft]),c_int(pft),c_int(nirb),*ccharnb("rhos"))
        # tau (vis+nir)
        iret = set_radparams_call(c_double(leaf_tauvis[ft]),c_int(pft),c_int(visb),*ccharnb("taul"))
        iret = set_radparams_call(c_double(leaf_taunir[ft]),c_int(pft),c_int(nirb),*ccharnb("taul"))
        iret = set_radparams_call(c_double(stem_tauvis[ft]),c_int(pft),c_int(visb),*ccharnb("taus"))
        iret = set_radparams_call(c_double(stem_taunir[ft]),c_int(pft),c_int(nirb),*ccharnb("taus"))
        # orientations
        iret = set_radparams_call(c_double(leaf_xl[ft]),c_int(pft),c_int(0),*ccharnb("xl"))
        iret = set_radparams_call(c_double(leaf_clumping_index[ft]),c_int(pft),c_int(0),*ccharnb("clumping_index"))
        
    # Process the core 2Stream parameters from parameters in file
    iret = param_prep_call(ci(n_pft))

    if(False):
        TestCrash()
    
    if(False):
        ParallelElementPerturbDist()

    if(False):
        SunFracTests()

    if(True):
        SingleElementPerturbTest()

    if(False):
        SerialParallelCanopyTest()

    plt.show()

def TestCrash():

    # This is used to diagnose a specific failure.  This is probably
    # reconstructed from the output dump of a failed solve.

    xmlfile = "f45error_elements.xml"
    xmlroot = ET.parse(xmlfile).getroot()
    print("\nOpenend: "+xmlfile)
    
    cosz = float(xmlroot.find('cosz').text.strip())
    ib = int(xmlroot.find('band_id').text.strip())
    #elem              = xmlroot.find('time_control')

    # Iterate through canopy layers
    areas = []
    print("Loading Layers")
    for can in xmlroot.iter('can'):
        print("canopy layer: {}".format(int(can.attrib['id'].strip())))
        # Iterate through elements in each layer
        can_id = int(can.attrib['id'].strip())
        for elem in can.iter('elem'):
            elem_id = int(elem.attrib['id'].strip())
            textlist = elem.text.split(',')
            pft  = int(textlist[0].strip())
            lai  = float(textlist[1].strip())
            sai  = float(textlist[2].strip())
            area = float(textlist[3].strip())

            areas.append(area)

    code.interact(local=dict(globals(), **locals()))

            
def SerialParallelCanopyTest():


    # Lets first construct a bunch of cohorts, 5 cohorts
    # equal area, but folding by 2 in LAI

    cohort_lai  = np.array([0.25,0.5,1.0,2.0,4.0])
    cohort_area = np.array([0.2,0.2,0.2,0.2,0.2])
    n_cohorts = len(cohort_lai)
    
    sai_frac = 0.1
    
    pft = 1
    
    # Serial approach: 5 layers with veg and ghost
    n_col = 2
    n_layer = 5
    iret = alloc_twostream_call(ci(n_layer),ci(n_col))

    #class cohort_type:
    #def __init__(self,n_vai,area_frac,lai,sai,pft)
        
    # Five elements (cohorts), each take up 20% of the space
    area_frac = 0.2
    serialc = []
    serialc.append(cohort_type(100,area_frac,cohort_lai[0],cohort_lai[0]*sai_frac,pft))
    serialc.append(cohort_type(100,area_frac,cohort_lai[1],cohort_lai[1]*sai_frac,pft))
    serialc.append(cohort_type(100,area_frac,cohort_lai[2],cohort_lai[2]*sai_frac,pft))
    serialc.append(cohort_type(100,area_frac,cohort_lai[3],cohort_lai[3]*sai_frac,pft))
    serialc.append(cohort_type(100,area_frac,cohort_lai[4],cohort_lai[4]*sai_frac,pft))

    parallelc = []
    parallelc.append(cohort_type(100,area_frac,cohort_lai[0],cohort_lai[0]*sai_frac,pft))
    parallelc.append(cohort_type(100,area_frac,cohort_lai[1],cohort_lai[1]*sai_frac,pft))
    parallelc.append(cohort_type(100,area_frac,cohort_lai[2],cohort_lai[2]*sai_frac,pft))
    parallelc.append(cohort_type(100,area_frac,cohort_lai[3],cohort_lai[3]*sai_frac,pft))
    parallelc.append(cohort_type(100,area_frac,cohort_lai[4],cohort_lai[4]*sai_frac,pft))

    # Setup serial canopy "s_elems"
    
    s_elems = []
    #s_elems.append([])
    
    n_vai = 100



    dvai = 0.05
    for i in range(n_layer):
        s_elems.append([])
        # Serial Setup
        ican = i+1
        icol = 1
        area = np.sum(cohort_area[i:])
        if(i==0):
            lai = cohort_lai[i]
        else:
            lai = cohort_lai[i]-cohort_lai[i-1]
        
        sai  = lai*sai_frac

        n_vai = int((lai+sai)/dvai)
        s_elems[i].append(elem_type(n_vai))
        
        s_elems[i][-1].lai  = lai
        s_elems[i][-1].sai  = sai
        s_elems[i][-1].area = area
        s_elems[i][-1].avai = np.linspace(0,lai+sai,num=n_vai)
        iret = setup_canopy_call(c_int(ican),c_int(icol),c_int(pft),c_double(area),c_double(lai),c_double(sai))

        icol = 2
        area = 1-np.sum(cohort_area[i:])
        s_elems[i].append(elem_type(1))
        s_elems[i][-1].lai  = 0.0
        s_elems[i][-1].sai  = 0.0
        s_elems[i][-1].area = area
        lai  = 0.0
        sai  = 0.0
        air_pft = 0
        iret = setup_canopy_call(c_int(ican),c_int(icol),c_int(air_pft),c_double(area),c_double(lai),c_double(sai))
        
    # Decide on a band:
    ib = visb
    
    cd_r_beam = c_double(-9.0)
    cd_r_diff_up = c_double(-9.0)
    cd_r_diff_dn = c_double(-9.0)
    cd_kb = c_double(-9.0)
    cd_kd = c_double(-9.0)
    cd_om = c_double(-9.0)
    cd_betad = c_double(-9.0)
    cd_betab = c_double(-9.0)
    cd_rd_abs_leaf = c_double(-9.0)
    cd_rb_abs_leaf = c_double(-9.0)
    cd_r_abs_stem  = c_double(-9.0)
    cd_r_abs_snow  = c_double(-9.0)
    cd_leaf_sun_frac = c_double(-9.0)

    cd_albedo_beam = c_double(-9.0)
    cd_albedo_diff = c_double(-9.0)
    cd_canabs_beam = c_double(-9.0)
    cd_canabs_diff = c_double(-9.0)
    cd_ffbeam_beam = c_double(-9.0)
    cd_ffdiff_beam = c_double(-9.0)
    cd_ffdiff_diff = c_double(-9.0)
    
    
    R_beam = 1.
    R_diff = 1.
    cosz   = np.cos(0.0)

    ground_albedo_diff = 0.3
    ground_albedo_beam = 0.3
    frac_snow = 0.0
    
    iret = grndsnow_albedo_call(c_int(visb),c_double(ground_albedo_diff),*ccharnb('albedo_grnd_diff'))
    iret = grndsnow_albedo_call(c_int(visb),c_double(ground_albedo_beam),*ccharnb('albedo_grnd_beam'))
    iret = grndsnow_albedo_call(c_int(nirb),c_double(ground_albedo_diff),*ccharnb('albedo_grnd_diff'))
    iret = grndsnow_albedo_call(c_int(nirb),c_double(ground_albedo_beam),*ccharnb('albedo_grnd_beam'))
    iret = canopy_prep_call(c8(frac_snow))
    iret = zenith_prep_call(c8(cosz))
    iret = solver_call(ci(ib),ci(normalized_boundary),c8(1.0),c8(1.0), \
                       byref(cd_albedo_beam),byref(cd_albedo_diff), \
                       byref(cd_canabs_beam),byref(cd_canabs_diff), \
                       byref(cd_ffbeam_beam),byref(cd_ffdiff_beam),byref(cd_ffdiff_diff))
    iret = setdown_call(ci(ib),c8(R_beam),c8(R_diff))
    
    for i in range(n_layer):
        
        ican = i+1
        icol = 1
        for iv in range(s_elems[i][0].n_vai):
            iret = getintens_call(ci(ican),ci(icol),ci(ib),c8(s_elems[i][0].avai[iv]),byref(cd_r_diff_dn),byref(cd_r_diff_up),byref(cd_r_beam))
            s_elems[i][0].r_dn[iv] = cd_r_diff_dn.value
            s_elems[i][0].r_up[iv] = cd_r_diff_up.value
            s_elems[i][0].r_b[iv] = cd_r_beam.value
            if(iv>0):
                s_elems[i][0].r_abs[iv-1] = (s_elems[i][0].r_dn[iv]-s_elems[i][0].r_dn[iv-1]) + \
                    (s_elems[i][0].r_up[iv-1]-s_elems[i][0].r_up[iv]) + \
                    (s_elems[i][0].r_b[iv]-s_elems[i][0].r_b[iv-1])

        icol=2					 
        for iv in range(s_elems[i][1].n_vai):
            iret = getintens_call(ci(ican),ci(icol),ci(ib),c8(s_elems[i][1].avai[iv]),byref(cd_r_diff_dn),byref(cd_r_diff_up),byref(cd_r_beam))
            s_elems[i][1].r_dn[iv] = cd_r_diff_dn.value
            s_elems[i][1].r_up[iv] = cd_r_diff_up.value
            s_elems[i][1].r_b[iv] = cd_r_beam.value
            print('air: {} {} {}'.format(ican,icol,cd_r_beam.value))
            if(iv>0):
                s_elems[i][1].r_abs[iv-1] = (s_elems[i][1].r_dn[iv]-s_elems[i][1].r_dn[iv-1]) + \
                    (s_elems[i][1].r_up[iv-1]-s_elems[i][1].r_up[iv]) + \
                    (s_elems[i][1].r_b[iv]-s_elems[i][1].r_b[iv-1])

    # Lets get the absorbed radiation from the cohorts
    
    #class cohort_type:
    #def __init__(self,n_vai,lai,sai):
        #self.n_vai = n_vai
        ##self.avai = np.zeros([n_vai])
        #dvai = (lai+sai/n_vai)
        #self.avai = np.linspace(dvai,lai+sai,num=n_vai)
        #self.rabs_leaf = np.zeros([n_vai])
        #self.rabs_stem = np.zeros([n_vai])

    for i in range(len(serialc)):
        for iv in range(serialc[i].n_vai):

            vai_bot = serialc[i].avai[iv]

            ican = np.sum(serialc[i].avai[iv]>(cohort_lai*(1+sai_frac)))
            if(ican>0):
                vai_above = cohort_lai[ican-1]*(1+sai_frac)
            else:
                vai_above = 0.
                
            vai_bot = serialc[i].avai[iv]-vai_above
            if(iv==0):
                vai_top = 0
            else:
                vai_top = np.max([0,serialc[i].avai[iv-1]-vai_above])

            #print(i,iv,serialc[i].avai[iv],vai_above,vai_bot,vai_top,ican,cohort_lai*(1+sai_frac))
            icol = 1  # b/c 2 is air
            iret = getabsrad_call(ci(ican+1),ci(icol),ci(ib),c8(vai_top),c8(vai_bot), \
                                  byref(cd_rd_abs_leaf),byref(cd_rb_abs_leaf),byref(cd_r_abs_stem), \
                                  byref(cd_r_abs_snow),byref(cd_leaf_sun_frac))
            serialc[i].rd_abs_leaf[iv] = cd_rd_abs_leaf.value
            serialc[i].rb_abs_leaf[iv] = cd_rb_abs_leaf.value
            serialc[i].r_abs_stem[iv] = cd_r_abs_stem.value
            serialc[i].sunfrac[iv] = cd_leaf_sun_frac.value



    # Plot out absorbances and sun fractions in cohorts only
    # ---------------------------------------------
        
    
    
    max_rd_abs_leaf = 0
    max_rb_abs_leaf = 0
    max_r_abs_stem  = 0
    max_r_abs       = 0
    maxlai          = 0
    max_sunfrac     = 0
    for i in range(n_cohorts):
        max_rd_abs_leaf = np.max([max_rd_abs_leaf,np.max(serialc[i].rd_abs_leaf) ])
        max_rb_abs_leaf = np.max([max_rb_abs_leaf,np.max(serialc[i].rb_abs_leaf) ])
        max_r_abs_stem  = np.max([max_r_abs_stem,np.max(serialc[i].r_abs_stem) ])
        max_r_abs       = np.max([max_r_abs,np.max(serialc[i].r_abs_stem+serialc[i].rd_abs_leaf+serialc[i].rb_abs_leaf) ])
        maxlai          = np.max([maxlai,np.max(serialc[i].avai) ])
        max_sunfrac     = np.max([max_sunfrac,np.max(serialc[i].sunfrac)])

    fig, axs = plt.subplots(ncols=n_cohorts,nrows=1,figsize=(6,3))
    ax1s = axs.reshape(-1)
    
    y0   = 0.1
    xpad = 0.1
    dx   = (1.0-2*xpad)/float(n_cohorts)
    dy   = 0.8
        
    ic=0
    x0 = xpad
    for i in range(n_cohorts):

        ax = ax1s[ic]
        ap = ax.plot(serialc[i].rd_abs_leaf+serialc[i].rb_abs_leaf+serialc[i].r_abs_stem ,serialc[i].avai)
        ax.set_ylim([0,maxlai])
        ax.invert_yaxis()
        ax.set_xlabel('[W/m2]')
        ax.set_xlim([0,max_r_abs])

        ax.set_title('Cohort {}'.format(i+1))
        if(i==0):
            
            ax.set_ylabel('Absorbed Radiation\nVAI [m2/m2]')
        else:
            ax.set_yticklabels([])
            
        ax.grid(True)
        ax.set_position([x0,y0,dx,dy])
        x0 = x0+dx
        ic=ic+1

    fig, axs = plt.subplots(ncols=n_cohorts,nrows=1,figsize=(6,3))
    ax1s = axs.reshape(-1)
    
    y0   = 0.1
    xpad = 0.1
    dx   = (1.0-2*xpad)/float(n_cohorts)
    dy   = 0.8

    # Sun fractions
    ic=0
    x0 = xpad
    for i in range(n_cohorts):

        ax = ax1s[ic]
        ap = ax.plot(serialc[i].sunfrac ,serialc[i].avai)
        ax.set_ylim([0,maxlai])
        ax.invert_yaxis()
        ax.set_xlabel('[m2/m2]')
        ax.set_xlim([0,max_sunfrac])

        ax.set_title('Cohort {}'.format(i+1))
        if(i==0):
            
            ax.set_ylabel('Sunlit fraction of leaves [m2/m2]')
        else:
            ax.set_yticklabels([])
            
        ax.grid(True)
        ax.set_position([x0,y0,dx,dy])
        x0 = x0+dx
        ic=ic+1

    dealloc_twostream_call()


        
    if(False):
        PlotRadMaps(s_elems,0,'Beam Radiation [W/m2]')
        PlotRadMaps(s_elems,1,'Downwelling Diffuse Radiation [W/m2]')
        PlotRadMaps(s_elems,2,'Upwelling Diffuse Radiation [W/m2]')
        
    # Setup paralell canopy p_elems
    p_elems = []
    iret = alloc_twostream_call(ci(1),ci(n_cohorts))
    # Only one layer, so just one append
    p_elems.append([])
    for i in range(n_cohorts):

        icol = i+1
        ican = 1
        # Parallel
        
        p_elems[0].append(elem_type(n_vai))
        lai = cohort_lai[i]
        sai = sai_frac * cohort_lai[i]
        area = cohort_area[i]
        p_elems[0][-1].lai  = lai
        p_elems[0][-1].sai  = sai
        p_elems[0][-1].area = area
        p_elems[0][-1].avai = np.linspace(0,cohort_lai[i]*(1.+sai_frac),num=n_vai)
        iret = setup_canopy_call(c_int(ican),c_int(icol),c_int(pft),c_double(area),c_double(lai),c_double(sai))

    iret = grndsnow_albedo_call(c_int(visb),c_double(ground_albedo_diff),*ccharnb('albedo_grnd_diff'))
    iret = grndsnow_albedo_call(c_int(visb),c_double(ground_albedo_beam),*ccharnb('albedo_grnd_beam'))
    iret = grndsnow_albedo_call(c_int(nirb),c_double(ground_albedo_diff),*ccharnb('albedo_grnd_diff'))
    iret = grndsnow_albedo_call(c_int(nirb),c_double(ground_albedo_beam),*ccharnb('albedo_grnd_beam'))
    iret = canopy_prep_call(c8(frac_snow))
    iret = zenith_prep_call(c8(cosz))
    iret = solver_call(ci(ib),ci(normalized_boundary),c8(1.0),c8(1.0), \
                       byref(cd_albedo_beam),byref(cd_albedo_diff), \
                       byref(cd_canabs_beam),byref(cd_canabs_diff), \
                       byref(cd_ffbeam_beam),byref(cd_ffdiff_beam),byref(cd_ffdiff_diff))
    iret = setdown_call(ci(ib),c8(R_beam),c8(R_diff))

    ican = 1
    for i in range(n_cohorts):
       icol = i+1
       for iv in range(p_elems[0][i].n_vai):
          iret = getintens_call(ci(ican),ci(icol),ci(ib),c8(p_elems[0][i].avai[iv]),byref(cd_r_diff_dn),byref(cd_r_diff_up),byref(cd_r_beam))
          p_elems[0][i].r_dn[iv] = cd_r_diff_dn.value
          p_elems[0][i].r_up[iv] = cd_r_diff_up.value
          p_elems[0][i].r_b[iv] = cd_r_beam.value
          if(iv>0):
             p_elems[0][i].r_abs[iv-1] = (p_elems[0][i].r_dn[iv]-p_elems[0][i].r_dn[iv-1]) + \
                (p_elems[0][i].r_up[iv-1]-p_elems[0][i].r_up[iv]) + \
                (p_elems[0][i].r_b[iv]-p_elems[0][i].r_b[iv-1])

    dealloc_twostream_call()
    if(True):
        PlotRadMaps(p_elems,0,'Beam Radiation [W/m2]')
        PlotRadMaps(p_elems,1,'Downwelling Diffuse Radiation [W/m2]')
        PlotRadMaps(p_elems,2,'Upwelling Diffuse Radiation [W/m2]')
				 
def SunFracTests():


    n_col    = 1
    n_layer  = 1
    iret = alloc_twostream_call(ci(n_layer),ci(n_col))

    ican = 1 # Single canopy layer
    icol = 1 # Single PFT
    pft  = 1 # Use PFT number 1
    area = 1.0  # Assume only 90% of the ground is covered
    lai  = 5.0  # LAI
    sai  = 0.5  # SAI
    vai = lai+sai
    iret = setup_canopy_call(c_int(1),c_int(1),c_int(pft),c_double(area),c_double(lai),c_double(sai))

    # Decide on a band:
    ib = visb
    cd_r_beam = c_double(-9.0)
    cd_r_diff_up = c_double(-9.0)
    cd_r_diff_dn = c_double(-9.0)
    cd_kb = c_double(-9.0)
    cd_kd = c_double(-9.0)
    cd_om = c_double(-9.0)
    cd_betad = c_double(-9.0)
    cd_betab = c_double(-9.0)

    R_beam = 1.
    R_diff = 0.
    cosz   = np.cos(0.0)
    n_vai  = 200
    n_cosz = 100

    dv = vai/n_vai
    vai_a  = np.linspace(dv,vai,num=n_vai)
    cosz_a = np.linspace(0,1.0,num=n_cosz)
    kb_a   = np.zeros([n_cosz])
    lsf_a  = np.zeros([n_cosz,n_vai])
    rbeamsf_a =  np.zeros([n_cosz,n_vai])
    rbeam_a =  np.zeros([n_cosz,n_vai])
    
    

    cd_rd_abs_leaf = c_double(-9.0)
    cd_rb_abs_leaf = c_double(-9.0)
    cd_r_abs_stem  = c_double(-9.0)
    cd_r_abs_snow  = c_double(-9.0)
    cd_leaf_sun_frac = c_double(-9.0)
    cd_albedo_beam = c_double(-9.0)
    cd_albedo_diff = c_double(-9.0)
    cd_canabs_beam = c_double(-9.0)
    cd_canabs_diff = c_double(-9.0)
    cd_ffbeam_beam = c_double(-9.0)
    cd_ffdiff_beam = c_double(-9.0)
    cd_ffdiff_diff = c_double(-9.0)
    cd_r_diff_dn   = c_double(-9.0)
    cd_r_diff_up   = c_double(-9.0)
    cd_r_beam      = c_double(-9.0)
    
    ground_albedo_diff = 0.3
    ground_albedo_beam = 0.3
    frac_snow = 0.5
    
    iret = grndsnow_albedo_call(c_int(ib),c_double(ground_albedo_diff),*ccharnb('albedo_grnd_diff'))
    iret = grndsnow_albedo_call(c_int(ib),c_double(ground_albedo_beam),*ccharnb('albedo_grnd_beam'))

    iret = canopy_prep_call(c8(frac_snow))

    for ic,cosz in enumerate(cosz_a):
        iret = zenith_prep_call(c8(cosz))
    
        iret = solver_call(ci(ib),ci(normalized_boundary),c8(1.0),c8(1.0), \
                           byref(cd_albedo_beam),byref(cd_albedo_diff), \
                           byref(cd_canabs_beam),byref(cd_canabs_diff), \
                           byref(cd_ffbeam_beam),byref(cd_ffdiff_beam),byref(cd_ffdiff_diff))

        iret = setdown_call(ci(ib),c8(R_beam),c8(R_diff))
        
        iret = getparams_call(ci(ican),ci(icol),ci(ib),byref(cd_kb), \
                              byref(cd_kd),byref(cd_om),byref(cd_betad),byref(cd_betab))

        kb_a[ic] = cd_betab.value

        for iv in range(n_vai):

            if(iv==0):
                vai_top = 0.
            else:
                vai_top = vai_a[iv-1]

            vai_bot = vai_a[iv]
            
            
            iret = getabsrad_call(ci(ican),ci(icol),ci(ib),c8(vai_top),c8(vai_bot), \
                                  byref(cd_rd_abs_leaf),byref(cd_rb_abs_leaf),byref(cd_r_abs_stem), \
                                  byref(cd_r_abs_snow),byref(cd_leaf_sun_frac))

            iret = getintens_call(ci(ican),ci(icol),ci(ib),c8(vai_bot),byref(cd_r_diff_dn), \
                              byref(cd_r_diff_up),byref(cd_r_beam))
            
            lsf_a[ic,iv] = cd_leaf_sun_frac.value
            
            #sun_area = (vai_bot - vai_top)*cd_leaf_sun_frac.value/cd_kb.value
            sun_area = (vai_bot - vai_top)*cd_kb.value
            rbeam_a[ic,iv] = cd_r_beam.value
            
            if(iv==0):
                rbeamsf_a[ic,iv] = R_beam*(1.0 - sun_area)
                #print(rbeamsf_a[ic,iv],sun_area,vai_bot,vai_top,cd_leaf_sun_frac.value,vai_a[iv])
                #exit(0)
            else:
                rbeamsf_a[ic,iv] = rbeamsf_a[ic,iv-1]*(1.0 - sun_area)
                #print(rbeamsf_a[ic,iv])
            
    fig, axs = plt.subplots(ncols=2,nrows=2,figsize=(9,5))
    ax1s = axs.reshape(-1)

    ic0 = [2,25,50,99]
    
    for ia,ax in enumerate(ax1s):

        #Plot LSF profiles at 4 different cosz's

        ap = ax.plot(lsf_a[ic0[ia],:],vai_a[:],rbeam_a[ic0[ia],:],vai_a[:])
        ax.invert_yaxis()
        ax.set_title('cos(z) = {:.2f}'.format(cosz_a[ic0[ia]]))
        ax.set_xlabel('[Sun Fraction]')
        ax.set_xlim([0,1])
        ax.grid(True)
        if(ia<2):
            ax.set_xlabel('')
            ax.set_xticklabels([])
        #if(ia==0):
        #    ax.set_ylabel('Absorbed Radiation\nVAI [m2/m2]')
        #else:
        #    ax.set_yticklabels([])
    plt.tight_layout()

    fig2, axs = plt.subplots(ncols=2,nrows=2,figsize=(9,5))
    ax1s = axs.reshape(-1)

    ic0 = [2,25,50,99]
    
    for ia,ax in enumerate(ax1s):

        #Plot LSF profiles at 4 different cosz's

        ap = ax.plot(rbeam_a[ic0[ia],:],vai_a[:],rbeamsf_a[ic0[ia],:],vai_a[:])
        ax.invert_yaxis()
        ax.set_title('cos(z) = {:.2f}'.format(cosz_a[ic0[ia]]))
        ax.set_xlabel('[Beam Fraction]')
        ax.set_xlim([0,1])
        ax.grid(True)
        if(ia<2):
            ax.set_xlabel('')
            ax.set_xticklabels([])
        #if(ia==0):
        #    ax.set_ylabel('Absorbed Radiation\nVAI [m2/m2]')
        #else:
        #    ax.set_yticklabels([])
    plt.tight_layout()

    
    dealloc_twostream_call()
            
        
def ParallelElementPerturbDist():


    # Lets first construct a bunch of cohorts, 5 cohorts
    # equal area, but folding by 2 in LAI

    cohort_lai  = np.array([0.25,0.5,1.0,2.0,4.0])
    cohort_area = np.array([0.9,0.19,0.19,0.19,0.19])
    n_cohorts = len(cohort_lai)
    
    sai_frac = 0.1
    
    pft = 1
    
    # Serial approach: 5 layers with veg and ghost
    n_col = n_cohorts+1
    n_layer = 1
    iret = alloc_twostream_call(ci(n_layer),ci(n_col))

    for icol in range(n_col-1):
        iret = setup_canopy_call(c_int(1),c_int(icol+1),c_int(pft), \
                                 c_double(cohort_area[icol]),c_double(cohort_lai[icol]),c_double(cohort_lai[icol]*sai_frac))
        
    # Add the air element
    iret = setup_canopy_call(c_int(1),c_int(n_col),c_int(0),c_double(1.0-np.sum(cohort_area)),c_double(0.0),c_double(0.0))

    num_params = 9
    paramsets = []
    
    labels = ["clumping_index","leaf_rhonir","leaf_rhovis","leaf_taunir","leaf_tauvis", \
              "stem_rhonir","stem_rhovis","stem_taunir","stem_tauvis"]
    
    ic = 0
    with open('albedo_callib_param_vals.csv', newline='') as csvfile:
        
        reader = csv.reader(csvfile, delimiter=',')
        next(reader, None)
        nsets=0
        for irow, rowtext in enumerate(reader):
            ic=ic+1
            if(ic==num_params):
                ic=0
                nsets=nsets+1
                
    with open('albedo_callib_param_vals.csv', newline='') as csvfile:       
        paramset = np.zeros([num_params,nsets])
        reader = csv.reader(csvfile, delimiter=',')
        next(reader, None)
        ic=0
        iset=0
        for irow, rowtext in enumerate(reader):
            paramset[ic,iset] = float(rowtext[3])
            ic=ic+1
            if(ic==num_params):
                ic=0
                iset=iset+1


    fig1, axs = plt.subplots(3,3,figsize=(9,7))
    ax1s = axs.reshape(-1)

    for ip,ax in enumerate(ax1s):
    
        ap = ax.hist(paramset[ip,:])
        #ax1.set_ylabel('Integrated VAI [m2/m2]')
        ax.set_title(labels[ip])
        ax.grid(True)
        
    plt.tight_layout()
    plt.show()
    dealloc_twostream_call()

    
def SingleElementPerturbTest():

    
    # ===================================================================================
    # In this test, we have a canopy that is constructed from a single cohort
    # and therefore a single element. The cohort does not cover all of the ground
    # so their is an air element in parallel with the leaf/stem element.

    ground_albedo_diff = 0.1
    ground_albedo_beam = 0.1
    veg_frac_snow = 0.0
    
    patch = patch_type(ground_albedo_diff,ground_albedo_beam)

    # Vegetation cohort
    area_frac = 0.9
    lai = 2.0
    sai = 0.5
    pft = 1
    air_pft = 0
    patch.cohorts.append(cohort_type(100,area_frac,lai,sai,pft))

    # Open space (air)
    patch.cohorts.append(cohort_type(100,1.0-area_frac,0.,0.,air_pft))

    
    n_col    = 2
    n_layer  = 1
    iret = alloc_twostream_call(ci(n_layer),ci(n_col))

    ican = 1 # Single canopy layer
    icol = 1 # Single PFT
    pft  = 1 # Use PFT number 1
    vai = lai+sai
    iret = setup_canopy_call(c_int(1),c_int(1),c_int(pft),c_double(area_frac),c_double(lai),c_double(sai))
    iret = setup_canopy_call(c_int(1),c_int(2),c_int(0),c_double(1.0-area_frac),c_double(0.0),c_double(0.0))

    # Decide on a band:

    ib = visb
    
    cd_r_beam = c_double(-9.0)
    cd_r_diff_up = c_double(-9.0)
    cd_r_diff_dn = c_double(-9.0)
    cd_kb = c_double(-9.0)
    cd_kd = c_double(-9.0)
    cd_om = c_double(-9.0)
    cd_betad = c_double(-9.0)
    cd_betab = c_double(-9.0)
    

    # Make parameter pertubations, bump up 50%
    pp_dict = {}
    pp_dict['Kb'] = 1.5*0.66118239744  #74   #*1.5
    pp_dict['Kd'] = 1.5*0.9063246621781269  #*1.5
    pp_dict['om'] = 1.5*0.17819999999999997  #*1.5
    pp_dict['betab'] = 1.5*0.48253004714288084  #*1.5
    pp_dict['betad'] = 1.5*0.5999777777777778  #*1.5

    R_beam = 1.0
    R_diff = 1.0
    cosz   = np.cos(0.0)
    n_vai  = 100
    vai_a  = np.linspace(0,vai,num=n_vai)

    dv = vai/n_vai
    
    r_diff_up = np.zeros(n_vai)   
    r_diff_dn = np.zeros(n_vai)
    r_beam    = np.zeros(n_vai)

    drdv_diff_up = np.zeros(n_vai-1) # Delta
    drdv_diff_dn = np.zeros(n_vai-1) # Delta
    drdv_ubeam    = np.zeros(n_vai-1) # Delta
    drdv_dbeam    = np.zeros(n_vai-1) # Delta
    
    p_r_diff_up = np.zeros([n_vai,len(pp_dict)])
    p_r_diff_dn = np.zeros([n_vai,len(pp_dict)])
    p_r_beam    = np.zeros([n_vai,len(pp_dict)])
    p_drdv_diff_up = np.zeros([n_vai-1,len(pp_dict)])
    p_drdv_diff_dn = np.zeros([n_vai-1,len(pp_dict)])
    p_drdv_ubeam    = np.zeros([n_vai-1,len(pp_dict)])
    p_drdv_dbeam    = np.zeros([n_vai-1,len(pp_dict)])

    cd_albedo_beam = c_double(-9.0)
    cd_albedo_diff = c_double(-9.0)
    cd_canabs_beam = c_double(-9.0)
    cd_canabs_diff = c_double(-9.0)
    cd_ffbeam_beam = c_double(-9.0)
    cd_ffdiff_beam = c_double(-9.0)
    cd_ffdiff_diff = c_double(-9.0)
    
   
    iret = grndsnow_albedo_call(c_int(ib),c_double(ground_albedo_diff),*ccharnb('albedo_grnd_diff'))
    iret = grndsnow_albedo_call(c_int(ib),c_double(ground_albedo_beam),*ccharnb('albedo_grnd_beam'))
    iret = canopy_prep_call(c8(veg_frac_snow))
    iret = zenith_prep_call(c8(cosz))

    iret = solver_call(ci(ib),ci(normalized_boundary),c8(1.0),c8(1.0), \
                       byref(cd_albedo_beam),byref(cd_albedo_diff), \
                       byref(cd_canabs_beam),byref(cd_canabs_diff), \
                       byref(cd_ffbeam_beam),byref(cd_ffdiff_beam),byref(cd_ffdiff_diff))

    iret = setdown_call(ci(ib),c8(R_beam),c8(R_diff))
    
    iret = getparams_call(ci(ican),ci(icol),ci(ib),byref(cd_kb), \
                          byref(cd_kd),byref(cd_om),byref(cd_betad),byref(cd_betab))

    #print(cd_kb.value,cd_kd.value,cd_om.value,cd_betad.value,cd_betab.value)
    #exit(0)
    
    
    for iv in range(n_vai):
        iret = getintens_call(ci(ican),ci(icol),ci(ib),c8(vai_a[iv]),byref(cd_r_diff_dn), \
                              byref(cd_r_diff_up),byref(cd_r_beam))
        
        r_beam[iv] = cd_r_beam.value
        r_diff_up[iv] = cd_r_diff_up.value
        r_diff_dn[iv] = cd_r_diff_dn.value

        if(iv>0):
            drdv_ubeam[iv-1] = -cd_om.value*cd_betab.value*(r_beam[iv]-r_beam[iv-1])/dv
            drdv_dbeam[iv-1] = -cd_om.value*(1.-cd_betab.value)*(r_beam[iv]-r_beam[iv-1])/dv
            drdv_diff_dn[iv-1] = -(r_diff_dn[iv]-r_diff_dn[iv-1])/dv
            drdv_diff_up[iv-1] = (r_diff_up[iv]-r_diff_up[iv-1])/dv
            
    # Redo the scattering with perturbations
    i = -1
    for key,val in pp_dict.items():
        i=i+1
        iret = canopy_prep_call(c8(veg_frac_snow))
        iret = zenith_prep_call(c8(cosz))
        iret = forceparam_call(c_int(ican),c_int(icol),ci(ib),c_double(val),*ccharnb(key))

        iret = solver_call(ci(ib),ci(normalized_boundary),c8(1.0),c8(1.0), \
                       byref(cd_albedo_beam),byref(cd_albedo_diff), \
                       byref(cd_canabs_beam),byref(cd_canabs_diff), \
                       byref(cd_ffbeam_beam),byref(cd_ffdiff_beam),byref(cd_ffdiff_diff))

        iret = setdown_call(ci(ib),c8(R_beam),c8(R_diff))
        
        for iv in range(n_vai):
            iret = getintens_call(ci(ican),ci(icol),ci(ib),c8(vai_a[iv]),byref(cd_r_diff_dn),byref(cd_r_diff_up),byref(cd_r_beam))

            #print(iv,i,cd_r_beam.value)
            p_r_beam[iv,i] = cd_r_beam.value
            p_r_diff_up[iv,i] = cd_r_diff_up.value
            p_r_diff_dn[iv,i] = cd_r_diff_dn.value

            if(iv>0):
                p_drdv_ubeam[iv-1] = -cd_om.value*cd_betab.value*(p_r_beam[iv]-p_r_beam[iv-1])/dv
                p_drdv_dbeam[iv-1] = -cd_om.value*(1.-cd_betab.value)*(p_r_beam[iv]-p_r_beam[iv-1])/dv
                p_drdv_diff_dn[iv-1] = -(p_r_diff_dn[iv]-p_r_diff_dn[iv-1])/dv
                p_drdv_diff_up[iv-1] = (p_r_diff_up[iv]-p_r_diff_up[iv-1])/dv


        fig1, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(6.5,5.5))

        ap = ax1.plot(r_beam,vai_a,p_r_beam[:,i],vai_a)
        first_color = ap[0].get_color()
        last_color = ap[-1].get_color()
        ax1.invert_yaxis()
        ax1.set_xlabel('')
        ax1.set_ylabel('Integrated VAI [m2/m2]')
        ax1.set_title('Beam Intensity [W/m2]')
        ax1.grid(True)

        ax2.plot(r_diff_dn,vai_a,p_r_diff_dn[:,i],vai_a)
        ax2.invert_yaxis()
        ax2.set_xlabel('')
        ax2.set_yticklabels('')
        ax2.set_ylabel('')
        ax2.set_title('Down Diffuse Intensity [W/m2] ')
        ax2.grid(True)
        
        ax3.plot(r_diff_up,vai_a,p_r_diff_up[:,i],vai_a)
        ax3.invert_yaxis()
        ax3.set_xlabel('')
        ax3.set_ylabel('Integrated VAI [m2/m2]')
        ax3.set_title('Up Diffuse Intensity [W/m2]')
        ax3.grid(True)

        ax4.axis("off")
        ax4.set_axis_off()

        if(ib==visb):
            band_name = "Visible"
        elif(ib==nirb):
            band_name = "Near Infrared"
        else:
            print("Unknown band")
            exit(2)
        
            
        param_str = r"""In-element Scattering Profiles
Broad band: {0}
$cos(\phi) = ${1:.2f}
$K_b = ${2:.2f}
$K_d = ${3:.2f} 
$\omega = ${4:.2f} 
$\beta_b = ${5:.2f}
$\beta_d = ${6:.2f}
$\alpha_{{gd}} = ${7:.2f}
$\alpha_{{gb}} = ${8:.2f}""".format(band_name,cosz,cd_kb.value,cd_kd.value,cd_om.value,cd_betab.value,cd_betad.value,ground_albedo_diff,ground_albedo_beam)    
        ax4.text(0.1, 0.5, param_str, horizontalalignment='left', \
                 verticalalignment='center', transform=ax4.transAxes,backgroundcolor=[1.0,1.0,1.0],fontsize=11,color=first_color)
        ax4.text(0.5,0.5,r"{0}={1:.2f}".format(key,val),color=last_color)
        plt.subplots_adjust(wspace=0.1, hspace=0.25)
        plt.tight_layout()
        plt.show()


    dealloc_twostream_call()
    

# Plotting Functions

    
def PlotRadMaps(elems,rtype,plt_title):    
        
    fig, ax = plt.subplots(ncols=1,nrows=1,figsize=(5,5))

    cmap = mpl.cm.Reds
    #code.interact(local=dict(globals(), **locals()))
    n_layer = len(elems)
    total_vai = 0
    for i in range(n_layer):
       max_vai = 0.
       for j in range(len(elems[i])):	 
          max_vai = np.max([max_vai,elems[i][j].lai+elems[i][j].sai])
       total_vai = total_vai + max_vai
    
    ax.set_ylim([0,total_vai])
    
    total_vai = 0
    rect = []
    rcolor = []
    for i in range(n_layer):  
		  #
        max_vai = 0.
        area_off = 0.
        for j in range(len(elems[i])):
            max_vai = np.max([max_vai,elems[i][j].lai+elems[i][j].sai])
        for j in range(len(elems[i])):
            for iv in range(elems[i][j].n_vai):
               if(rtype==0):
                  rel_intense = np.max([0,elems[i][j].r_b[iv]])
               elif(rtype==1):
                  rel_intense = np.max([0,elems[i][j].r_dn[iv]])
               elif(rtype==2):
                  rel_intense = np.max([0,elems[i][j].r_up[iv]])
               
               if(iv==0):
                  yoff = total_vai
                  dvai = elems[i][j].avai[iv]
               else:
                  yoff = total_vai+elems[i][j].avai[iv-1]
                  dvai = elems[i][j].avai[iv]-elems[i][j].avai[iv-1]
               rect.append(mpl.patches.Rectangle((area_off,yoff),elems[i][j].area,dvai))
               rcolor.append(rel_intense)
            area_off = area_off + elems[i][j].area

        total_vai = total_vai + max_vai
				
        # Air
        #rel_intense = np.max([0,np.min([1.,elems[1][i].r_dn[0]/R_diff])])
        #rel_intense = np.max([0,elems[1][i].r_dn[0]])
        #if(rtype==0):
        #    rel_intense = np.max([0,elems[1][i].r_b[0]])
        #elif(rtype==1):
        #    rel_intense = np.max([0,elems[1][i].r_dn[0]])
        #elif(rtype==2):
        #    rel_intense = np.max([0,elems[1][i].r_up[0]])

        
        #rect.append(mpl.patches.Rectangle((elems[0][i].area,total_vai),(1.-elems[0][i].area),(elems[0][i].lai+elems[0][i].sai))) #,color = [rel_intense,0.5,0.5]))
        #rcolor.append(rel_intense)
        
    p = mpl.collections.PatchCollection(rect,cmap = cmap,alpha = 1.0)
    p.set_array(rcolor)
    im = ax.add_collection(p)

        
        #code.interact(local=dict(globals(), **locals()))
        
    ax.invert_yaxis()
    ax.set_ylabel('Integrated Vegetated Area Index')
    ax.set_xlabel('Ground Area Fraction')
    ax.set_title(plt_title)  #)
    plt.colorbar(im)
    plt.show()
    
def PlotRadLines():
    
    fig, axs = plt.subplots(ncols=2,nrows=n_layer,figsize=(8,8))
    ax1s = axs.reshape(-1)
    ic=0
    y0   = 0.9
    ypad = 0.1
    dy = (y0-ypad)/n_layer
    xpad = 0.1
    xwid = 1-2*xpad
    
    for i in range(n_layer):

        ax = ax1s[ic]
        ap = ax.plot(elems[0][i].r_dn,elems[0][i].avai)
        ax.set_ylim([np.min(elems[0][i].avai),np.max(elems[0][i].avai)])
        ax.invert_yaxis()
        ax.set_xlabel('')
        ax.set_xlim([0,R_diff])
        ax.set_ylabel('VAI [m2/m2]')
        if(i==0):
            ax.set_title('Diffuse Down Intensity [W/m2]')
        if(i!=n_layer-1):
            ax.set_xticklabels([])
        ax.grid(True)
        y0 = y0-dy
        x0 = xpad
        dx = 0.4
        #dx = elems[0][i].area*(1-2*xpad)
        ax.set_position([x0,y0,dx,dy])
        ic=ic+1
        
        ax = ax1s[ic]
        ap = ax.plot([elems[1][i].r_dn[0],elems[1][i].r_dn[-1]],[0,1])
        ax.invert_yaxis()
        ax.set_xlabel('')
        ax.set_xlim([0,R_diff])
        if(i==0):
            ax.set_title('Diffuse Down Intensity [W/m2]')
        if(i!=n_layer-1):
            ax.set_xticklabels([])
        ax.set_ylabel('')
        ax.set_yticklabels([])
        ax.set_yticks([])
        ax.set_ylim([0,1])
        x0 = xpad+dx
        dx=0.4
        #dx = elems[1][i].area*(1-2*xpad)
        ax.set_position([x0,y0,dx,dy])
        ax.grid(True)
        ic=ic+1


    
# =======================================================================================
# This is the actual call to main

if __name__ == "__main__":
    main(sys.argv)
