import numpy as np
import os
import sys
import getopt
import code  # For development: code.interact(local=locals())
import time

day_per_year = 365.0

class pft_bc_type:

    def __init__(self):

        # Initialize a dictionary of parameters for any pft
        self.pft_bc_dic = {}


def DailyCFromUnitGPPAR(leaf_area,AGB):
    
    # -----------------------------------------------------------------------------------
    # This routine estimates Net Daily Carbon Gains (GPP-AR) by estimating
    # a mean canopy GPP per leaf area per year, and by estimating
    # a mean autotrophic respiration per kilogram per year, both from literature.
    # Thus to scale to a plant, the plant's leaf area and total biomass are needed.
    #
    # THese numbers are taken from Chambers et al. 2004
    # from ZF2 Manaus Brazil
    # -----------------------------------------------------------------------------------
        
    kg_per_Mg  = 1000.0
    m2_per_ha  = 10000.0
        
    site_AGB   = 151.35 # MgC/ha
    site_NPP   = 9.0    # MgC/ha/yr
    site_AR    = 21.0   # MgC/ha/yr
    site_LAI   = 4.7    # m2/m2

    #site_Rleaf = 9.8    # MgC/ha/yr
    #site_Rwood = 4.2    # MgC/ha/yr
    #site_Rroot = 5.5    # MgC/ha/yr


    GPP_per_larea_yr = kg_per_Mg * (site_NPP + site_AR) / \
                       site_LAI / m2_per_ha
    AR_per_kg_yr   = kg_per_Mg * site_AR / site_AGB / \
                     m2_per_ha

    GPP = 100.8*GPP_per_larea_yr * leaf_area / day_per_year
    AR  = AR_per_kg_yr * AGB     / day_per_year

    NetDailyC = GPP - AR
        
    return NetDailyC


def DailyCFromCArea(presc_npp_p1,c_area,phen_type,leaf_status):
    
    # -----------------------------------------------------------------------------------
    # This method was provided by Charlie Koven via is inferences from the PPA
    # literature.  Here, net daily carbon [kg] is based on one of two excluding
    # parmaters (NPP per crown area per year), for plants that are either in
    # the upper canopy (access to sunlight) or in the understory (low sunlight)
    #
    # c_area, footprint of the crown area [m2].
    # presc_npp_p1, npp generated per crown area         [kgC/m2/yr]
    # -----------------------------------------------------------------------------------

    if( (phen_type == 1) or (leaf_status ==2)):
        NetDailyC = presc_npp_p1 * c_area / day_per_year        
    else:
        NetDailyC = 0.0

    return NetDailyC


def DailyCNPFromCArea(presc_npp_p1,presc_nflux_p1, \
                      presc_pflux_p1,c_area,phen_type,leaf_status):
    
    # -----------------------------------------------------------------------------------
    # This method was provided by Charlie Koven via is inferences from the PPA
    # literature.  Here, net daily carbon [kg] is based on one of two excluding
    # parmaters (NPP per crown area per year), for plants that are either in
    # the upper canopy (access to sunlight) or in the understory (low sunlight)
    #
    # c_area, footprint of the crown area [m2].
    # presc_npp_canopy, npp generated per crown area in canopy         [kgC/m2/yr]
    # presc_npp_understory, npp generated per crown area in understory [kgC/m2/yr]
    # presc_nflux_p1, Nitrogen flux per crown area                     [kgN/m2/yr]
    # presc_pflux_p1, Phosphorous flux per crown area                  [kgP/m2/yr]
    # -----------------------------------------------------------------------------------

    if( (phen_type == 1) or (leaf_status ==2)):
        NetDailyC = presc_npp_p1 * c_area / day_per_year        
        NetDailyN = presc_nflux_p1 * c_area / day_per_year
        NetDailyP = presc_pflux_p1 * c_area / day_per_year
    else:
        NetDailyC = 0.0
        NetDailyN = 0.0
        NetDailyP = 0.0

    return NetDailyC, NetDailyN, NetDailyP


def DailyCNPFromStorageSinWave(doy,store_c,presc_npp_p1, \
                               presc_nflux_p1,presc_pflux_p1,c_area,presc_npp_amp, \
                               phen_type, leaf_status):


    # This method is supposed to simulate a seasonal cycle of NPP
    # In some cases we pass negative daily carbon gain to the allocation model
    # however, we have to be careful to not make negative gains larger
    # than available storage in those cases.  This is not necessarily the most
    # realistic model, but its important to test that the parteh algorithms can handle
    # these stressfull negative gain conditions.

    doy0=0.0

    sin_func = np.sin( (doy-doy0)/366.0 * 2.0 * np.pi )
    
    #if (sin_func>0.0):
    #    NetDailyC = sin_func * presc_npp_p1 * c_area / day_per_year
    #else:
    #    NetDailyC = -np.minimum( -neg_store_frac * sin_func * presc_npp_p1* c_area / day_per_year, 0.98* np.float(store_c))

    NetDailyC = (presc_npp_amp * sin_func * presc_npp_p1 + presc_npp_p1) * c_area/day_per_year

    # This is a fail-safe, for large negatives, cant be larger than storage
    
    if (NetDailyC < 0.0):
        NetDailyC = -np.minimum(-NetDailyC,0.98* np.float(store_c))

    #print("sin_func: {}, NetDailyC: {}, store_c: {}, c_area :{}".format(sin_func,NetDailyC,store_c,c_area))

    if( (phen_type == 1) or (leaf_status ==2)):
        NetDailyN = presc_nflux_p1 * c_area / day_per_year
        NetDailyP = presc_pflux_p1 * c_area / day_per_year
    else:
        NetDailyN = 0.0
        NetDailyP = 0.0
        NetDailyC = 0.0
    
    return NetDailyC, NetDailyN, NetDailyP


def DeciduousPhenology(doy, target_leaf_c, store_c, phen_type):

    # Time leaf-on with rising NPP
    leaf_on_doy  = np.int(366.0 * 0.01)

    leaf_off_doy = np.int(366.0 * 0.55)

    if ( doy==leaf_on_doy):
        flush_c = np.minimum(store_c,target_leaf_c * 0.5)
    else:
        flush_c = 0.0

    if ( doy==leaf_off_doy):
        drop_frac_c = 1.0
    else:
        drop_frac_c = 0.0

    if(doy>=leaf_on_doy and doy<leaf_off_doy):
        leaf_status = 2          # Leaves are on
    else:
        leaf_status = 1          # Leaves are off

    if(phen_type==1):
        flush_c     = 0.0
        drop_frac_c = 0.0
        leaf_status = 2

    return  flush_c, drop_frac_c, leaf_status



