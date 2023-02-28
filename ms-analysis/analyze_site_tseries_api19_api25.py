import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.dates as mdates
import datetime
import sys
import code  # For development: code.interact(local=locals())  code.interact(local=dict(globals(), **locals()))
import argparse
import math
from scipy.io import netcdf as nc
import cftime as cf
import os
import nc_time_axis
from scipy import interpolate
import pandas as pd
from mpl_toolkits.axes_grid1 import make_axes_locatable


xr.set_options(enable_cftimeindex=True)

# Use a window of this number of indices wide when performing rolling averages
roll_w_size = 16*12    # 1.5 years
#roll_w_size = 1

#years = [1600,2100]
years = [1,500]
#years=[2000,2020]
m2_per_ha = 10000.0
g_per_Mg  = 1000000.0
g_per_kg  = 1000.0
kg_per_g  = 0.001
m2_per_cm2 = 1.0/(100.0*100.0)
ha_per_m2 = 1.e-4
sec_per_day = 86400
day_per_sec = 1.0/86400.0
day_per_year = 365
sec_per_year = 86400*365

dbhcoord = 1
depthcoord = 2

# This routine evaluates the time signatures
# and makes sure that the data is annual. It also
# converts the time-stamps into the nc_time_axis format.
# Further, if we have a break in the time (which is typical
# with spin-ups), it identifies if that is true, and records
# the index


def monthly_to_annual(array):
    """ calculate annual mena from monthly data, using unequal month lengths fros noleap calendar.  
    originally written by Keith Lindsay."""
    mon_day  = xr.DataArray(np.array([31.,28.,31.,30.,31.,30.,31.,31.,30.,31.,30.,31.]), dims=['month'])
    mon_wgt  = mon_day/mon_day.sum()
    return (array.rolling(time=12, center=False) # rolling
            .construct("month") # construct the array
            .isel(time=slice(11, None, 12)) # slice so that the first element is [1..12], second is [13..24]
            .dot(mon_wgt, dims=["month"]))

def PrepNCTime(hist_ds):
    
    cdmo = [0, 31, 59, 90, 120, 151, 182, 212, 243, 273, 304, 334]
        
    yrs = np.floor(hist_ds.mcdate.data/10000.0)
    #mos = np.floor((hist_ds.mcdate.data-yrs*10000.0)/100.0)
    #dys = np.floor(hist_ds.mcdate.data-(yrs*10000.0 + mos*100.0))
    mos = 6*np.ones(yrs.shape)
    dys = np.zeros(yrs.shape)

    
    #code.interact(local=dict(globals(), **locals()))

    decyear = np.asarray([ (yrs[i] + float(cdmo[int(mos[i])-1] + dys[i])/365) for i in range(len(yrs))])

    # Lets find if there are any large breaks in the data
    yrdiff = decyear[1:] - decyear[:-1]
    maxdiff = np.max(yrdiff)
    meddiff = np.median(yrdiff)
    maxid   = np.argmax(yrdiff)+1

    if( np.abs(meddiff-1.0) < 0.001 ):
        # Annual, do nothing
        print('Identified dataset as annually averaged, no upscaling required')
        hist_dsy = hist_ds
    elif( np.abs(1./meddiff-12.0) < 0.01):
        print('Identified dataset as monthly averaged, please run upscale utility')
        exit()

    split_time=False
    maxid=None
    if(maxdiff>10.0*meddiff):
        split_time=True
        maxid = np.argmax(yrdiff)+1
        
    return decyear, split_time, maxid




# ========================================================================================
#                                        Main
# ========================================================================================

def main(argv):


    # Check versions, because this matters

    print(xr.__version__)
    print(pd.__version__)

    #    code.interact(local=dict(globals(), **locals()))
    
    parser = argparse.ArgumentParser(description='Parse command line arguments to this script.')
    parser.add_argument('--input_hists', dest='histfiles', type=str, help="comma delimited list of files", required=True)
    parser.add_argument('--labels',dest='labels',type=str,help="comma delimted list of labels",required=False)
    parser.add_argument('--lat',dest='lat',type=float,required=False)
    parser.add_argument('--lon',dest='lon',type=float,required=False)

    
    args = parser.parse_args()

    
    # Find all files with the prefixes, loop through and pre-process datasets

    
    hists_ds = LoadHists(args.histfiles)

    n_ds = len(hists_ds)

    if(args.labels is None):
        letters=['A','B','C','D','E','F','G']
        if(n_ds>len(letters)):
            print('Dang lots of datasets, add more label letters')
            exit(2)
        labels = letters[:n_ds]
    else:
        lab_splits = args.labels.split(',')
        if(len(lab_splits) != n_ds):
            print('Make sure you have the same number of labels')
            print('specified, as you do datasets')
            exit(2)
        else:
            labels=[]
            for i in range(n_ds):
                labels.append(lab_splits[i].strip())
        
    
    # Get rid of that pesky grid dimension

    if(args.lat is None):
        grdid = 0

        for i in range(n_ds):
            hists_ds[i] = hists_ds[i].sel(lndgrid=0)
        
    else:
        print("havent added grid search yet")
        exit(0)
        #for var in hist_ds.variables:
        #if('lndgrid' in var.dims):
        #    for id, dimname in enumerate(var.dims):
        #        if(dimname=='lndgrid'):
        #            lndid = id
                


    # Prep the time-series arrays, and check to see if they
    # are split (ie for a restart partially through at new time)

    tlists = []

    for i in range(n_ds):

       
        #dyear1, split_time1, split_id1 = PrepNCTime(hists_ds[i])
        #tlists.append(dyear1)
        #hists_ds[i]['time'] = hists_ds[i].time_bounds.mean(dim="hist_interval")
        
        hists_ds[i]['delta_sec'] = (hists_ds[i].time_bounds.isel(hist_interval=1) - \
                                    hists_ds[i].time_bounds.isel(hist_interval=0))/ \
                                    np.timedelta64(1, 's')
        
        #code.interact(local=dict(globals(), **locals()))
        #datetimeindex = hists_ds[i].indexes['time'].to_datetimeindex()
        #hists_ds[i]['time'] = datetimeindex
        
        #if(i==0):
        #    split_time=split_time1
        #    split_id=split_id1
        #else:
        #    if(split_id != split_id1):
        #        print("If you have more than 1 dataset,")
        #        print("and 1 of them is discontinuous (split), then they must")
        #        print("both be, and have the same time signatures")
        #        print(split_id)
        #        print(split_id1)
        #        exit(2)


    
                
            
    for i in range(n_ds):
        hists_ds[i] = CreateNewDataArrays(hists_ds[i])
        hist_ds = hists_ds[i]
        #code.interact(local=dict(globals(), **locals()))

    tlists = []
    for i in range(n_ds):
        tlists.append(hists_ds[i]['time'])

    if(True):
        for i in range(n_ds):
            PlotPFTSeries(tlists[i],hists_ds[i])
           
    if(False):
        levscls = np.insert(hists_ds[0].fates_levscls.data,0 , 0.0, axis=0)
        vlists = []
        for i in range(n_ds):
            vlists.append(Prep2DPlantStocks(hists_ds[i]))
            Plot2DTSeries(tlists[i],vlists[i],levscls,dbhcoord)

    if(False):
        levscls = np.insert(hists_ds[0].fates_levscls.data,0 , 0.0, axis=0)
        vlists = []
        for i in range(n_ds):
            vlists.append(PrepNEfficiency2D(hists_ds[i]))
            Plot2DTSeries(tlists[i],vlists[i],levscls,dbhcoord)

    
        
    # Evaluate key carbon states in both FATES, ELM and both together
    vlists = []
    for i in range(n_ds):
        vlists.append(PrepCStateOut(hists_ds[i]))
    PlotTSeries(tlists,vlists,labels)

    # Evaluate canopy versus understory variables
    #vlists = []
    #for i in range(n_ds):
    #    vlists.append(PrepCanUnder(hists_ds[i]))
    #PlotTSeries(tlists,vlists,labels)

    
    # Evaluate discrepancies between stocks and integrated fluxes
    if(False):
        vlists = []
        for i in range(n_ds):
            vlists.append(PrepDiscrepancies(hists_ds[i]))
        PlotTSeries(tlists,vlists,labels)
    
    if(True):
        vlists = []
        for i in range(n_ds):
            vlists.append(PrepResp(hists_ds[i]))
        PlotTSeries(tlists,vlists,labels)
    
    # Evaluate 2D soil variables

    if(False):
        levdcmp = np.insert(-hists_ds[0].levdcmp.data,0 , 0.0, axis=0)
        vlists = []
        for i in range(n_ds):
            vlists.append(Prep2DSoilLittCStocks(hists_ds[i]))
            Plot2DTSeries(tlists[i],vlists[i],levdcmp,depthcoord)

    
    # Check to see if we have some nutrient variables active (this is flawed)
    use_np = False
    for i in range(n_ds):
        if(("FPI_N" in hists_ds[0]) and ("FPI_P" in hists_ds[0])):
            use_np = True
            if(i>0):
                if(not(("FPI_N" in hists_ds[i]) and ("FPI_P" in hists_ds[i]))):
                    print("Inconsitent existence of variables in datasets")
                    exit(2)

    #use_np=False
    if(use_np):


        #vlists = []
        #for i in range(n_ds):
        #    vlists.append(PrepIntNRates(hists_ds[i]))
        #PlotTSeries(tlists,vlists,labels)

        vlists = []
        for i in range(n_ds):
            vlists.append(PrepMeanNRates(hists_ds[i]))
        PlotTSeries(tlists,vlists,labels)

        vlists = []
        for i in range(n_ds):
            vlists.append(PrepMeanPRates(hists_ds[i]))
        PlotTSeries(tlists,vlists,labels)
        
        vlists = []
        for i in range(n_ds):
            vlists.append(PrepMineralNRates(hists_ds[i]))
        PlotTSeries(tlists,vlists,labels)


        vlists = []
        for i in range(n_ds):
            vlists.append(PrepNPStates(hists_ds[i]))
        PlotTSeries(tlists,vlists,labels)
        
        #        vlists = []
        #for i in range(n_ds):
        #    vlists.append(PrepMineralNRatios(hists_ds[i]))
        #StackedAreaPlots(vlists,labels)


    

    #BA_SCLS(time, fates_levscls, lndgrid) ;
    
    code.interact(local=dict(globals(), **locals()))
        
# =======================================================================================
# This routine calculates variables that are products of other raw variables found
# in the dataset
# =======================================================================================
        
def CreateNewDataArrays(hist_dsy):


    # Note about units:
    #
    # All FATES output is kg/m2/s
    # ELM units appear to be mostly g/m2/s

    # NBP [kgC/m2/s]
    hist_dsy["NBP"] = (hist_dsy["FATES_NEP"] + \
                      hist_dsy["FATES_SEEDS_IN_EXTERN_EL"].isel(fates_levelem=0) - \
                    # hist_dsy["FATES_FIRE_CLOSS"])*hist_dsy["FATES_FRACTION"] - \
                      hist_dsy["FATES_FIRE_CLOSS"]) - \
                      hist_dsy["SOM_C_LEACHED"]*kg_per_g

    # CO2 in PPM
    hist_dsy["CO2_PPM"] = 1.e6*hist_dsy["PCO2"]/hist_dsy["PBOT"]

    hist_dsy["UNIT_NUPTAKE"] = (hist_dsy["FATES_NH4UPTAKE"]+hist_dsy["FATES_NO3UPTAKE"]) / hist_dsy["FATES_FROOTC"]
    
    
    fines = hist_dsy["FATES_LITTER_AG_FINE_EL"]+hist_dsy["FATES_LITTER_BG_FINE_EL"]

    # Total un-fragmented litter kgC/m2
    hist_dsy["FATES_TOTLITC"] = hist_dsy["FATES_LITTER_AG_FINE_EL"].isel(fates_levelem=0) + \
                                hist_dsy["FATES_LITTER_BG_FINE_EL"].isel(fates_levelem=0) + \
                                hist_dsy["FATES_LITTER_AG_CWD_EL"].isel(fates_levelem=0) + \
                                hist_dsy["FATES_LITTER_BG_CWD_EL"].isel(fates_levelem=0)
                                  
    # Total litter (fragmented and unfragmented in HLM) [kg/m2]
    hist_dsy["TOTLITC_BOTH"] = hist_dsy["TOTLITC"]*kg_per_g + hist_dsy["FATES_TOTLITC"]

    #code.interact(local=dict(globals(), **locals()))

        
    # ELM Net C Flux Rate [kgC/m2/s]
    hist_dsy["NETCFLUX_ELM"] = hist_dsy["FATES_LITTER_OUT_EL"].isel(fates_levelem=0) - \
                               hist_dsy["HR"]*kg_per_g - \
                               hist_dsy["SOM_C_LEACHED"]*kg_per_g

    # Total ELM Carbon [kgC/m2]
    hist_dsy["TOTC_ELM"] = hist_dsy["TOTSOMC"]*kg_per_g + \
                           hist_dsy["TOTLITC"]*kg_per_g #+ \
#                           hist_dsy["COL_CTRUNC"]*kg_per_g


    # Net Allocated Carbon (Ie NPP-efflux) [kgC/m2/s]
    hist_dsy["FATES_NAC"] = hist_dsy["FATES_NPP"] -  hist_dsy["FATES_EXCESS_RESP"]

    
    
    # kgC/m2 Total FATES Carbon
    hist_dsy["TOTC_FATES"] = hist_dsy["FATES_VEGC"] + hist_dsy["FATES_SEED_BANK"] + hist_dsy["FATES_TOTLITC"]


    # Total combined carbon [kgC/m2]
    hist_dsy["TOTC_BOTH"] = hist_dsy["TOTC_ELM"] + hist_dsy["TOTC_FATES"]

    # FATES Side Carbon Balance
    # Net FATES Production [kgC/m2/s]

    hist_dsy["NFP"] = hist_dsy["FATES_NPP"] + \
                  hist_dsy["FATES_SEEDS_IN_EXTERN_EL"].isel(fates_levelem=0) - \
                  hist_dsy["FATES_FIRE_CLOSS"] - \
                  hist_dsy["FATES_LITTER_OUT_EL"].isel(fates_levelem=0)

    # Fates C balance error [kgC/m2]
    hist_dsy["FCE"] = np.cumsum(hist_dsy["NFP"]*hist_dsy['delta_sec']) - (hist_dsy["TOTC_FATES"] - hist_dsy["TOTC_FATES"].isel(time=0))

    hist_dsy["FCE1"] = np.cumsum(hist_dsy["FATES_NPP"]*hist_dsy['delta_sec'])
    hist_dsy["FCE2"] = np.cumsum(hist_dsy["FATES_SEEDS_IN_EXTERN_EL"].isel(fates_levelem=0)*hist_dsy['delta_sec'])
    hist_dsy["FCE3"] = np.cumsum(hist_dsy["FATES_FIRE_CLOSS"]*hist_dsy['delta_sec'])
    hist_dsy["FCE4"] = np.cumsum(hist_dsy["FATES_LITTER_OUT"]*hist_dsy['delta_sec'])

    #FATES_GPP_PF
    #FATES_VEGC_PF
    #FATES_L2FR_CLSZPF
    #FATES_NPP_SZPF
    #FATES_DDBH_CANOPY_SZPF
    #FATES_DDBH_USTORY_SZPF
    #FATES_CROWNAREA_PF

    #FATES_NPP_PF/FATES_NPLANT_PF
    #FATES_VEGC_PF
    #FATES_MORTALITY_PF, FATES_NPLANT_PF
    #FATES_DDBH_CANOPY_SZPF
    #FATES_DDBH_USTORY_SZPF
   

    #FATES_RECL2FR_CANOPY_PF
    #FATES_RECL2FR_USTORY_PF
    
    
    # Total carbon error [kgC/m2]
    #hist_dsy["TCE"] = np.cumsum(hist_dsy["NBP"]*hist_dsy['delta_sec'])-(hist_dsy["TOTC_BOTH"] - hist_dsy["TOTC_BOTH"].isel(time=0))

    # ELM C balance error [kgC/m2]
    hist_dsy["ECE"] = np.cumsum(hist_dsy["NETCFLUX_ELM"]*hist_dsy['delta_sec']) - (hist_dsy["TOTC_ELM"] - hist_dsy["TOTC_ELM"].isel(time=0))

    use_n = False
    if(("ACTUAL_IMMOB" in hist_dsy) and ("POTENTIAL_IMMOB" in hist_dsy)):
        hist_dsy["FPI_N"] = hist_dsy["ACTUAL_IMMOB"]/hist_dsy["POTENTIAL_IMMOB"]
        use_n = True

    use_p = False
    if(("ACTUAL_IMMOB_P" in hist_dsy) and ("POTENTIAL_IMMOB_P" in hist_dsy)):
        hist_dsy["FPI_P"] = hist_dsy["ACTUAL_IMMOB_P"]/hist_dsy["POTENTIAL_IMMOB_P"]
        use_p = True

    # These are all in units of grams, m2 and s
    # =================================================
    if(use_n):
    
        hist_dsy["FATES_TOTLITN"] = (hist_dsy["FATES_LITTER_AG_FINE_EL"].isel(fates_levelem=1) + \
                                     hist_dsy["FATES_LITTER_BG_FINE_EL"].isel(fates_levelem=1) + \
                                     hist_dsy["FATES_LITTER_AG_CWD_EL"].isel(fates_levelem=1) + \
                                     hist_dsy["FATES_LITTER_BG_CWD_EL"].isel(fates_levelem=1))

        hist_dsy["TOTLITN_BOTH"] = hist_dsy["TOTLITN"] + g_per_kg*hist_dsy["FATES_TOTLITN"]
        
        hist_dsy["ALL_NMIN"] = hist_dsy["GROSS_NMIN"] + hist_dsy["NDEP_TO_SMINN"] + hist_dsy["NFIX_TO_SMINN"]

        hist_dsy["N_LITT_FLUX"] = hist_dsy["FATES_LITTER_OUT_EL"].isel(fates_levelem=1) #- hist_dsy["FATES_NEFFLUX"]

        hist_dsy["FATES_NUPTAKE"] = hist_dsy["FATES_NH4UPTAKE"]+hist_dsy["FATES_NO3UPTAKE"]
        
    if(use_p):

        hist_dsy["FATES_TOTLITP"] =(hist_dsy["FATES_LITTER_AG_FINE_EL"].isel(fates_levelem=2) + \
                                    hist_dsy["FATES_LITTER_BG_FINE_EL"].isel(fates_levelem=2) + \
                                    hist_dsy["FATES_LITTER_AG_CWD_EL"].isel(fates_levelem=2) + \
                                    hist_dsy["FATES_LITTER_BG_CWD_EL"].isel(fates_levelem=2))
        
        hist_dsy["TOTLITP_BOTH"] = g_per_kg*hist_dsy["FATES_TOTLITP"] + hist_dsy["TOTLITP"]

        # kg/m2/s
        hist_dsy["P_LITT_FLUX"] = hist_dsy["FATES_LITTER_OUT_EL"].isel(fates_levelem=2) #-hist_dsy["FATES_PEFFLUX"]
        
    
    return hist_dsy


def PrepCanUnder(hist_dsy):

    varlist = []
    hist_dsy["L2FR_CANOPY"] = (hist_dsy["FATES_L2FR_CANOPY_SZPF"].sum(dim='fates_levpft') * \
                               hist_dsy["FATES_FROOTMAINTAR_CANOPY_SZ"]).sum(dim='fates_levscls') \
                               / hist_dsy["FATES_FROOTMAINTAR_CANOPY_SZ"].sum(dim='fates_levscls')

    hist_dsy["L2FR_USTORY"] = (hist_dsy["FATES_L2FR_USTORY_SZPF"].sum(dim='fates_levpft') * \
                               hist_dsy["FATES_FROOTMAINTAR_USTORY_SZ"]).sum(dim='fates_levscls') \
                               / hist_dsy["FATES_FROOTMAINTAR_USTORY_SZ"].sum(dim='fates_levscls')
    
    return varlist

# =======================================================================================
# This routine preps output variables for the carbon state
# =======================================================================================  

def PrepCStateOut(hist_dsy):
        
    varlist = []

   

    
    # THE NPP IN THE HISTORY DOES NOT CONTAIN EXCESS RESPIRATION IN IT
    hist_dsy["FATES_TOT_NPP"] = hist_dsy["FATES_NPP"] - hist_dsy["FATES_EXCESS_RESP"]    # 'kg m-2 s-1' -  kg m-2 s-1
    hist_dsy["FROOTC_TURN_EFF"] = (hist_dsy["FATES_FROOTCTURN_USTORY_SZ"].sum(dim='fates_levscls') + \
                                   hist_dsy["FATES_FROOTCTURN_CANOPY_SZ"].sum(dim='fates_levscls')) / \
                                   hist_dsy["FATES_TOT_NPP"]
    
    
    #varlist.append(varlist_type('HR','[kgC/m2/yr]',hist_dsy["HR"].rolling(time=roll_w_size).mean()*sec_per_year*kg_per_g))
    #varlist.append(varlist_type('Frag. Flux C','[kgC/m2/yr]',hist_dsy["FATES_LITTER_OUT_EL"].isel(fates_levelem=0).rolling(time=roll_w_size).mean()*sec_per_year))
    #varlist.append(varlist_type('log(abs(NBP))','[kgC/m2/yr]',np.log(np.abs(hist_dsy["NBP"].rolling(time=roll_w_size).mean()*sec_per_year))))
    #CO2_PPM


     
    #varlist.append(varlist_type('Free N Fixation','[gN/m2/yr]',hist_dsy["NFIX_TO_SMINN"].rolling(time=roll_w_size).mean()*sec_per_year))
    varlist.append(varlist_type('CO2','[ppm]',hist_dsy["CO2_PPM"].rolling(time=roll_w_size).mean()))

    varlist.append(varlist_type('GPP','[kgC/m2/yr]',hist_dsy["FATES_GPP"].rolling(time=roll_w_size).mean()*sec_per_year))
    #varlist.append(varlist_type('TRIM','[-]',hist_dsy["FATES_TRIMMING"].rolling(time=roll_w_size).mean()))
    #varlist.append(varlist_type('Sym. N Fixation','[gN/m2/yr]',hist_dsy["FATES_SYMNFIX"].rolling(time=roll_w_size).mean()*1000.0*sec_per_year))
    
    #code.interact(local=dict(globals(), **locals()))
    if("FATES_STOREC_TF" in hist_dsy.data_vars):
        hist_dsy["STORECN_TF"] = hist_dsy["FATES_STOREC_TF"]/hist_dsy["FATES_STOREN_TF"]
    else:
        hist_dsy["STORECN_TF"] = hist_dsy["FATES_STOREC_TFRAC"]/hist_dsy["FATES_STOREN_TFRAC"]
    
    
    hist_dsy["CUE"] = hist_dsy["FATES_TOT_NPP"]/hist_dsy["FATES_GPP"]
    da_masked = xr.ufuncs.isinf(hist_dsy["CUE"])
    hist_dsy["CUE"].data[da_masked] = np.nan
    #code.interact(local=dict(globals(), **locals()))
    hist_dsy["FATES_NUP_EFF"] = (hist_dsy["FATES_NH4UPTAKE"]+hist_dsy["FATES_NH4UPTAKE"])*sec_per_year/hist_dsy["FATES_FROOTC"]
    
    #varlist.append(varlist_type('N Uptake Efficiency','[kgN/kgC/yr]',hist_dsy["FATES_NUP_EFF"].rolling(time=roll_w_size).mean()))
    
    hist_dsy["FATES_CALLOC_EFF"] = hist_dsy["FATES_EXCESS_RESP"] / hist_dsy["FATES_TOT_NPP"]
    hist_dsy["FATES_NALLOC_EFF"] = hist_dsy["FATES_NEFFLUX"] / (hist_dsy["FATES_NH4UPTAKE"]+hist_dsy["FATES_NO3UPTAKE"])
    hist_dsy["FATES_PALLOC_EFF"] = hist_dsy["FATES_PEFFLUX"] / hist_dsy["FATES_PUPTAKE"]

    #varlist.append(varlist_type('Overflow Respiration','[kgC m-2 yr-1]',hist_dsy["FATES_EXCESS_RESP"].rolling(time=roll_w_size).mean()*sec_per_year))
    varlist.append(varlist_type('NPP','[kgC/m2/yr]',hist_dsy["FATES_TOT_NPP"].rolling(time=roll_w_size).mean()*sec_per_year))
    varlist.append(varlist_type('Basal Area','[m2/ha]',hist_dsy["FATES_BASALAREA_SZ"].sum(dim='fates_levscls').rolling(time=roll_w_size).mean()*10000.0))
    varlist.append(varlist_type('LAI','[m2/m2]',hist_dsy["TLAI"].rolling(time=roll_w_size).mean()))
    #varlist.append(varlist_type('AGB','[kgC/m2]',hist_dsy["FATES_VEGC_ABOVEGROUND"].rolling(time=roll_w_size).mean()))

    
    varlist.append(varlist_type('CUE','[-]',hist_dsy["CUE"].rolling(time=roll_w_size).mean()))
    varlist.append(varlist_type('$\lambda$ (fine-root fraction)','[-]',hist_dsy["FATES_L2FR"].rolling(time=roll_w_size).mean()))
    ##varlist.append(varlist_type('$f_{CN(N)}$','[log(gC/gC/gN/gN)]',np.log(hist_dsy["STORECN_TF"].rolling(time=roll_w_size).mean()),mincap=-1,maxcap=1))
    ##varlist.append(varlist_type('C Allocation Inefficiency','[/]',hist_dsy["FATES_CALLOC_EFF"].rolling(time=roll_w_size).mean()))
    #varlist.append(varlist_type('N Allocation Inefficiency','[/]',hist_dsy["FATES_NALLOC_EFF"].rolling(time=roll_w_size).mean()))

    varlist.append(varlist_type('Plant Unit N Uptake','[gN/gC/year]',hist_dsy["UNIT_NUPTAKE"].rolling(time=roll_w_size).mean()*sec_per_year  ))

    
    hist_dsy["SMIN_N"] = hist_dsy["SMIN_NO3"]+hist_dsy["SMIN_NH4"]
    varlist.append(varlist_type('Mineralized NH4+NO3','[gN/m2]',hist_dsy["SMIN_N"].rolling(time=roll_w_size).mean()))
    #varlist.append(varlist_type('Mineralized NO3','[gN/m2]',hist_dsy["SMIN_NO3"].rolling(time=roll_w_size).mean()))
    # USING TOTVEGC AS PROXY FOR AGB: NOTE THE AGB FRACTION OF 0.6 !!!
    ##varlist.append(varlist_type('Organic Soil N','[gN/m2]',hist_dsy["TOTSOMN"].rolling(time=roll_w_size).mean()))
    
    
   
    
    
    #area_si_age
    # mean age of forest

    #fates_history_ageclass_bin_edges = 0, 1, 2, 5, 10, 20, 50
    # age_bin_centers = 0.5,1.5,3.5,7.5,15,35,+
    
    #code.interact(local=dict(globals(), **locals()))
    
    #hist_dsy["FATES_PATCHAREA_AP"]
    #hist_dsy["FATES_PATCHAREA_AP"].dims = ('time', 'fates_levage')



    
    #varlist.append(varlist_type('Total Litter C','[kgC/m2]',hist_dsy["TOTLITC_BOTH"].rolling(time=roll_w_size).mean()))
    #varlist.append(varlist_type('Total Litter N','[gN/m2]',hist_dsy["TOTLITN_BOTH"].rolling(time=roll_w_size).mean()))
    #varlist.append(varlist_type('Organic Soil C','[kgC/m2]',hist_dsy["TOTSOMC"].rolling(time=roll_w_size).mean()*kg_per_g))
    #varlist.append(varlist_type('NAC','[kgC/m2/yr]',hist_dsy["FATES_NAC"].rolling(time=roll_w_size).mean()*sec_per_year))
    #varlist.append(varlist_type('Plant C Store Frac.','[gC/gC]',hist_dsy["FATES_STOREC_TF"].rolling(time=roll_w_size).mean()))

    hist_dsy["CO_PER_PA"] = hist_dsy["FATES_NCOHORTS"]/hist_dsy["FATES_NPATCHES"]
    #varlist.append(varlist_type('Cohorts per Patch','[/]',hist_dsy["CO_PER_PA"].rolling(time=roll_w_size).mean()))
    #varlist.append(varlist_type('N patches','[/]',hist_dsy["FATES_NPATCHES"].rolling(time=roll_w_size).mean()))
    
    #code.interact(local=dict(globals(), **locals()))

    #varlist.append(varlist_type('$f_{N(N)}$','[gN/gN]',hist_dsy["FATES_STOREN_TF"].rolling(time=roll_w_size).mean()))
    #varlist.append(varlist_type('$f_{C}$','[gC/gC]',hist_dsy["FATES_STOREC_TF"].rolling(time=roll_w_size).mean()))
    
    #varlist.append(varlist_type('Plant C Store Frac.','[gC/gC]',hist_dsy["FATES_STOREC_TF"].rolling(time=roll_w_size).mean()))
    #varlist.append(varlist_type('Plant N Store Frac.','[gC/gC]',hist_dsy["FATES_STOREN_TF"].rolling(time=roll_w_size).mean()))
    #varlist.append(varlist_type('Plant P Store Frac.','[gC/gC]',hist_dsy["FATES_STOREP_TF"].rolling(time=roll_w_size).mean()))

    

    #hist_dsy["STOREPC_TF"] = hist_dsy["FATES_STOREP_TF"]/hist_dsy["FATES_STOREC_TF"]
    #varlist.append(varlist_type('Plant P/C Store Frac.','[gP/gP/gC/gC]',hist_dsy["STOREPC_TF"].rolling(time=roll_w_size).mean()))
    
    
    

    #varlist.append(varlist_type('$\lambda$recc','[-]',hist_dsy["FATES_RECL2FR_CANOPY_PF"].mean(dim='fates_levpft').rolling(time=roll_w_size).mean()))
    #varlist.append(varlist_type('$\lambda$recu','[-]',hist_dsy["FATES_RECL2FR_USTORY_PF"].mean(dim='fates_levpft').rolling(time=roll_w_size).mean()))
    #varlist.append(varlist_type('Storage Respiration','[kgC m-2 yr-1]',hist_dsy["FATES_EXCESS_RESP"].rolling(time=roll_w_size).mean()*sec_per_year))

    #varlist.append(varlist_type('Total Litter N','[gN/m2]',hist_dsy["TOTLITN_BOTH"].rolling(time=roll_w_size).mean()))
    #varlist.append(varlist_type('Mineralized NH4','[gN/m2]',hist_dsy["SMIN_NH4"].rolling(time=roll_w_size).mean()))


    # Uptake efficiency [kg/kg/year]
    
    
    
    #varlist.append(varlist_type('P Allocation Loss Fraction','[/]',hist_dsy["FATES_PALLOC_EFF"].rolling(time=roll_w_size).mean()))

    #varlist.append(varlist_type('Root Turnover Fraction','[/]',hist_dsy['FROOTC_TURN_EFF'].rolling(time=roll_w_size).mean()))
    
    #varlist.append(varlist_type('Frac. Immob. (N)','[]',hist_dsy["FPI_N"].rolling(time=roll_w_size).mean()))
    
    #varlist.append(varlist_type('L2FR EMA','[/]',hist_dsy["FATES_L2FR_EMA"].rolling(time=roll_w_size).mean()))

    #hist_dsy["DELTA_L2FR"] = hist_dsy["FATES_L2FR"] - hist_dsy["FATES_L2FR_EMA"]
     
    #varlist.append(varlist_type('L2FR-L2FR_EMA','[/]',hist_dsy["DELTA_L2FR"]))
    
    return varlist



def PrepNEfficiency2D(hist_dsy):
    # NEED/FNRTC  # gN/gC/day
    # NNEED_SCPF:units = "kgN d-1 ha-1" ;
    # FNRTC:units = "kgC ha-1" ;
    
    #hist_dsy["NNEED_PER_FNRTC_SCPF"] = hist_dsy["NNEED_SCPF"] / hist_dsy["FNRTC_SCPF"]

    #hist_dsy["NUPTAKE_PER_FNRTC_SCPF"] = (hist_dsy["NH4UPTAKE_SCPF"]+hist_dsy["NO3UPTAKE_SCPF"] )/ hist_dsy["FNRTC_SCPF"]
    
    #hist_dsy["NALLOC_PER_FNRTC_SCPF"] = (hist_dsy["NH4UPTAKE_SCPF"]+hist_dsy["NO3UPTAKE_SCPF"]-hist_dsy["NEFFLUX_SCPF"]) / hist_dsy["FNRTC_SCPF"]
    
    #hist_dsy["FNRTC_PLANT"] = hist_dsy["FNRTC_SCPF"] / hist_dsy["NPLANT_SCPF"]
    #hist_dsy["NNEED_PLANT"] = 1000.*hist_dsy["NNEED_SCPF"] / hist_dsy["NPLANT_SCPF"]
    #hist_dsy["NUPTAKE_PLANT"] = 1000.*(hist_dsy["NH4UPTAKE_SCPF"]+hist_dsy["NO3UPTAKE_SCPF"] ) / hist_dsy["NPLANT_SCPF"]
    
    #hist_dsy["NNUPTAKENNEED"] = (hist_dsy["NH4UPTAKE_SCPF"]+hist_dsy["NO3UPTAKE_SCPF"] )/hist_dsy["NNEED_SCPF"]

    

    
    # Storage as a function of FRNRTC

    
    

    
    
    varlist = []
#    varlist.append(varlist_type('FNRTC','[kgC]',hist_dsy["FNRTC_PLANT"]))
#    varlist.append(varlist_type('N NEED','[kgN/day]',hist_dsy["NNEED_PLANT"]))
#    varlist.append(varlist_type('N UPTAKE','[gN/day]',hist_dsy["NUPTAKE_PLANT"],use_log=True))
    
#    varlist.append(varlist_type('N UPTAKE / FNRTC','[gN/gC/day]',hist_dsy["NUPTAKE_PER_FNRTC_SCPF"]))

    #varlist.append(varlist_type('N UPTAKE / N NEED','[gN/gN/day]',hist_dsy["NNUPTAKENNEED"],use_log=True,mincap=0.001,maxcap=1000.0))
    #varlist.append(varlist_type('N Need / FNRTC','[gN/gC/day]',hist_dsy["NUPTAKE_PER_FNRTC_SCPF"]))
    #varlist.append(varlist_type('N Store Frac.','[gN/gN]',hist_dsy["STOREN_TF_SCPF"]))
    #varlist.append(varlist_type('L2FR','[-]',hist_dsy["FATES_L2FR_SZPF"]))
    
    #varlist.append(varlist_type('N USAGE / FNRTC','[gN/gC/day]',hist_dsy["NALLOC_PER_FNRTC_SCPF"]))

    return varlist

    
# =======================================================================================
# This routine prepares output variables for the discrepancy between stocks
# and integrated fluxes
# =======================================================================================        

def PrepDiscrepancies(hist_dsy):

    varlist = []
    varlist.append(varlist_type('Integrated Fates \n Carbon Discrepancy','[kgC/m2]',hist_dsy["FCE"]))
    varlist.append(varlist_type('Integrated ELM \n Carbon Discrepancy','[kgC/m2]',hist_dsy["ECE"]))
    varlist.append(varlist_type('Total Carbon \n Discrepancy','[kgC/m2]',hist_dsy["TCE"].rolling(time=roll_w_size).mean()))
    #varlist.append(varlist_type('FCE1','[kgC/m2]',hist_dsy["FCE1"]))
    #varlist.append(varlist_type('FCE2','[kgC/m2]',hist_dsy["FCE2"]))
    #varlist.append(varlist_type('FCE3','[kgC/m2]',hist_dsy["FCE3"]))
    #varlist.append(varlist_type('FCE4','[kgC/m2]',hist_dsy["FCE4"]))
    
    #varlist.append(varlist_type('Total C Bal. Disc.','[kgC/m2]',hist_dsy["TCE"].rolling(time=roll_w_size).mean()))

    return(varlist)


def PrepResp(hist_dsy):
        
        
    varlist = []
    # Fine-Root Respiration Absolute (maintenance + growth)
    # FATES_GROWTH_RESP', units='kg m-2 s-1
    # FATES_LEAFMAINTAR
    # FATES_FROOTMAINTAR
    # FATES_CROOTMAINTAR
    # FATES_LSTEMMAINTAR

    # 'kg m-2 s-1',
    hist_dsy["FNRT_COST"] =  (hist_dsy["FATES_FROOTMAINTAR"]+hist_dsy["FATES_FROOTCTURN_USTORY_SZ"].sum(dim='fates_levscls') + \
        hist_dsy["FATES_FROOTCTURN_CANOPY_SZ"].sum(dim='fates_levscls'))/hist_dsy["FATES_NPP"]

    
    hist_dsy["MR_TOTAL"] = hist_dsy["FATES_FROOTMAINTAR"]+hist_dsy["FATES_LSTEMMAINTAR"]+hist_dsy["FATES_CROOTMAINTAR"]+hist_dsy["FATES_LEAFMAINTAR"]
    hist_dsy["FNRT_MR_FRAC"] = hist_dsy["FATES_FROOTMAINTAR"] / hist_dsy["MR_TOTAL"]
    hist_dsy["LEAF_MR_FRAC"] = hist_dsy["FATES_LEAFMAINTAR"] / hist_dsy["MR_TOTAL"]
    hist_dsy["LSTEM_MR_FRAC"] = hist_dsy["FATES_LSTEMMAINTAR"] / hist_dsy["MR_TOTAL"]
    hist_dsy["CROOT_MR_FRAC"] = hist_dsy["FATES_CROOTMAINTAR"] / hist_dsy["MR_TOTAL"]
    
    varlist.append(varlist_type('FNRT MR','[kgC/m2/yr]',hist_dsy["FATES_FROOTMAINTAR"].rolling(time=roll_w_size).mean()*sec_per_year ))
    varlist.append(varlist_type('FNRT MR Fraction','[/]',hist_dsy["FNRT_MR_FRAC"].rolling(time=roll_w_size).mean()))
    varlist.append(varlist_type('LEAF MR Fraction','[/]',hist_dsy["LEAF_MR_FRAC"].rolling(time=roll_w_size).mean()))
    varlist.append(varlist_type('STEM MR Fraction','[/]',hist_dsy["LSTEM_MR_FRAC"].rolling(time=roll_w_size).mean()))
    varlist.append(varlist_type('CSRT MR Fraction','[/]',hist_dsy["CROOT_MR_FRAC"].rolling(time=roll_w_size).mean()))
    
    hist_dsy["FNRT_GR_FRAC"] = hist_dsy["FATES_FROOT_ALLOC"] / (hist_dsy["FATES_FROOT_ALLOC"] + hist_dsy["FATES_LEAF_ALLOC"] + hist_dsy["FATES_SEED_ALLOC"] + hist_dsy["FATES_STEM_ALLOC"] + hist_dsy["FATES_CROOT_ALLOC"] + hist_dsy["FATES_STORE_ALLOC"])

    hist_dsy["FNRT_GR"] = hist_dsy["FATES_GROWTH_RESP"] * hist_dsy["FNRT_GR_FRAC"]

    varlist.append(varlist_type('FNRT COST FRAC','[/]',hist_dsy["FNRT_COST"].rolling(time=roll_w_size).mean() ))
    varlist.append(varlist_type('FNRT GR Fraction','[/]',hist_dsy["FNRT_GR_FRAC"].rolling(time=roll_w_size).mean()))
    varlist.append(varlist_type('GPP','[kgC/m2/yr]',hist_dsy["FATES_GPP"].rolling(time=roll_w_size).mean()*sec_per_year))
    
    return varlist
    
# =======================================================================================
# This routine prepares 2D output variables for the discrepancy between stocks
# and integrated fluxes
# =======================================================================================

def Prep2DSoilLittCStocks(hist_dsy):
        
        
    varlist = []
    varlist.append(varlist_type('Met Soil C','[gc/m3]',hist_dsy["SOIL1C_vr"]))
    varlist.append(varlist_type('Cel Soil C','[gc/m3]',hist_dsy["SOIL2C_vr"]))
    varlist.append(varlist_type('Lig Soil C','[gc/m3]',hist_dsy["SOIL3C_vr"]))
    
    varlist.append(varlist_type('LITTER 1 C','[gc/m3]',hist_dsy["LITR1C_vr"]))
    varlist.append(varlist_type('LITTER 2 C','[gc/m3]',hist_dsy["LITR2C_vr"]))
    varlist.append(varlist_type('LITTER 3 C','[gc/m3]',hist_dsy["LITR3C_vr"]))
    
    return varlist

def Prep2DPlantStocks(hist_dsy):
        
    
    varlist = []
    varlist.append(varlist_type('Basal Area by dbh','[m2/ha]',hist_dsy["FATES_BASALAREA_SZ"]*10000.0))
    varlist.append(varlist_type('AGB by dbh','[MgC/ha]',hist_dsy["FATES_VEGC_ABOVEGROUND_SZ"]*10.0))

    #[pft,sclass,layer] = UnwrapCLSZPF(hist_dsy.dims['fates_levclscpf'], \
    #             hist_dsy.dims['fates_levpft'], \
    #             hist_dsy.dims['fates_levscls'], \
    #             hist_dsy.dims['fates_levcan'])
    #

    da4 = UnstackCLSZPF(hist_dsy["FATES_L2FR_CLSZPF"], \
                        hist_dsy.dims['fates_levpft'], \
                        hist_dsy.dims['fates_levscls'], \
                        hist_dsy.dims['fates_levcan'])

    hist_dsy["FATES_L2FR_CANOPY_SZPF"] = da4.sel(clscpf_clmap=1).sel(clscpf_ftmap=1)
    hist_dsy["FATES_L2FR_USTORY_SZPF"] = da4.sel(clscpf_clmap=2).sel(clscpf_ftmap=1)
    
    varlist.append(varlist_type(r'$\lambda$ Canopy','[-]',hist_dsy["FATES_L2FR_CANOPY_SZPF"] ))
    varlist.append(varlist_type(r'$\lambda$ Understory','[-]',hist_dsy["FATES_L2FR_USTORY_SZPF"] ))
    #code.interact(local=dict(globals(), **locals()))
    return varlist


def UnstackCLSZPF(da_s,npft,nsclass,nlayer):

    ntime = len(da_s.time)
    nclscpf = len(da_s.fates_levclscpf)
        
    da_u = xr.DataArray( \
        np.zeros((ntime,npft,nsclass,nlayer)), \
        coords=[("time", da_s.time), \
                ("clscpf_ftmap", np.array(range(1,npft+1))),\
                ("clscpf_scmap", np.array(range(1,nsclass+1))), \
                ("clscpf_clmap", np.array(range(1,nlayer+1)))],  )

    for ft in range(npft):
        for cl in range(nlayer):
            for sc in range(nsclass):
                i = ft*nsclass*nlayer + sc*nlayer + cl
                da_u[:,ft,sc,cl] = da_s[:,i]
    return(da_u)


        
def PrepMeanPRates(hist_dsy):

    varlist = []
    varlist.append(varlist_type('P Deposition ','[gP/m2/yr]',hist_dsy["PDEP_TO_SMINP"].rolling(time=roll_w_size).mean()*sec_per_year))
    varlist.append(varlist_type('Plant P Efflux ','[gP/m2/yr]',hist_dsy["FATES_PEFFLUX"].rolling(time=roll_w_size).mean()*g_per_kg*sec_per_year))
    varlist.append(varlist_type('P Litterfall','[gP/m2/yr]',hist_dsy["P_LITT_FLUX"].rolling(time=roll_w_size).mean()*g_per_kg*sec_per_year))
    varlist.append(varlist_type('Plant P Uptake','[gP/m2/yr]',hist_dsy["FATES_PUPTAKE"].rolling(time=roll_w_size).mean()*g_per_kg*sec_per_year))
    varlist.append(varlist_type('P Supplementation','[gP/m2/yr]',hist_dsy['SUPPLEMENT_TO_SMINP'].rolling(time=roll_w_size).mean()*sec_per_year))
    varlist.append(varlist_type('SolP/DT','[gP/m2/yr]',hist_dsy["SOLUTIONP"].rolling(time=roll_w_size).mean()/1800.0*sec_per_year))
    return varlist

def PrepMeanNRates(hist_dsy):


    # Mineralized N Balance
    varlist = []
    

    hist_dsy["NET_N_GAINS"] = hist_dsy["NDEP_TO_SMINN"]*sec_per_year + \
                              hist_dsy["NFIX_TO_SMINN"]*sec_per_year + \
                              hist_dsy["FATES_SEEDS_IN_EXTERN_EL"].isel(fates_levelem=1)*g_per_kg*sec_per_year - \
                              (hist_dsy["SMIN_NO3_LEACHED"]+hist_dsy["SMIN_NO3_RUNOFF"])*sec_per_year - \
                              (hist_dsy["F_N2O_NIT"]+hist_dsy["F_DENIT"])*sec_per_year

    # Kg/m2/s
    #hist_dsy["CANFRAC"] = hist_dsy["PFTcanopycrownarea"].isel(fates_levpft=0)/hist_dsy["PFTcrownarea"].isel(fates_levpft=0)

    #
    
    
    hist_dsy["UNIT_NFIX"] = hist_dsy["NFIX_TO_SMINN"]/(hist_dsy["FATES_NPP"]*1000)
    
    #varlist.append(varlist_type('N Deposition ','[gN/m2/yr]',hist_dsy["NDEP_TO_SMINN"].rolling(time=roll_w_size).mean()*sec_per_year))
    varlist.append(varlist_type('N Fixation','[gN/m2/yr]',hist_dsy["NFIX_TO_SMINN"].rolling(time=roll_w_size).mean()*sec_per_year))
    #varlist.append(varlist_type('N Net Seed Influx','[gN/m2/yr]',hist_dsy["SEEDS_IN_EXTERN_EL"].isel(fates_levelem=1).rolling(time=roll_w_size).mean()*ha_per_m2*day_per_year*g_per_kg))

    varlist.append(varlist_type('Unit N Fixation','[gN/gC]',hist_dsy["UNIT_NFIX"].rolling(time=roll_w_size).mean()))
    
    varlist.append(varlist_type('Plant N Efflux ','[gN/m2/yr]',hist_dsy["FATES_NEFFLUX"].rolling(time=roll_w_size).mean()*g_per_kg*sec_per_year))
    varlist.append(varlist_type('N Litterfall','[gN/m2/yr]',hist_dsy["N_LITT_FLUX"].rolling(time=roll_w_size).mean()*g_per_kg*sec_per_year))
    varlist.append(varlist_type('Plant NH4 Uptake','[gN/m2/yr]',hist_dsy["FATES_NH4UPTAKE"].rolling(time=roll_w_size).mean()*g_per_kg*sec_per_year  ))

    #hist_dsy["UNIT_NUPTAKE"] = (hist_dsy["FATES_NH4UPTAKE"]+hist_dsy["FATES_NO3UPTAKE"]) / hist_dsy["FATES_FROOTC"]
    varlist.append(varlist_type('Plant UNIT UPTAKE','[gN/s/gC]',hist_dsy["UNIT_NUPTAKE"].rolling(time=roll_w_size).mean()))
    #varlist.append(varlist_type('Plant NO3 Uptake','[gN/m2/yr]',hist_dsy["FATES_NO3UPTAKE"].rolling(time=roll_w_size).mean()*g_per_kg*sec_per_year  ))
    
    varlist.append(varlist_type('NO3 Leach+RO Loss','[gN/m2/yr]',(hist_dsy["SMIN_NO3_LEACHED"]+hist_dsy["SMIN_NO3_RUNOFF"]).rolling(time=roll_w_size).mean()*sec_per_year))
    varlist.append(varlist_type('N2O+N2 Losses','[gN/m2/yr]',(hist_dsy["F_N2O_NIT"]+hist_dsy["F_DENIT"]).rolling(time=roll_w_size).mean()*sec_per_year))
    varlist.append(varlist_type('Nitrification Rate','[gN/m2/yr]',hist_dsy["F_NIT"].rolling(time=roll_w_size).mean()*sec_per_year))
    varlist.append(varlist_type('N Supplementation','[gN/m2/yr]',hist_dsy['SUPPLEMENT_TO_SMINN'].rolling(time=roll_w_size).mean()*sec_per_year))
    varlist.append(varlist_type('Gross N Mineralization ','[gN/m2/yr]',hist_dsy["GROSS_NMIN"].rolling(time=roll_w_size).mean()*sec_per_year))
    varlist.append(varlist_type('N Immobilization','[gN/m2/yr]',hist_dsy["ACTUAL_IMMOB"].rolling(time=roll_w_size).mean()*sec_per_year))
    
    
    
    return varlist


def PrepMineralNRates(hist_dsy):

    varlist = []
    varlist.append(varlist_type('Frac. Immob. (N)','[]',hist_dsy["FPI_N"].rolling(time=roll_w_size).mean()))
    #varlist.append(varlist_type('Frac. Immob. (P)','[]',hist_dsy["FPI_P"].rolling(time=roll_w_size).mean()))
    varlist.append(varlist_type('Immob/ Min. Gain ','[-]',hist_dsy["ACTUAL_IMMOB"].rolling(time=roll_w_size).mean()/hist_dsy["ALL_NMIN"].rolling(time=roll_w_size).mean()))
    varlist.append(varlist_type('(N20+N2)/ Min. Gain','[-]',(hist_dsy["F_N2O_NIT"].rolling(time=roll_w_size).mean()+hist_dsy["F_DENIT"].rolling(time=roll_w_size).mean())/hist_dsy["ALL_NMIN"].rolling(time=roll_w_size).mean()))
    varlist.append(varlist_type('(Leach+RO)/ Min. Gain','[-]',(hist_dsy["SMIN_NO3_LEACHED"].rolling(time=roll_w_size).mean()+hist_dsy["SMIN_NO3_RUNOFF"].rolling(time=roll_w_size).mean())/hist_dsy["ALL_NMIN"].rolling(time=roll_w_size).mean()))
    varlist.append(varlist_type('Plant N Uptake/ Min. Gain','[-]',hist_dsy["FATES_NUPTAKE"].rolling(time=roll_w_size).mean()*g_per_kg   /hist_dsy["ALL_NMIN"].rolling(time=roll_w_size).mean()))
    varlist.append(varlist_type('Nitrification/ Min. Gain','[-]',hist_dsy["F_NIT"].rolling(time=roll_w_size).mean()/hist_dsy["ALL_NMIN"].rolling(time=roll_w_size).mean()))
    
    return varlist
    
def PrepMineralNRatios(hist_dsy):

    varlist = []
    varlist.append(varlist_type('Immob/ Min. Gain ','[-]',hist_dsy["ACTUAL_IMMOB"].rolling(time=roll_w_size).mean()/hist_dsy["ALL_NMIN"].rolling(time=roll_w_size).mean()))
    varlist.append(varlist_type('(N20+N2)/ Min. Gain','[-]',(hist_dsy["F_N2O_NIT"].rolling(time=roll_w_size).mean()+hist_dsy["F_DENIT"].rolling(time=roll_w_size).mean())/hist_dsy["ALL_NMIN"].rolling(time=roll_w_size).mean()))
    varlist.append(varlist_type('(Leach+RO)/ Min. Gain','[-]',(hist_dsy["SMIN_NO3_LEACHED"].rolling(time=roll_w_size).mean()+hist_dsy["SMIN_NO3_RUNOFF"].rolling(time=roll_w_size).mean())/hist_dsy["ALL_NMIN"].rolling(time=roll_w_size).mean()))
    varlist.append(varlist_type('Plant N Uptake/ Min. Gain','[-]',hist_dsy["FATES_NUPTAKE"].rolling(time=roll_w_size).mean()*g_per_kg/hist_dsy["ALL_NMIN"].rolling(time=roll_w_size).mean()))

    return varlist



def PrepNPStates(hist_dsy):

    # Nutrient States
    varlist = []
    varlist.append(varlist_type('Organic Soil N','[gN/m2]',hist_dsy["TOTSOMN"].rolling(time=roll_w_size).mean()))
    varlist.append(varlist_type('Total Litter N','[gN/m2]',hist_dsy["TOTLITN_BOTH"].rolling(time=roll_w_size).mean()))
    varlist.append(varlist_type('Mineralized NH4','[gN/m2]',hist_dsy["SMIN_NH4"].rolling(time=roll_w_size).mean()))
    varlist.append(varlist_type('Mineralized NO3','[gN/m2]',hist_dsy["SMIN_NO3"].rolling(time=roll_w_size).mean()))

    varlist.append(varlist_type('Organic Soil P','[gP/m2]',hist_dsy["TOTSOMP"].rolling(time=roll_w_size).mean()))
    varlist.append(varlist_type('Total Litter P','[gP/m2]',hist_dsy["TOTLITP_BOTH"].rolling(time=roll_w_size).mean()))
    varlist.append(varlist_type('Solution Min P','[gP/m2]',hist_dsy["SOLUTIONP"].rolling(time=roll_w_size).mean()))

    
    
    varlist.append(varlist_type('Labile P','[gP/m2]',hist_dsy["LABILEP"].rolling(time=roll_w_size).mean()))
    varlist.append(varlist_type('Secondary P','[gP/m2]',hist_dsy["SECONDP"].rolling(time=roll_w_size).mean()))
    varlist.append(varlist_type('Occluded P','[gP/m2]',hist_dsy["OCCLP"].rolling(time=roll_w_size).mean()))
    
    #    varlist.append(varlist_type('Mineralized P','[gP/m2]',hist_dsy["SMINP"].rolling(time=roll_w_size).mean()))
#    varlist.append(varlist_type('Organic Soil P','[gP/m2]',hist_dsy["TOTSOMP"].rolling(time=roll_w_size).mean()))
#    varlist.append(varlist_type('Total Litter P','[gP/m2]',hist_dsy["TOTLITP_BOTH"].rolling(time=roll_w_size).mean()))
#    varlist.append(varlist_type('Plant N Store','[gN/m2]',hist_dsy["STOREN"].rolling(time=roll_w_size).mean()))
#    varlist.append(varlist_type('Plant N Store Frac.','[gN/gN]',hist_dsy["STOREN_TF"].rolling(time=roll_w_size).mean()))
 #   varlist.append(varlist_type('Plant P Store Frac.','[gP/gP]',hist_dsy["STOREP_TF"].rolling(time=roll_w_size).mean()))
    

    
    
    
    return varlist
    

            
                       



class varlist_type:

    def __init__(self, varname, varunits, data, *args, **kwargs):

        self.name = varname
        self.units = varunits
        self.data  = data
        self.use_log = False
        self.maxcap = None
        self.mincap = None
        
        use_log = kwargs.get('use_log', None)
        maxcap = kwargs.get('maxcap', None)
        mincap = kwargs.get('mincap', None)
        
        if(not(use_log is None)):
            self.use_log = use_log

        if(not(maxcap is None)):
            self.maxcap = maxcap
            
        if(not(mincap is None)):
            self.mincap = mincap
            
        
def StackedAreaPlots(vlists,tlists,labels,splitid): 


    # PLot Time-series data from an xarray dataset
    # There should be an axis break at split_id

    matplotlib.rc('axes.formatter', useoffset=False)

    nsets = len(vlists)
    nvars = len(vlists[0])

    # Choose a qualitative colormap

    #'Accent', 'Dark2', 'Set2'
    #code.interact(local=dict(globals(), **locals()))
    colormap = plt.cm.Set3
    
    for id in range(nsets):

        stfig = plt.figure(figsize=(8, 8))
        stackarray = np.empty([nvars,len(tlists[id])])

        for iv in range(nvars):
            stackarray[iv,:] = vlists[id][iv].data[:]

        plt.stackplot(tlists[id],stackarray)

        plt.show()
        #code.interact(local=dict(globals(), **locals()))




        
def Plot2DTSeries(tlist,varlist_obj,levcoord,levtype): 

    # Create a 2D "pcolor" plot using time as the
    # horizontal (x) axis. Note that it is expected
    # that coordinates are provided for the 2D variable
    # of interest. All values must have coordinate info
    # on their boundaries, and thus the coordiantes must
    # be of dimensions +1 larger than the centers. For
    # time, we fudge this, and just remove the last(first?)
    # value, but for the vertical dimension, you must pass
    # in the boundary array.


    #       J  F   M   A   M    J    J    A    S    O    N    D
    dmo  = [0, 31, 28, 31,  30,  31,  31,  30,  31,  30,  31,  30]
    cdmo = [0, 31, 59, 90, 120, 151, 182, 212, 243, 273, 304, 334]


    ntime=len(tlist.data)
    decyears = [tlist.data[j].year+(tlist.data[j].dayofyr/365.0)+tlist.data[j].hour/(365.0*24.0) for j in range(ntime)]


        
    ni = len(levcoord)
    nj = ntime

    X = np.zeros(shape=(ni,nj),dtype=np.float64)
    Y = np.zeros(shape=(ni,nj),dtype=np.float64)
    
    # (X[i+1, j], Y[i+1, j])      (X[i+1, j+1], Y[i+1, j+1])
    #                   +--------+
    #                   | C[i,j] |
    #                   +--------+
    # (X[i, j],   Y[i, j])        (X[i, j+1], Y[i, j+1]),
    #
    # Time (X, j-index) ----> 
    # Vertical (depth) (Y, i-index)

    for i in range(ni):
        for j in range(nj):
            X[i,j] = decyears[j]
            Y[i,j] = levcoord[i]

    nvars = len(varlist_obj)

    if(nvars==1):
        nc = 1
        nr = 1
        fig, axs = plt.subplots(ncols=1,nrows=1,figsize=(7,7))
    elif(nvars==2):
        nc = 2
        nr = 1
        fig, axs = plt.subplots(ncols=2,nrows=1,figsize=(8,4))
    elif(nvars>2):
        nc = 2
        nr = int(np.ceil(float(nvars)/2.0))
        fig, axs = plt.subplots(ncols=2,nrows=nr,figsize=(7,nr*3))
        
    ax1s = axs.reshape(-1)

    for iv in range(nvars):

        ax1 = ax1s[iv]
        
        irow = int(np.floor(float(iv+.01)/float(nc))+1)
        icol = np.mod(iv,nc)+1

        ax1.grid('on')
            
        var = varlist_obj[iv].data[:-1,:].values.transpose()

        if(not(varlist_obj[iv].maxcap == None)):
            vmax = varlist_obj[iv].maxcap
        else:
            vmax = np.nanpercentile(var,97)

        if(not(varlist_obj[iv].mincap == None)):
            vmin = varlist_obj[iv].mincap
        else:
            vmin = np.nanpercentile(var,3)

#        smalloffset = 1e-20
        if(varlist_obj[iv].use_log==True):
            var = np.log(var)#+smalloffset)
#            vmax = np.log(vmax)#+smalloffset)
#            vmin = np.log(vmin)#+smalloffset)

        if(iv>1):
            vmax = 0.7
            vmin = 0.1
        pc1 = ax1.pcolor(X,Y,var,cmap='Oranges',vmin=vmin,vmax=vmax)
        #pc1 = ax1.pcolor(X,Y,var,cmap='Oranges')
        ax1.set_title('{} {}'.format(varlist_obj[iv].name,varlist_obj[iv].units))
            
        if(irow<nr):
            ax1.set_xticklabels('')
        else:
            ax1.set_xlabel('Year')

        if(icol>1):
            ax1.set_yticklabels('')
        else:
            #ax1.set_ylabel('Depth [m]')
            if(levtype==dbhcoord):
                ax1.set_ylabel('DBH [cm]')
            if(levtype==depthcoord):
                ax1.set_ylabel('Depth [m]')
            

        ax1.grid('on')
        fig.colorbar(pc1,ax=ax1)

    fig.show()
    #code.interact(local=dict(globals(), **locals()))

    
def MinMaxCDateTime(tlists):

    nsets  = len(tlists)
    mindates = []
    maxdates = []
    for i in range(nsets):
        mindates.append(tlists[i][0])
        maxdates.append(tlists[i][-1])
        

    minid = np.argmin(mindates)
    maxid = np.argmax(maxdates)
    
    return tlists[minid][0], tlists[maxid][-1]


def PlotPFTSeries(tlist,hist_ds):

    # PLot Time-series data from an xarray dataset
    # There should be an axis break at split_id
    font = {'weight' : 'normal', 'size'   : 14}
    matplotlib.rc('font', **font)
    matplotlib.rc('axes.formatter', useoffset=False)

    #FATES_NPP_PF/FATES_NPLANT_PF
    #FATES_VEGC_PF
    
    #idc = [0,0,2,2,4,4,6,6,8,8,10,10,12,12,14,14]
    #idl = [0,1,0,1,0,1,0,1,0,1, 0, 1, 0, 1, 0, 1]

    #idc = [0,0,6,6,10,10]
    #idl = [0,1,0,1,0,1]

    # Good for just two
    colormap = plt.cm.tab20b
    idc = [0,2,9,15]
    #labels = ["A","B","C","D"]
    labels = ['non-fixer','fixer']
    #fig, axs = plt.subplots(ncols=1,nrows=2,figsize=(5,8))
    fig, axs = plt.subplots(ncols=2,nrows=3,figsize=(8,7))

    axs = axs.reshape(-1)
    
    megag_per_kg = 0.001
    m2_per_ha = 10000.0
    allom_agb_frac = 0.6
    #code.interact(local=dict(globals(), **locals()))
    ax=axs[0]
    #ax=axs
    numpft = len(hist_ds['fates_levpft'].data)
    for i in range(numpft):
        ax.plot(tlist.data,allom_agb_frac*hist_ds["FATES_VEGC_PF"].isel(fates_levpft=i).data, \
                label=labels[i],color=colormap.colors[idc[i]])

    ax.set_ylabel("Total Above \n Ground Biomass \n[kgC/m2]")
    ax.set_xlabel("Year")
    ax.legend()
    ax.grid('on')

    #da4 = UnstackCLSZPF(hist_ds["FATES_L2FR_CLSZPF"], \
    #                    hist_ds.dims['fates_levpft'], \
    #                    hist_ds.dims['fates_levscls'], \
    #                    hist_ds.dims['fates_levcan'])

    #"clscpf_ftmap"
    #sum(dim='clscpf_scmap','clscpf_clmap')
    #code.interact(local=dict(globals(), **locals()))
    
    #hist_ds["FATES_L2FR_PF1"] = da4.sel(clscpf_ftmap=1).mean(dim=('clscpf_scmap','clscpf_clmap'))
    #hist_ds["FATES_L2FR_PF2"] = da4.sel(clscpf_ftmap=2).mean(dim=('clscpf_scmap','clscpf_clmap'))

    #hist_ds["FATES_L2FR_USTORY_SZPF"] = da4.sel(clscpf_clmap=2).sel(clscpf_ftmap=1)
    
    hist_ds["FATES_L2FR_PF1"] = hist_ds["FATES_L2FR_CANOPY_REC_PF"].sel(fates_levpft=1)
    hist_ds["FATES_L2FR_PF2"] = hist_ds["FATES_L2FR_CANOPY_REC_PF"].sel(fates_levpft=2)
    
    #ax=axs[1]
    #FATES_RECL2FR_CANOPY_PF
    #FATES_RECL2FR_USTORY_PF
    ax=axs[1]
    #for i in range(numpft):
        #ax.plot(tlist.data,hist_ds["FATES_RECL2FR_USTORY_PF"].isel(fates_levpft=[i]).rolling(time=roll_w_size).mean(), \
        #        color=colormap.colors[idc[i]])
    ax.plot(tlist.data,hist_ds["FATES_L2FR_PF1"].rolling(time=roll_w_size).mean(),color=colormap.colors[idc[0]])
    ax.plot(tlist.data,hist_ds["FATES_L2FR_PF2"].rolling(time=roll_w_size).mean(),color=colormap.colors[idc[1]])

    #FATES_L2FR_CLSZPF
    #code.interact(local=dict(globals(), **locals()))   
    #sz_zeros = [i*len(hist_ds.fates_levpft) for i in range(len(hist_ds.fates_levscls))]
    
    #ax.plot(tlist.data,hist_ds["FATES_DDBH_USTORY_SZPF"].isel(fates_levscpf=sz_zeros))
    ax.set_ylabel(r"$\lambda$")
    #ax.set_ylabel(r"Understory Recruit $\lambda$")
    #ax.set_ylabel(r"Recruit Growth Increment")
    ax.set_xlabel("Year")
    ax.grid('on')
    
    hist_ds["FATES_CUE_PF"] = hist_ds['FATES_NPP_PF'] / hist_ds['FATES_GPP_PF']
    hist_ds["SYM_AFRAC"] = hist_ds["FATES_NFIX_SYM"] / (hist_ds["FATES_NH4UPTAKE"]+hist_ds["FATES_NO3UPTAKE"]+ hist_ds["FATES_NFIX_SYM"])
    hist_ds["SYM_PFRAC"] = hist_ds["FATES_NFIX_SYM"] / (0.001*hist_ds["NFIX_TO_SMINN"]+hist_ds["FATES_NFIX_SYM"])

    # 
    # FATES_SYMNFIX = 'kg m-2 s-1'
    
    ax=axs[2]
    ax.plot(tlist.data,hist_ds["SYM_AFRAC"].rolling(time=roll_w_size).mean(),color=colormap.colors[idc[2]])
    ax.set_ylabel("Symbiotic Fixation,\n Fraction of \nAcquisition ")
    ax.set_xlabel("Year")
    ax.grid('on')

    ax=axs[3]
    ax.plot(tlist.data,hist_ds["SYM_PFRAC"].rolling(time=roll_w_size).mean(),color=colormap.colors[idc[2]])
    ax.set_ylabel("Symbiotic Fixation,\n Fraction of \nGeneration")
    ax.set_xlabel("Year")
    ax.grid('on')

    ax=axs[4]
    hist_ds["SMIN_N"]=hist_ds["SMIN_NO3"]+hist_ds["SMIN_NH4"]
    ax.plot(tlist.data,hist_ds["SMIN_N"].rolling(time=roll_w_size).mean(),color=colormap.colors[idc[2]])
    ax.set_ylabel("Total Mineralized \nNitrogen in \nSoil [gN/m2]")
    ax.set_xlabel("Year")
    ax.grid('on')

    ax=axs[5]

    ax.plot(tlist.data,hist_ds["UNIT_NUPTAKE"].rolling(time=roll_w_size).mean()*sec_per_year,color=colormap.colors[idc[2]])
    ax.set_ylabel("Plant Unit \nN Uptake \n[gN/gC/year]")
    ax.set_xlabel("Year")
    ax.grid('on')
    
    #    varlist.append(varlist_type('Plant Unit N Uptake','[gN/gC/year]',hist_dsy["UNIT_NUPTAKE"].rolling(time=roll_w_size).mean()*sec_per_year  ))
    
    #varlist.append(varlist_type('CUE','[-]',hist_dsy["CUE"].rolling(time=roll_w_size).mean()))
    #varlist.append(varlist_type('Sym. N Fixation','[gN/m2/yr]',hist_dsy["FATES_SYMNFIX"].rolling(time=roll_w_size).mean()*1000.0*sec_per_year))
    #(hist_dsy["FATES_NH4UPTAKE"]+hist_dsy["FATES_NO3UPTAKE"])

    

    #FATES_DDBH_CANOPY_SZPF
    #FATES_DDBH_USTORY_SZPF

    
    plt.tight_layout()
    fig.show()
    code.interact(local=dict(globals(), **locals()))    


def PlotTSeries(tlists,vlists,labels): 


    # PLot Time-series data from an xarray dataset
    # There should be an axis break at split_id

    matplotlib.rc('axes.formatter', useoffset=False)

    nsets = len(vlists)
    nvars = len(vlists[0])

    # Choose a qualitative colormap

    #'Accent', 'Dark2', 'Set2'
    #code.interact(local=dict(globals(), **locals()))
    #colormap = plt.cm.Accent
    #colormap = plt.cm.Dark2
    #colormap = plt.cm.Pastel1
    colormap = plt.cm.tab20b
    
    if(nvars==1):
        nc = 1
        nr = 1
        fig, axs = plt.subplots(ncols=1,nrows=1,figsize=(4,4))
    elif(nvars==2):
        nc = 2
        nr = 1
        fig, axs = plt.subplots(ncols=2,nrows=1,figsize=(5,3))
    elif(nvars==4):
        nc = 2
        nr = 2
        fig, axs = plt.subplots(ncols=2,nrows=2,figsize=(5,4))
    elif(nvars>2):
        nc = 3
        nr = int(np.ceil(float(nvars)/3.0))
        fig, axs = plt.subplots(ncols=3,nrows=nr,figsize=(7.5,nr*2.1))
        
    ax1s = axs.reshape(-1)

    decyears = []
    for id in range(nsets):
        #code.interact(local=dict(globals(), **locals()))
        ntime=len(tlists[id].data)
        decyears.append([tlists[id].data[j].year+(tlists[id].data[j].dayofyr/365.0)+tlists[id].data[j].hour/(365.0*24.0) for j in range(ntime)])


    # SPINUP SETTINGS
    idc = [2,4,15,2,9]
    idl = [0,1,0,1]
    lstyle = ['solid','dashed','dotted']
    #cos = []
    #cos.append((0.2,0.2,0.6))
    #cos.append((0.6,0.2,0.2))
    #cos.append((0.6,0.6,0.6))
    #cos.append((0.6,0.6,0.9))
    #cos.append((0.9,0.6,0.6))

    
    #code.interact(local=dict(globals(), **locals()))
    cos = []
    #cos.append((colormap.colors[idc[0]]))
    #cos.append((colormap.colors[idc[1]]))
    #cos.append((colormap.colors[idc[2]]))

    cos.append((0.7,0.7,0.9))
    cos.append((0.3,0.3,0.7))
    cos.append((0.9,0.7,0.7))
    cos.append((0.7,0.3,0.3))
    
    
    #idc = [0,0,2,2,4,4,6,6,8,8,10,10,12,12,14,14]
    #idl = [0,1,0,1,0,1,0,1,0,1, 0, 1, 0, 1, 0, 1]

    #idc = [0,0,6,6,10,10]
    #idl = [0,1,0,1,0,1]

    # Good for just two 
    #idc = [0,8,14]
    #idc = [2,9,15,2,9,15,2,9,15]
    #idl = [0,0,0,1,1,1,2,2,2]
    
    #idc = [0,10,0,10]
    #idl = [0,0,1,1]

    
    
    #lstyle = ['solid','dashed','dotted']
    
    for iv in range(nvars):

        ax1 = ax1s[iv]
        irow = np.floor(float(iv+.01)/float(nc))+1
        ax1.set_ylabel(vlists[0][iv].units)

        if(vlists[0][iv].name[0]=='$'):
            ax1.set_title(r'{}'.format(vlists[0][iv].name))
        else:
            ax1.set_title(vlists[0][iv].name)


        vmax = -1e10
        vmin = 1e10
        for id in range(nsets):
            #code.interact(local=dict(globals(), **locals()))


            #data = monthly_to_annual(vlists[id][iv].data)
            #color=cos[id], \
                
            ax1.plot(decyears[id],vlists[id][iv].data, \
                     label=labels[id],color=cos[id], \
                     linestyle=lstyle[idl[id]])

            dat = vlists[id][iv].data
            if(not(vlists[id][iv].maxcap == None)):
                vmax = np.max([vmax,vlists[id][iv].maxcap])
            else:
                vmax = np.max([vmax,np.max(dat)+0.05*(np.max(dat)-np.min(dat))])
            
            if(not(vlists[id][iv].mincap == None)):
                vmin = np.min([vmin,vlists[0][iv].mincap])
            else:
                vmin = np.min([vmin,np.min(dat)-0.05*(np.max(dat)-np.min(dat))])
            
        ax1.set_xlim([years[0],years[1]])

        

        ax1.set_ylim([vmin,vmax])
        
        #colormap.colors[idc[id]], \
        if(irow!=nr):
            ax1.set_xticklabels('')
        else:
            ax1.set_xlabel('Year')
            
        ax1.grid('on')
        #ax1.set_xlim([mintime,maxtime])

        if(iv==0):
            ax1.legend()#loc='upper left')
            
        #if(iv==3):
        #    ax1.set_ylim([-0.5,0.5])

    
    #code.interact(local=dict(globals(), **locals()))

    #Turn off axes for unused panels
    for iv in range(nvars,len(ax1s)):
        #for id in range(nsets):
        #    ax1s[iv].plot(decyears[id],vlists[id][iv].data, \
        #                  label=labels[id],color=cos[id], \
        #                  linestyle=lstyle[idl[id]])
            
        #ax1s[iv].legend(loc='upper right') 
        ax1s[iv].set_axis_off()
        
        
    plt.tight_layout()
    fig.show()
    #code.interact(local=dict(globals(), **locals()))

    

def LoadHists(histfilepref):

    # First, lets find some commas

    histsplits=histfilepref.split(',')
    nfiles = len(histsplits)

    histds = []

    for i in range(nfiles):

        histfilepref=histsplits[i]
        
        # If the input path argument is a directory
        # try to use all the nc files in that directory
        if os.path.isdir(histfilepref):
            catstr = histfilepref+"/*nc"
            os.system("ncrcat {} {}".format(catstr,histfile))
        
        # If the input path argument is a single
        # file, then we assume that is the only file
        # we want, and it has already been concatenated
        elif os.path.isfile(histfilepref):
            # No concatenation necessary
            catstr = ''
            histfile = histfilepref

        # If the input path is neither a file or
        # a directory, it must be a partial path
        # string, and we will try just adding a wild-card
        else:
            catstr = histfilepref+"*nc"
            os.system("ncrcat {} {}".format(catstr,histfile))

        dstemp = xr.open_dataset(histfile)

        dstemp['time'] = dstemp.time_bounds.mean(dim="hist_interval")
        #code.interact(local=dict(globals(), **locals()))
        #histds.append(dstemp.sel(time=slice('0050-01-01','0500-01-01'))  )
        print('trimming time record, if necessary')
        histds.append(dstemp.where(dstemp['time.year'] >= years[0],drop=True).where(dstemp['time.year'] <= years[1], drop=True))
        print('done')
        
    return histds

    

# =======================================================================================
# This is the actual call to main

if __name__ == "__main__":
    main(sys.argv)
