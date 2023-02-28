import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib
import sys
import code  # For development: code.interact(local=locals())  code.interact(local=dict(globals(), **locals()))
import argparse
from scipy.io import netcdf as nc
import cftime as cf
import os
from scipy import interpolate
import nc_time_axis
import pandas as pd

m2_per_ha = 10000.0
g_per_Mg  = 1000000.0
m2_per_cm2 = 1.0/(100.0*100.0)
numyears = 14   # Number of years in last set of data to evaluate (match met cycle)



# This script expects two file arguments, a parameter file and
# an output file.  It is expected that these files will be dimensioned
# by ensemble member, and that the order of the two is the same.
# The output from the history file will be compared against
# observations, and a fitness function will be calculated.

# Fitness function
# RRMSE of the mean values over a comparison window
# For:
#
# LAI
# AGB
# BA
# GPP
# Growth Increment
# Mortality Rate



# This class holds the ensemble dimensioned data
# for each parameter (instantiated for each parameter)
# and each darray* holds different ensemble members

class param_field:

    def __init__(self,name,optval,darrayf,darrayn,darrayo):

        self.name   = name
        self.optval = optval
        self.darrayf = darrayf  # Filtered out from meta analysis
        self.darrayn = darrayn  # Not in most fit group of members
        self.darrayo = darrayo  # Is in most fit group of members


# ========================================================================================
#                                        Main
# ========================================================================================

def main(argv):


    # Establish file paths
    # -----------------------------------------------------------------------------------
    
    histfile0 = 'Data/fates_e3sm_fullmodel_bci_parameter_ensemble_1pft_190329_multiinst_576inst_544db3b_0bc7a5d.h0.ensemble.sofar.nc'
    histfile1 = 'Data/fates_e3sm_fullmodel_bci_parameter_ensemble_1pft_190329_multiinst_576inst_544db3b_0bc7a5d.h1.ensemble.sofar.nc'
    
        
    paramfile = 'Data/fates_params_default_106ac7a_mod1PFT_exp1.ensemble.c190329.nc'
    lai_datafilename = 'Data/LAI.csv'
    gpp_datafilename = 'Data/gpp_gcm2y.txt'
    leaf_trait_spreadsheet = 'Data/Panama-TraitData-Ely/2016ENSO_Panama_CHN_data.csv'

    paramfile_out = 'Data/fates_params_default_106ac7a_mod1PFT_exp1.optimal.c201021.nc'

    
    # Perform an evaluation of existing trait data
    # which we can then use to filter acceptable ensemble members
    # -----------------------------------------------------------------------------------

    k_ely_lt_data = pd.read_csv(leaf_trait_spreadsheet)

    nc_all = 1.0/k_ely_lt_data['C:N']
    nc_all = nc_all[nc_all>0]

    sla_all = k_ely_lt_data['SLA_cm2_g']*m2_per_cm2
    sla_all = sla_all[sla_all>0]

    xlow = 0.15
    xhi  = 0.85

    nc_sorted = np.sort(nc_all)
    ilow_nc = int(len(nc_all)*xlow)
    ihi_nc  = int(len(nc_all)*xhi)
    
    sla_sorted = np.sort(sla_all)
    ilow_sla = int(len(sla_all)*xlow)
    ihi_sla  = int(len(sla_all)*xhi)

    sla_median = sla_sorted[int(len(sla_all)*0.5)]
    print(sla_median)
    
    
    fig3,ax2 = plt.subplots(ncols=1,nrows=1,figsize=(5,4))

    #frac = [float(i)/float(len(nc_all)) for i in range(len(nc_all))]
    #ax1.plot(frac,nc_sorted,color=[0.3,0.3,0.7])
    #ax1.axvline(x=xlow,color='k')
    #ax1.text(1.05*xlow,1.5*nc_sorted[ilow_nc],'{:.3f}\n[g/g]'.format(nc_sorted[ilow_nc]))
    #ax1.axvline(x=xhi,color='k')
    #ax1.text(1.05*xhi,0.8*nc_sorted[ihi_nc],'{:.3f}\n[g/g]'.format(nc_sorted[ihi_nc]))
    #ax1.set_title('Leaf N:C')
    #ax1.set_ylabel('[g/g]')
    #ax1.set_xlabel('rank')
    
    frac = [float(i)/float(len(sla_all)) for i in range(len(sla_all))]
    ax2.plot(frac,sla_sorted,color=[0.3,0.3,0.7])
    ax2.axvline(x=xlow,color='k')
    ax2.text(1.05*xlow,1.5*sla_sorted[ilow_sla],'{:.3f}\n[m2/g]'.format(sla_sorted[ilow_sla]))
    ax2.axvline(x=xhi,color='k')
    ax2.text(1.05*xhi,0.8*sla_sorted[ihi_sla],'{:.3f}\n[m2/g]'.format(sla_sorted[ihi_sla]))
    ax2.axvline(x=0.5,color='k')
    ax2.text(0.55,0.65*sla_sorted[int(len(sla_all)*0.5)],'{:.3f}\n[m2/g]'.format(sla_median))
    
    ax2.set_title('SLA')
    ax2.set_ylabel('[m2/g]')
    ax2.set_xlabel('rank')
    code.interact(local=dict(globals(), **locals()))

    nc_low = nc_sorted[ilow_nc]
    nc_hi = nc_sorted[ihi_nc]

    sla_low = sla_sorted[ilow_sla]
    sla_hi  = sla_sorted[ihi_sla]




    
    # Open the model ensemble output.  Assume that the output is already
    # Concatenated into 1 file, where each output variable has
    # a "record" dimension, as well as a time dimension
    # -----------------------------------------------------------------------------------


    # Raw un-modified datasets
    
    h0_ds = xr.open_dataset(histfile0)
    h1_ds = xr.open_dataset(histfile1)
    param_ds = xr.open_dataset(paramfile)
   
    #code.interact(local=dict(globals(), **locals()))

    # Trim the datasets down so we are only working with the last several years
    # Suggest using a number of years that matches your meteorological record

    endyear = int(h0_ds.time.dt.year[-1])+1
    begyear = endyear - numyears
    
    h0_ds14 = h0_ds.sel(time=h0_ds.time.dt.year.isin([i for i in range(begyear,endyear)]))
    h1_ds14 = h1_ds.sel(time=h1_ds.time.dt.year.isin([i for i in range(begyear,endyear)]))



    # Mask of records with reasonable C:N ratios (based on analysis of step 1)

    nc_rec_flt_in = np.logical_and(param_ds.fates_prt_nitr_stoich_p1.data[:,0,0] >= nc_low, \
                                   param_ds.fates_prt_nitr_stoich_p1.data[:,0,0] < nc_hi)
    
    nc_rec_flt_out = np.logical_or(param_ds.fates_prt_nitr_stoich_p1.data[:,0,0] < nc_low, \
                                   param_ds.fates_prt_nitr_stoich_p1.data[:,0,0] >= nc_hi)
    
    # Indices of records with reasonable C:N ratios
    nc_rec_idx_in  = np.where(nc_rec_flt_in)[0]
    nc_rec_idx_out = np.where(nc_rec_flt_out)[0]
    
    # Split the datasets into "in filter" and "out filter"
    h0_ds14if = h0_ds14.sel(record=nc_rec_flt_in)
    h0_ds14of = h0_ds14.sel(record=~nc_rec_flt_in)

    h1_ds14if = h1_ds14.sel(record=nc_rec_flt_in)
    h1_ds14of = h1_ds14.sel(record=~nc_rec_flt_in)


    n_ens = len(h1_ds14.record)
    n_ens_of = len(h1_ds14of.record)
    n_ens_if = len(h1_ds14if.record)
    


    
    # Load up the inventory benchmark data
    fp_ib = nc.netcdf_file('Data/census_bmks_lscaled_allyears_1pft_v2bci_201015.nc', 'r', mmap=False)

    
    # LAI
    # -----------------------------------------------------------------------------------
    
    # Get the observed LAI (BCI 01/09/2015 - 
    
    lai_data = pd.read_csv(lai_datafilename)
    lai_obs = lai_data['LAI '][:].mean()
    #lai_std = lai_data['LAI '][:].std()
    

    # Get the mean LAI (results in shape (ens,lat)
    
    lai_ens = h0_ds14if.TLAI.mean(dim="time").data       # Perform mean
    lai_ens = lai_ens.reshape(np.prod(lai_ens.shape))  # Compress to 1 dimension
    


    # AGB
    # -----------------------------------------------------------------------------------

    agb_obs = 136.0

    agb_ens = h0_ds14if.AGB.mean(dim="time").data * m2_per_ha / g_per_Mg
    agb_ens = agb_ens.reshape(np.prod(agb_ens.shape))


    
    # GPP
    # -----------------------------------------------------------------------------------

   
    gpp_gcm2y = np.ma.masked_array(np.loadtxt(gpp_datafilename))
    gpp_gcm2y = np.ma.masked_equal(gpp_gcm2y, 0.)
    gpp_obs = gpp_gcm2y.mean()   #Units: gC m-2 yr-1
    #gpp_std = gpp_gcm2y.std()

    gpp_ens = h0_ds14if.GPP.mean(dim="time").data * 86400 * 365
    gpp_ens = gpp_ens.reshape(np.prod(gpp_ens.shape))


    # For size-structured metrics, we will be using
    # some interpolations to match rates at the right size
    # -----------------------------------------------------------------------------------

    # Uppber bounds for shared/comparison size coordinates
    szu_shared = [30, 70]

    # Comparison points for mortality
    szm_shared = [5.5, 30] 

    # Comparison points for growth Increment
    szi_shared = [7.5, 12.5, 40]
     
    # Native size coordinate for ensemble members
    sz_ens  = h1_ds14if.fates_levscls.data[1:]

    # Size coordinate for inventory set
    sz_obs = fp_ib.variables['dclass'].data
        
    # Upper-bound points for size-classes in history ensemble
    szu_ens = np.zeros(len(sz_ens))
    for i in range(len(szu_ens)-1):
        szu_ens[i] = sz_ens[i+1]
    szu_ens[-1] = 1000.0

    # Upper-bound points for size-classes in observations
    szu_obs = np.zeros(len(sz_obs))
    for i in range(len(sz_obs)-1):
        szu_obs[i] = sz_obs[i+1]
    szu_obs[-1] = 1000.0

    
    # Mid-points for size classes in history ensemble
    szm_ens = np.zeros(len(sz_ens)-1)
    for i in range(len(szm_ens)):
        szm_ens[i] = 0.5*(sz_ens[i+1]+sz_ens[i])

    # Mid-points for size-classes in observations
    szm_obs = np.zeros(len(sz_obs)-1)
    for i in range(len(sz_obs)-1):
        szm_obs[i] = 0.5*(sz_obs[i+1]+sz_obs[i])
        
        
    # Basal Area
    # -----------------------------------------------------------------------------------

    # Mean over all census (censai?) (>1cm)
    basal_area_by_size_census = fp_ib.variables['basal_area_by_size_census'].data
    ba_obs = np.cumsum(basal_area_by_size_census[:,:,1].mean(axis=0))
    f_obs_ba_interp = interpolate.interp1d(szu_obs,ba_obs,'cubic')
    ba_obs_rgrd  = f_obs_ba_interp(szu_shared)
    
    ba_ens = h1_ds14if.BA_SCLS.mean(dim="time").data[:,1:,:]
    ba_ens = np.cumsum(ba_ens.reshape(ba_ens.shape[0],ba_ens.shape[1]),axis=1)

    f_ens_interp = interpolate.interp1d(szu_ens,ba_ens,'cubic')
    ba_ens_rgrd = f_ens_interp(szu_shared)




    # Growth Rate
    # Move these to the centers of their bins, and then interpolate to the shared coordinates
    # Trim off the highest rate, and then find midpoints
    # -----------------------------------------------------------------------------------


    # From Condit et al. 2017 [cm/yr]:
    #                      1995   2000  2005  2010  2015
    # Fitted Sap (6-8):     0.05,  0.058, 0.067, 0.078, 0.089 = 0.068
    # Fitted Mid (10-15):   0.14,  0.10,  0.101, 0.10,  0.116 = 0.1114
    # Fitted Large (30-50): 0.276, 0.283, 0.288, 0.294, 0.300 = 0.2882
    #
    # Mortality Rates %/yr
    #                       4.7    4.65   4.61   4.57   4.54  = 0.04614
    #                       2.49   2.63   2.78   2.92   3.07  = 0.02778
    #
    # Census (@7.5)  = 0.12427398
    #        (@12.5) = 0.201
    #        (@40)   = 0.512688
    
    #inc_obs      = fp_ib.variables['growth_increment_by_size_census'].data[1:,:,1].mean(axis=0)[:-1]

    inc_obs = [0.068, 0.1114, 0.2882]
    szi_obs = [7.5, 12.5, 40]
    
    #f_obs_inc_interp = interpolate.interp1d(szi_obs,inc_obs,'cubic')
    #inc_obs_rgrd = f_obs_inc_interp(szi_shared)

    inc_obs_rgrd = inc_obs
    
    inc_ds = (h1_ds14if.DDBH_CANOPY_SCLS + h1_ds14if.DDBH_UNDERSTORY_SCLS ) \
             / ( h1_ds14if.NPLANT_CANOPY_SCLS  + h1_ds14if.NPLANT_UNDERSTORY_SCLS)

    
    inc_ens = inc_ds.mean(dim="time").data[:,1:-1,:]  # Remove the first index because not in census and last
                                                      # Index because no midpoint
    inc_ens = inc_ens.reshape(inc_ens.shape[0],inc_ens.shape[1])

    inc_ens_rgrd = np.zeros([inc_ens.shape[0],len(szi_shared)])

    for i in range(inc_ens.shape[0]):
        mask = ~np.isnan(inc_ens[i,:])
        if(sum(mask)>0):
            f_ens_inc_interp = interpolate.interp1d(szm_ens[mask],inc_ens[i,mask],'cubic',bounds_error=False)
            inc_ens_rgrd[i,:] = f_ens_inc_interp(szi_shared)
        else:
            inc_ens_rgrd[i,:] = np.nan

    

    # Mortality Rate
    # -----------------------------------------------------------------------------------

    #mort_obs      = fp_ib.variables['mortality_rate_by_size_census'].data[1:,:,1].mean(axis=0)[:-1]

    mort_obs = [0.04614,  0.02778]
    szm_obs  = [5.5, 30]
    
    #f_obs_mort_interp = interpolate.interp1d(szm_obs,mort_obs,'cubic')
    #mort_obs_rgrd = f_obs_mort_interp(szm_shared)

    mort_obs_rgrd = mort_obs
    
    mort_ds = (h1_ds14if.MORTALITY_UNDERSTORY_SCLS + h1_ds14if.MORTALITY_CANOPY_SCLS) \
              /( h1_ds14if.NPLANT_CANOPY_SCLS  + h1_ds14if.NPLANT_UNDERSTORY_SCLS)
    
    
    #    code.interact(local=dict(globals(), **locals()))
    
    mort_ens = mort_ds.mean(dim="time").data[:,1:-1,:]
    mort_ens = mort_ens.reshape(mort_ens.shape[0],mort_ens.shape[1])

    #inc_obs_rgrd[1] = 0.068
    #inc_obs_rgrd[2] = 0.1114
    #inc_obs_rgrd[3] = 0.2882

    
    mort_ens_rgrd = np.zeros([mort_ens.shape[0],len(szm_shared)])
    for i in range(mort_ens.shape[0]):
        mask = ~np.isnan(mort_ens[i,:])
        if(sum(mask)>0):
            f_ens_mort_interp = interpolate.interp1d(szm_ens[mask],mort_ens[i,mask],'cubic',bounds_error=False)
            mort_ens_rgrd[i,:] = f_ens_mort_interp(szm_shared)
        else:
            mort_ens_rgrd[i,:] = np.nan



    # Calculate the differences and difference scores
            
    d_lai    = lai_ens - lai_obs
    d_gpp    = gpp_ens - gpp_obs
    d_agb    = agb_ens - agb_obs
    nd2_lai    = d_lai**2.0 / np.nanvar(d_lai)
    nd2_gpp    = d_gpp**2.0 / np.nanvar(d_gpp)
    nd2_agb    = d_agb**2.0 / np.nanvar(d_agb)
    
    d_inc  = np.zeros([n_ens_if,len(szi_shared)])
    d_mort = np.zeros([n_ens_if,len(szm_shared)])
    d_ba   = np.zeros([n_ens_if,len(szu_shared)])
    
    nd2_inc  = np.zeros([n_ens_if,len(szi_shared)])
    nd2_mort = np.zeros([n_ens_if,len(szm_shared)])
    nd2_ba   = np.zeros([n_ens_if,len(szu_shared)])
    
    for i in range(len(szi_shared)):
        d_inc[:,i]  = inc_ens_rgrd[:,i]  - inc_obs_rgrd[i]
        nd2_inc[:,i] = d_inc[:,i]**2.0 / np.nanvar(d_inc[:,i])
    
    for i in range(len(szm_shared)):
        d_mort[:,i] = mort_ens_rgrd[:,i] - mort_obs_rgrd[i]
        nd2_mort[:,i] = d_mort[:,i]**2.0 / np.nanvar(d_mort[:,i])
    
    for i in range(len(szu_shared)):
        d_ba[:,i]  = ba_ens_rgrd[:,i]  - ba_obs_rgrd[i]
        nd2_ba[:,i] = d_ba[:,i]**2.0 / np.nanvar(d_ba[:,i])


    #code.interact(local=dict(globals(), **locals()))
        
    # Find the member with the lowest difference score
    # -----------------------------------------------------------------------------------

    # REMOVING AGB !!!!
    
    nd2_all = nd2_lai+nd2_gpp + \
              np.sum(nd2_inc[:,:],axis=1)  + \
              np.sum(nd2_mort[:,:],axis=1) + \
              np.sum(nd2_ba[:,:],axis=1)

    # This is the id in the parameter file optimum
    opt_ens_fid = nc_rec_idx_in[np.nanargmin(nd2_all)]
    
    # This is the id in the dataset with optimum
    opt_ens_id  = np.nanargmin(nd2_all)
    

    
    # This is the list of the top 5% optimal
    nd2_allnn = np.zeros(nd2_all.shape)
    nd2_allnn = nd2_all[:]
    nd2_allnn[np.isnan(nd2_allnn)] = 100000.0
    
    ens_ord_ids = np.argsort(nd2_allnn)

    # Number of records in the top 5% optimized and filtered
    n5_opt_ens = int(0.05*len(nd2_all))


    filt_ens_ids      = nc_rec_idx_out[:]
    opt_ens_ord_ids   = nc_rec_idx_in[ens_ord_ids[:n5_opt_ens]]
    noopt_ens_ord_ids = nc_rec_idx_in[ens_ord_ids[n5_opt_ens:]]
    
    

    
    #code.interact(local=dict(globals(), **locals()))
    
    # Histograms, Raw means
    
    fig0,((ax1,ax2,ax3,ax4),(ax5,ax6,ax7,ax8),(ax9,ax10,ax11,ax12)) = plt.subplots(ncols=4,nrows=3,figsize=(10,7))

    
    ax1.hist(lai_ens,bins=20)
    ax1.set_title('LAI ')
    ax1.set_xlabel('m2/m2')
    ax1.axvline(x=lai_obs,color='k')
        
    ax2.hist(gpp_ens,bins=20)
    ax2.set_title('GPP ')
    ax2.set_xlabel('gC/m2/yr')
    ax2.axvline(x=gpp_obs,color='k')
    
    ax3.hist(agb_ens,bins=20)
    ax3.set_title('AGB ')
    ax3.set_xlabel('MgC/ha')
    ax3.axvline(x=agb_obs,color='k')
    
    ax4.axis("off")

    ax5.hist(mort_ens_rgrd[:,0],bins=20)
    ax5.set_title('Mortality @ 5.5cm')
    ax5.set_xlabel('/yr')
    ax5.set_xlim([0,0.35])
    ax5.axvline(x=mort_obs_rgrd[0],color='k')

    ax6.hist(mort_ens_rgrd[:,1],bins=20)
    ax6.set_title('Mortality @ 30cm')
    ax6.set_xlabel('/yr')
    ax6.axvline(x=mort_obs_rgrd[1],color='k')
    ax6.set_xlim([0,0.35])
    
    ax10.hist(np.log(inc_ens_rgrd[:,0]),bins=20)
    ax10.set_title('Growth Increment @ 7.5cm')
    ax10.set_xlabel('log(cm/yr)')
    ax10.axvline(x=np.log(inc_obs_rgrd[0]),color='k')
    ax10.set_xlim([-5,2])
    
    ax11.hist(np.log(inc_ens_rgrd[:,1]),bins=20)
    ax11.set_title('Growth Increment @ 12.5cm')
    ax11.set_xlabel('log(cm/yr)')
    ax11.axvline(x=np.log(inc_obs_rgrd[1]),color='k')
    ax11.set_xlim([-5,2])
    
    ax12.hist(np.log(inc_ens_rgrd[:,2]),bins=20)
    ax12.set_title('Growth Increment @ 40cm')
    ax12.set_xlabel('log(cm/yr)')
    ax12.axvline(x=np.log(inc_obs_rgrd[2]),color='k')
    ax12.set_xlim([-5,2])

    plt.tight_layout()

    
    # Rank plots

    #fig2,((ax1,ax2),(ax3,ax4),(ax5,ax6)) = plt.subplots(ncols=2,nrows=3,figsize=(8,7))
    
    #ax1.semilogy(np.sort(nd2_lai)[::-1])
    #ax1.set_title('Ranked Ensemble LAI Differences')
    #ax1.set_ylabel('s^2/var(s)')
    #ax2.semilogy(np.sort(nd2_agb)[::-1])
    #ax2.set_title('Ranked Ensemble AGB Differences')
    #ax2.set_ylabel('s^2/var(s)')
    #ax3.semilogy(np.sort(nd2_gpp)[::-1])
    #ax3.set_title('Ranked Ensemble GPP Differences')
    #ax3.set_ylabel('s^2/var(s)')
    #ax4.semilogy(np.sort(nd2_mort3)[::-1])
    #ax4.semilogy(np.sort(nd2_mort10)[::-1])
    #ax4.semilogy(np.sort(nd2_mort27)[::-1])
    #ax4.semilogy(np.sort(nd2_mort60)[::-1])
    #ax4.set_title('Ranked Ensemble Mortality Differences')
    #ax4.set_ylabel('s^2/var(s)')
    #ax5.semilogy(np.sort(nd2_inc3)[::-1])
    #ax5.semilogy(np.sort(nd2_inc10)[::-1])
    #ax5.semilogy(np.sort(nd2_inc27)[::-1])
    #ax5.semilogy(np.sort(nd2_inc60)[::-1])
    #ax5.set_title('Ranked Ensemble Increment Differences')
    #ax5.set_ylabel('s^2/var(s)')
    #ax6.semilogy(np.sort(nd2_all)[::-1])
    #ax6.set_title('Ranked Ensemble Combined Differences')
    #ax6.set_ylabel('s^2/var(s)')
    #plt.tight_layout()
   

    # -----------------------------------------------------------------------------------
    # Verification
    #    
    # Plots that evaluate the optimal member (sanity check and cross check)
    # Plot the BA distribution with size
    # -----------------------------------------------------------------------------------
    
   

    # Truncate ensemble data to match observation endpoint

    # Basal Area is already a cumulative sum

    for id in range(0): #n5_opt_ens):
    
        opt_ens_id = ens_ord_ids[id]

    
        n_sz_opt = np.sum(sz_ens < np.max(sz_obs))
        ba_opt = np.zeros(n_sz_opt+1)
        sz_opt = np.zeros(n_sz_opt+1)
        ba_opt[:n_sz_opt] = ba_ens[opt_ens_id,:n_sz_opt]
        sz_opt[:n_sz_opt] = sz_ens[:n_sz_opt]
        sz_opt[n_sz_opt]  = sz_obs[-1]
        ba_opt[n_sz_opt]  = ba_ens[opt_ens_id,-1]

    
        
        n_szi_opt = np.sum(szm_ens < np.max(szi_obs))
        szi_opt   = np.zeros(n_szi_opt + 1)
        szi_opt[:n_szi_opt] = szm_ens[:n_szi_opt]
        szi_opt[n_szi_opt] = szi_obs[-1]

        inc_opt   = np.zeros(n_szi_opt + 1)
        f_inc = interpolate.interp1d(szm_ens[:n_szi_opt+1],inc_ens[opt_ens_id,:n_szi_opt+1],'cubic',bounds_error=False)
        inc_opt[:] = f_inc(szi_opt)

        n_szm_opt = np.sum(szm_ens < np.max(szm_obs))
        szm_opt   = np.zeros(n_szm_opt + 1)
        szm_opt[:n_szm_opt] = szm_ens[:n_szm_opt]
        szm_opt[n_szm_opt] = szm_obs[-1]

        mort_opt = np.zeros(n_szm_opt + 1)
        f_mort = interpolate.interp1d(szm_ens[:n_szm_opt+1],mort_ens[opt_ens_id,:n_szm_opt+1],'cubic',bounds_error=False)
        mort_opt[:] = f_mort(szm_opt)


        fig2,((ax1,ax2),(ax3,ax4)) = plt.subplots(ncols=2,nrows=2,figsize=(7,6))

        ax1.scatter(sz_obs,ba_obs,label='Census',marker ='o',color = [0.4,0.4,0.4])
        ax1.plot(sz_opt,ba_opt,label='Simulation')
        ax1.set_ylabel('Cumulative Basal Area')
        ax1.legend(loc='upper left')
        ax1.grid('on')
        ax1.set_xlabel('DBH')
        ax1.set_ylim([0,1.05*np.nanmax([np.nanmax([ba_obs]),np.nanmax(ba_opt)])])
        
        ax2.scatter(szi_obs, inc_obs,marker ='o',color = [0.4,0.4,0.4])
        ax2.plot(szi_opt, inc_opt)
        ax2.set_ylabel('Mean Growth Increment [cm/yr]')
        ax2.grid('on')
        ax2.set_xlabel('DBH')
        ax2.set_ylim([0,1.05*np.nanmax([np.nanmax([inc_obs]),np.nanmax(inc_opt)])])
    
        ax3.scatter(szm_obs, mort_obs,marker ='o',color = [0.4,0.4,0.4])
        ax3.plot(szm_opt, mort_opt)
        ax3.set_ylabel('Mean Mortality Rate [/yr]')
        ax3.grid('on')
        ax3.set_xlabel('DBH')
        ax3.set_ylim([0,1.05*np.nanmax([np.nanmax([mort_obs]),np.nanmax(mort_opt)])])
        
        ax4.axis('off')
        ax4.text(0.15,0.65,'   Opt.      Obs.')
        ax4.text(0.05,0.5,'LAI    {:.2f}      {:.2f} [m2/m2]'.format(lai_ens[opt_ens_id],lai_obs))
        ax4.text(0.05,0.3,'GPP   {:.2f}      {:.2f} [kgC/m2/yr]'.format(gpp_ens[opt_ens_id]/1000.0,gpp_obs/1000.0))
        ax4.text(0.05,0.1,'AGB  {:.2f}   {:.2f} [kgC/ha]'.format(agb_ens[opt_ens_id],agb_obs))

        ax4.text(0.10,0.9,'Ensemble Member Id: {}'.format(opt_ens_id))
        
        plt.tight_layout()
    
        plt.show()
        

    
    # Identify the record dimension, and the pft dimension
    # Find variance over the record dimension, if any
    # save that as a histogram.  We keep everything
    # really simple and just use the first index of all
    # other dimensions to test variability


    # Break up the parameters into 3 groups
    # 1) Filtered out by data
    # 2) Inside filter requirement, but not in the top 5% optimized
    # 3) Top 5% optimized
    

    param_fields = []
    nbins = 20
    blacklist = ['fates_leaf_slamax','fates_allom_d2ca_coefficient_max']


    darray = np.zeros(n_ens)
    for varkey in param_ds.variables:
        
        var = param_ds.variables[varkey]
        ndims = len(var.dims)

        if(var.data.dtype == np.dtype('f')):
            if(ndims==1):
                darray[:] = var.data[:]                
            elif(ndims==2):
                darray[:] = var.data[:,0]
            elif(ndims==3):
                darray[:] = var.data[:,0,0]
            else:
                darray[:] = var.data[:,0,0,0]
                
            if(np.var(darray)>0 and not(varkey in blacklist)):
                print('parameter: {} has variance'.format(varkey))
                param_fields.append(param_field(varkey, nbins, \
                                                darray[filt_ens_ids], \
                                                darray[noopt_ens_ord_ids], \
                                                darray[opt_ens_ord_ids]))


                

    # Histogram plots of the parameters from the ensemble, breaking
    # it up into filtered, non-optimal, and 5% optimal groups

    fig3, axs = plt.subplots(ncols=4,nrows=3,figsize=(10,9))

    axs = axs.reshape(-1)
    
    for vid in range(len(axs)):

        if(vid<len(param_fields)):
            axs[vid].hist([param_fields[vid].darrayf, \
                           param_fields[vid].darrayn, \
                           param_fields[vid].darrayo], \
                          bins=20,stacked=True,color= [[0.8,0.8,0.8],[0.5,0.5,0.5],[0.3,0.3,0.8]])
            axs[vid].set_title(param_fields[vid].name.split('fates_')[1])
        else:
            axs[vid].axis("off")

    plt.tight_layout()


    # co-variance scatter-plots of parameter combinations

    
    nparms = len(param_fields)
    fig4, axs = plt.subplots(ncols = nparms,nrows = nparms, figsize=(10,10))
    fig4.subplots_adjust(hspace=0,wspace=0,left=0.05,right=0.95,bottom=0.05,top=0.95)
    
    axs = axs.reshape(-1)

    for p1 in range(nparms):

        vid = p1*nparms
        axs[vid].axis("off")
        
        for p2 in range(1,p1+1):
            vid = p1*nparms+p2
            axs[vid].axis("off")
            
        for p2 in range(p1+1,nparms):
            vid = p1*nparms+p2
            axs[vid].scatter(param_fields[p1].darrayf,param_fields[p2].darrayf,marker ='o',color=[0.8,0.8,0.8],s=4)
            axs[vid].scatter(param_fields[p1].darrayn,param_fields[p2].darrayn,marker ='o',color=[0.5,0.5,0.5],s=4)
            axs[vid].scatter(param_fields[p1].darrayo,param_fields[p2].darrayo,marker ='o',color=[0.3,0.3,0.8],s=4)
            axs[vid].set_xticklabels([])
            axs[vid].set_yticklabels([])
            
    #plt.tight_layout()

    for p1 in range(nparms):

        vid = p1*nparms
        axs[vid].text(0.0,0.25,param_fields[p1].name.split('fates_')[1], \
                      size=9,rotation=-35,horizontalalignment='left',verticalalignment='center')

    
    plt.show()



    # Write out a new parameter file with just one record
    
    opt_ens_id = nc_rec_idx_in[224]     #ens_ord_ids[id]

    param_ds_out = param_ds.sel(record=opt_ens_id)


    
    code.interact(local=dict(globals(), **locals()))
    
    param_ds_out.to_netcdf(paramfile_out,mode='w')

    
    exit(0)
     
# =======================================================================================
# This is the actual call to main

if __name__ == "__main__":
    main(sys.argv)
