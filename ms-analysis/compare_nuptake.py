import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import sys
import code  # For development: code.interact(local=locals())  code.interact(local=dict(globals(), **locals()))
import math
import os
from scipy import interpolate
from mpl_toolkits.axes_grid1 import make_axes_locatable

sec_per_year = 86400.0*365.0

font = {'weight' : 'normal', 'size'   : 13}
matplotlib.rc('font', **font)
matplotlib.rc('axes.formatter', useoffset=False)

    
# How much N do you get from paying 1kgC/day?

# Or how much does it cost to acquire 1 gN/m2/year

# Fixation should be more expensive per unit N, instead of transfering
# existing Nitrogen, it is sourcing it from a place that requires more energy

def UnitNfix(s_fix,a_fix,t_soil):


    # To fix N, you need to pay the surcharge on resp, but you also need to
    # to consider the sunk cost?
    
    # t_soil is in Celsius
    
    # N fixation parameters from Houlton et al (2008) and Fisher et al (2010)
    tfrz  = 274
    
    #a_fix = -3.62    # a parameter from Houlton et al. 2010 (a = -3.62 +/- 0.52)
    b_fix = 0.27     # b parameter from Houlton et al. 2010 (b = 0.27 +/-0.04)
    c_fix = 25.15    # c parameter from Houlton et al. 2010 (c = 25.15 +/- 0.66)

    # This is the unit carbon cost for nitrogen fixation. It is temperature dependant [kgC/kgN]
    unit_ccost_nfix = s_fix * (np.exp(a_fix + b_fix * t_soil * (1.0 - 0.5 * (t_soil) / c_fix)) - 2.)

    n_per_c_per_sec = 1/unit_ccost_nfix

    return(n_per_c_per_sec)


def UnitNUptake(vmax_n,tau_fr,t_soil):


    # tau_fr [years] turnover timescale

    gr_resp = 0.1
    q10_mr  = 1.5
    #fnrt_nc_ratio = 0.024
    fnrt_nc_ratio = 0.012
    red_factor = 1.0
    base_mr_20 = 2.52e-06  # Base MR rate @20C  gC/gN/s
    #tau_fr = 3.0
    
    # This actual = potential acquisition

    # MR

    tcsoil = q10_mr**((t_soil - 20.0)/10.0)
    fnrt_mr_rate = fnrt_nc_ratio * base_mr_20 * tcsoil * red_factor
    
    # Carbon required to maintain 1 kg of fine-root
    # Turnover replacement and respiration
    c_cost_per_sec = 1/(tau_fr*sec_per_year) + fnrt_mr_rate


    n_per_c_per_sec = vmax_n / c_cost_per_sec
    

    return(n_per_c_per_sec)
    

def main():



    # Create a list of temperatures


    t_soil = np.linspace(-10,50,61)
    n_cost_fix0 = np.zeros(t_soil.shape)
    n_cost_fix1 = np.zeros(t_soil.shape)
    n_cost_fix2 = np.zeros(t_soil.shape)
    
    n_cost_upt0 = np.zeros(t_soil.shape)
    n_cost_upt1 = np.zeros(t_soil.shape)
    n_cost_upt2 = np.zeros(t_soil.shape)
    
    vmax = [5.e-9,17.5e-9,15.e-9]
    tau_fr = [1.,4.,4.]
    
    #s_fix = -6.25    # s parameter from FUN model (fisher et al 2010)

    #s_fix = [-6.25,-23.0,-21.0]
    s_fix = [-6.25,-25.25]
    a_fix = [-3.62,-3.62]
    
    for i,t in enumerate(t_soil):
        n_cost_fix0[i] = UnitNfix(s_fix[0],a_fix[0],t)
        n_cost_fix1[i] = UnitNfix(s_fix[1],a_fix[1],t)
        #n_cost_fix2[i] = UnitNfix(s_fix[2],t)
        n_cost_upt0[i] = UnitNUptake(vmax[0],tau_fr[0],t)
        n_cost_upt1[i] = UnitNUptake(vmax[1],tau_fr[1],t)
        n_cost_upt2[i] = UnitNUptake(vmax[2],tau_fr[2],t)


    fig0 = plt.figure(figsize=(5.5, 4.5))
    plt.plot(t_soil,n_cost_fix0,label=r'fixation',color=[0.0,0.0,0.0])
    #plt.plot(t_soil,n_cost_fix0,label=r'fixation $a_{f1}$=-6.5',color=[0.0,0.0,0.0])
    #plt.plot(t_soil,n_cost_fix1,label=r'fixation $a_{f1}$=-25.0',color=[0.0,0.0,0.0],linestyle='--')
    plt.plot(t_soil,n_cost_upt0,label=r'uptake $\nu_{max}=5e^{-9}$, $\tau_{fr}$ = 1',color=[0.8,0.7,0.4])
    plt.plot(t_soil,n_cost_upt1,label=r'uptake $\nu_{max}=1.75e^{-8}$, $\tau_{fr}$ = 4',color=[0.8,0.7,0.4],linestyle='--')
    #plt.plot(t_soil,n_cost_upt2,label=r'uptake vmax=1.5e-8, $\tau_{fr}$ = 4',color=[0.2,0.2,0.2])
    
    plt.legend(fontsize=12)
    plt.grid('on')
    plt.ylim([0.0,0.7])
    plt.xlim([10,50])
    plt.xlabel('Temperature [C]')
    plt.ylabel('Potential Acquisition\nEfficiency [gN/gC]')
    plt.tight_layout()
    #plt.gca().set_rasterized(True)
    fig0.show()
    plt.savefig('fix_vs_uptake_cost_v3.eps')
    code.interact(local=dict(globals(), **locals()))    
    
    
if __name__ == "__main__":
    main()
