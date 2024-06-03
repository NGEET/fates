"""Utility functions for fuel functional unit tests
"""
import os
import math
import pandas as pd
import numpy as np
import xarray as xr
import matplotlib
import matplotlib.pyplot as plt

def plot_fuel_dat(run_dir, out_file, save_figs, plot_dir):
    """Plot output associated with fuel tests

    Args:
        run_dir (str): run directory
        out_file (str): output file
        save_figs (bool): whether or not to save the figures
        plot_dir (str): plot directory
    """

    fuel_dat = xr.open_dataset(os.path.join(run_dir, out_file))
    
    plot_NI_dat(fuel_dat, save_figs, plot_dir)
    plot_barchart(fuel_dat, 'fuel_loading', 'Fuel loading', 'kgC m$^{-2}$', save_figs, plot_dir)
    plot_barchart(fuel_dat, 'frac_loading', 'Fractional fuel loading', '0-1', save_figs, plot_dir)
    
def plot_barchart(fuel_dat, var, varname, units, save_figs, plot_dir):
    litter_classes = ['twigs', 'small branches', 'large branches', 'dead leaves', 'live grass']

    fuel_models = ['Fuel model ' + str(m) for m in np.unique(fuel_dat.fuel_model)]
    data_dict = {'{}'.format(lc): fuel_dat.isel(litter_class=i)[var].values for i, lc in enumerate(litter_classes)}
    
    fig, ax = plt.subplots()
    bottom = np.zeros(len(fuel_models))
    for litter_class, dat in data_dict.items():
        p = ax.bar(fuel_models, dat, 0.5, label=litter_class, bottom=bottom)
        bottom += dat
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.75, box.height])
    plt.ylabel(f'{varname} ({units})', fontsize=11)
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    
def plot_NI_dat(fuel_dat, save_figs, plot_dir):
    """Plot output for Nesterov index

    Args:
        fuel_dat (Xarray Dataset): output fuel data
        save_figs (bool): whether or not to save the figures
        plot_dir (str): plot directory
    """
    
    plt.figure()
    fuel_dat.NI.plot()
    plt.xlabel('Time', fontsize=11)
    plt.ylabel('Nesterov Index', fontsize=11)
    
    if save_figs:
        fig_name = os.path.join(plot_dir, "Nesterov_plot.png")
        plt.savefig(fig_name)