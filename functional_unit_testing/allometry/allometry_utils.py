"""Utility functions for allometry functional unit tests
"""
import os
import math
import pandas as pd
import numpy as np
import xarray as xr
import matplotlib
import matplotlib.pyplot as plt
from utils import get_color_pallete, round_up

def plot_allometry_var(data, varname, units, save_fig, plot_dir=None):
    """Plot an allometry variable

    Args:
        data (xarray DataArray): the data array of the variable to plot
        var (str): variable name (for data structure)
        varname (str): variable name for plot labels
        units (str): variable units for plot labels
        save_fig (bool): whether or not to write out plot
        plot_dir (str): if saving figure, where to write to
    """
    data_frame = pd.DataFrame({'dbh': np.tile(data.dbh, len(data.pft)),
                       'pft': np.repeat(data.pft, len(data.dbh)),
                       data.name: data.values.flatten()})

    maxdbh = data_frame['dbh'].max()
    maxvar = round_up(data_frame[data.name].max())

    colors = get_color_pallete()

    plt.figure(figsize=(7, 5))
    axis = plt.subplot(111)
    axis.spines["top"].set_visible(False)
    axis.spines["bottom"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["left"].set_visible(False)

    axis.get_xaxis().tick_bottom()
    axis.get_yaxis().tick_left()

    plt.xlim(0.0, maxdbh)
    plt.ylim(0.0, maxvar)

    plt.yticks(fontsize=10)
    plt.xticks(fontsize=10)

    inc = (int(maxvar) - 0)/20
    for i in range(0, 20):
        y_val = 0.0 + i*inc
        plt.plot(range(math.floor(0), math.ceil(maxdbh)),
                  [y_val] * len(range(math.floor(0), math.ceil(maxdbh))),
                  "--", lw=0.5, color="black", alpha=0.3)

    plt.tick_params(bottom=False, top=False, left=False, right=False)

    pfts = np.unique(data_frame.pft.values)
    for rank, pft in enumerate(pfts):
        dat = data_frame[data_frame.pft == pft]
        plt.plot(dat.dbh.values, dat[data.name].values, lw=2, color=colors[rank],
                 label=pft)

    plt.xlabel('DBH (cm)', fontsize=11)
    plt.ylabel(f'{varname} ({units})', fontsize=11)
    plt.title(f"Simulated {varname} for input parameter file", fontsize=11)
    plt.legend(loc='upper left', title='PFT')

    if save_fig:
        fig_name = os.path.join(plot_dir, f"allometry_plot_{data.name}.png")
        plt.savefig(fig_name)

def plot_total_biomass(data, save_fig, plot_dir):
    """Plot two calculations of total biomass against each other

    Args:
        data (xarray DataSet): the allometry dataset
    """
    data_frame = pd.DataFrame({'dbh': np.tile(data.dbh, len(data.pft)),
                       'pft': np.repeat(data.pft, len(data.dbh)),
                       'total_biomass_parts': data.total_biomass_parts.values.flatten(),
                       'total_biomass_tissues': data.total_biomass_tissues.values.flatten()})

    colors = get_color_pallete()

    plt.figure(figsize=(7, 5))
    axis = plt.subplot(111)
    axis.spines["top"].set_visible(False)
    axis.spines["bottom"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["left"].set_visible(False)

    axis.get_xaxis().tick_bottom()
    axis.get_yaxis().tick_left()

    maxbiomass = np.maximum(data_frame['total_biomass_parts'].max(),
                            data_frame['total_biomass_tissues'].max())

    plt.xlim(0.0, maxbiomass)
    plt.ylim(0.0, maxbiomass)

    plt.yticks(fontsize=10)
    plt.xticks(fontsize=10)
    plt.tick_params(bottom=False, top=False, left=False, right=False)

    pfts = np.unique(data_frame.pft.values)
    for rank, pft in enumerate(pfts):
        data = data_frame[data_frame.pft == pft]
        plt.scatter(data.total_biomass_parts.values, data.total_biomass_parts.values,
                 color=colors[rank], label=pft)

    plt.xlabel('Total biomass (kgC) from parts', fontsize=11)
    plt.ylabel('Total biomass (kgC) from tissues', fontsize=11)
    plt.title("Simulated total biomass for input parameter file", fontsize=11)
    plt.legend(loc='upper left', title='PFT')

    if save_fig:
        fig_name = os.path.join(plot_dir, "allometry_plot_total_biomass_compare.png")
        plt.savefig(fig_name)

def plot_allometry_dat(run_dir, out_file, save_figs, plot_dir):
    """Plots all allometry plots

    Args:
        run_dir (str): run directory
        out_file (str): output file name
        save_figs (bool): whether or not to save the figures
        plot_dir (str): plot directory to save the figures to
    """

    # read in allometry data
    allometry_dat = xr.open_dataset(os.path.join(run_dir, out_file))

    plot_dict = {
      'height': {
        'varname': 'height',
        'units': 'm',
        },
      'bagw': {
        'varname': 'aboveground biomass',
        'units': 'kgC',
        },
      'blmax': {
        'varname': 'maximum leaf biomass',
        'units': 'kgC',
        },
      'crown_area': {
        'varname': 'crown area',
        'units': 'm$^2$',
        },
      'sapwood_area': {
        'varname': 'sapwood area',
        'units': 'm$^2$',
        },
      'bsap': {
        'varname': 'sapwood biomass',
        'units': 'kgC',
        },
      'bbgw': {
        'varname': 'belowground biomass',
        'units': 'kgC',
        },
      'fineroot_biomass': {
        'varname': 'fineroot biomass',
        'units': 'kgC',
        },
      'bstore': {
        'varname': 'storage biomass',
        'units': 'kgC',
        },
      'bdead': {
        'varname': 'deadwood biomass',
        'units': 'kgC',
        },
      'total_biomass_parts': {
        'varname': 'total biomass (calculated from parts)',
        'units': 'kgC',
        },
      'total_biomass_tissues': {
        'varname': 'total biomass (calculated from tissues)',
        'units': 'kgC',
        },

    }
    for plot, attributes in plot_dict.items():
        plot_allometry_var(allometry_dat[plot], attributes['varname'],
                           attributes['units'], save_figs, plot_dir)

    plot_total_biomass(allometry_dat, save_figs, plot_dir)
