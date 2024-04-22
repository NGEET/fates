"""Utility functions for allometry functional unit tests
"""
import os
import math
import pandas as pd
import numpy as np
import xarray as xr
import matplotlib
import matplotlib.pyplot as plt
from utils import get_color_palette, round_up

def blank_plot(x_max, x_min, y_max, y_min, draw_horizontal_lines=False):
    """Generate a blank plot with set attributes

    Args:
        x_max (float): maximum x value
        x_min (float): minimum x value
        y_max (float): maximum y value
        y_min (float): minimum y value
        draw_horizontal_lines (bool, optional): whether or not to draw horizontal
        lines across plot. Defaults to False.
    """

    plt.figure(figsize=(7, 5))
    axis = plt.subplot(111)
    axis.spines["top"].set_visible(False)
    axis.spines["bottom"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["left"].set_visible(False)

    axis.get_xaxis().tick_bottom()
    axis.get_yaxis().tick_left()

    plt.xlim(0.0, x_max)
    plt.ylim(0.0, y_max)

    plt.yticks(fontsize=10)
    plt.xticks(fontsize=10)

    if draw_horizontal_lines:
        inc = (int(y_max) - y_min)/20
        for i in range(0, 20):
            plt.plot(range(math.floor(x_min), math.ceil(x_max)),
                      [0.0 + i*inc] * len(range(math.floor(x_min), math.ceil(x_max))),
                      "--", lw=0.5, color="black", alpha=0.3)

    plt.tick_params(bottom=False, top=False, left=False, right=False)

    return plt

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

    max_dbh = data_frame['dbh'].max()
    max_var = round_up(data_frame[data.name].max())

    blank_plot(max_dbh, 0.0, max_var, 0.0, draw_horizontal_lines=True)

    pfts = np.unique(data_frame.pft.values)
    colors = get_color_palette(len(pfts))
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

    max_biomass = np.maximum(data_frame['total_biomass_parts'].max(),
                            data_frame['total_biomass_tissues'].max())

    blank_plot(max_biomass, 0.0, max_biomass, 0.0, draw_horizontal_lines=False)

    pfts = np.unique(data_frame.pft.values)
    colors = get_color_palette(len(pfts))
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
