"""Utility functions for allometry functional unit tests
"""
import os
import math
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

from utils import get_color_pallete

def plot_quadratic_dat(run_dir, out_file, save_figs, plot_dir):
    """Reads in and plots quadratic formula test output

    Args:
        run_dir (str): run directory
        out_file (str): output file
        save_figs (bool): whether or not to save the figures
        plot_dir (str): plot directory
    """

    # read in quadratic data
    quadratic_dat = xr.open_dataset(os.path.join(run_dir, out_file))

    # plot output
    plot_quad_and_roots(quadratic_dat.a.values, quadratic_dat.b.values,
                     quadratic_dat.c.values, quadratic_dat.root1.values,
                     quadratic_dat.root2.values)
    if save_figs:
      fig_name = os.path.join(plot_dir, "quadratic_test.png")
      plt.savefig(fig_name)

def plot_quad_and_roots(a_coeff, b_coeff, c_coeff, root1, root2):
    """Plots a set of quadratic formulas (ax**2 + bx + c) and their two roots

    Args:
        a (float array): set of a coefficients
        b (float array): set of b coefficients
        c (float array): set of b coefficients
        r1 (float array): set of first real roots
        r2 (float array): set of second real roots
    """

    colors = get_color_pallete()

    plt.figure(figsize=(7, 5))
    x_vals = np.linspace(-10.0, 10.0, num=20)

    for i in range(len(a_coeff)):
        y_vals = a_coeff[i]*x_vals**2 + b_coeff[i]*x_vals + c_coeff[i]
        plt.plot(x_vals, y_vals, lw=2, color=colors[i])
        plt.scatter(root1[i], root2[i], color=colors[i], s=50)
        plt.axhline(y=0.0, color='k', linestyle='dotted')
