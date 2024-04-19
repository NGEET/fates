"""Utility functions for allometry functional unit tests
"""
import os
import math
import xarray as xr
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from utils import get_color_pallete

def plot_quadratic_dat(run_dir, out_file, save_figs, plot_dir):

    # read in quadratic data
    quadratic_dat = xr.open_dataset(os.path.join(run_dir, out_file))

    # plot output
    PlotQuadAndRoots(quadratic_dat.a.values, quadratic_dat.b.values,
                     quadratic_dat.c.values, quadratic_dat.root1.values,
                     quadratic_dat.root2.values)

def PlotQuadAndRoots(a, b, c, r1, r2):

    colors = get_color_pallete()

    fig, axs = plt.subplots(ncols=1, nrows=1, figsize=(8,8))
    x = np.linspace(-10.0, 10.0, num=20)

    for i in range(0, len(a)):
      y = a[i]*x**2 + b[i]*x + c[i]
      plt.plot(x, y, lw=2, color=colors[i])
      plt.scatter(r1[i], r2[i], color=colors[i], s=50)
      plt.axhline(y=0.0, color='k', linestyle='dotted')
