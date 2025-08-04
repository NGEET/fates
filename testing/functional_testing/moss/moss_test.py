"""
Concrete class for running the moss functional tests for FATES.
"""

import os
import xarray as xr
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from utils import round_up, truncate
from utils_plotting import blank_plot, get_color_palette
from functional_class import FunctionalTest


class MossTest(FunctionalTest):
    """Moss test class"""

    name = "moss"

    def __init__(self, test_dict):
        super().__init__(
            MossTest.name,
            test_dict["test_dir"],
            test_dict["test_exe"],
            test_dict["out_file"],
            test_dict["use_param_file"],
            test_dict["other_args"],
        )
        self.plot = True

    def plot_output(self, run_dir: str, save_figs: bool, plot_dir: str):
        """Plots all moss plots

        Args:
            run_dir (str): run directory
            out_file (str): output file name
            save_figs (bool): whether or not to save the figures
            plot_dir (str): plot directory to save the figures to
        """
        print("Making plots")

        # read in moss data
        moss_dat = xr.open_dataset(os.path.join(run_dir, self.out_file))

        # Make plots from 1-d variables
        var_list = ["out_dlgf"]
        for var in var_list:
            self.plot_1d(
                moss_dat[var],
                save_figs,
                plot_dir,
            )

        # Make plots from variables with outputs dimensioned: leaf area index x moss biomass (before)
        var_list = ["out_al", "out_algf"]
        for var in var_list:
            self.plot_moss_lai_x_mossbiomass(
                moss_dat[var],
                save_figs,
                plot_dir,
            )

        # Make plots from variables with outputs dimensioned: assimilation x moss biomass (before)
        var_list = ["out_resp", "out_mort", "out_biomass_delta"]
        for var in var_list:
            self.plot_moss_assim_x_mossbiomass(
                moss_dat[var],
                save_figs,
                plot_dir,
            )

    @staticmethod
    def plot_1d(
        data: xr.DataArray, save_fig: bool, plot_dir: str
    ):
        """Plot a 1-d variable

        Args:
            data (xarray DataArray): the data array of the variable to plot
            save_fig (bool): whether or not to write out plot
            plot_dir (str): if saving figure, where to write to
        """
        dim0 = data.dims[0]

        plt.plot(
            data[dim0].values,
            data.values,
            lw=2,
        )

        x_label = f"{data[dim0].attrs['long_name']} ({data[dim0].attrs['units']})"
        varname = data.attrs["long_name"]
        y_label = f"{varname} ({data.attrs['units']})"

        plt.xlabel(x_label, fontsize=11)
        plt.ylabel(y_label, fontsize=11)
        plt.title(f"Simulated {varname}", fontsize=11)

        if save_fig:
            fig_name = os.path.join(plot_dir, f"moss_plot_{data.name}.png")
            plt.savefig(fig_name)
        plt.close()

    def plot_moss_lai_x_mossbiomass(
        self, *args,
    ):
        """
        Plot a variable with outputs dimensioned: leaf area index (lai) x moss biomass.

        This method takes a variable (typically an xarray DataArray) that is dimensioned by
        leaf area index and moss biomass, and generates a plot where the X axis is leaf area index
        and the legend represents different moss biomass values.

        Args:
            *args: Arguments passed to plot_dim0x_dim1legend, typically including:
                data (xarray.DataArray): The variable to plot.
                save_fig (bool): Whether or not to save the figure.
                plot_dir (str): Directory to save the figure to.

        Returns:
            None. Displays or saves the plot depending on the save_fig argument.
        """
        self.plot_dim0x_dim1legend(*args, dim0="lai", dim1="moss_biomass")

    def plot_moss_assim_x_mossbiomass(
        self, *args,
    ):
        """
        Plot a variable with outputs dimensioned: assimilation (assim) x moss biomass.

        This method takes a variable (typically an xarray DataArray) that is dimensioned by
        moss assimilation and moss biomass, and generates a plot where the X axis is assimilation
        and the legend represents different moss biomass values.

        Args:
            *args: Arguments passed to plot_dim0x_dim1legend, typically including:
                data (xarray.DataArray): The variable to plot.
                save_fig (bool): Whether or not to save the figure.
                plot_dir (str): Directory to save the figure to.

        Returns:
            None. Displays or saves the plot depending on the save_fig argument.
        """
        self.plot_dim0x_dim1legend(*args, dim0="assim", dim1="moss_biomass")

    @staticmethod
    def plot_dim0x_dim1legend(
        data: xr.DataArray, save_fig: bool, plot_dir: str, dim0: str, dim1: str
    ):
        """Plot a variable with outputs dimensioned: (X axis, legend item)

        Args:
            data (xarray DataArray): the data array of the variable to plot
            save_fig (bool): whether or not to write out plot
            plot_dir (str): if saving figure, where to write to
        """
        data_frame = pd.DataFrame(
            {
                dim0: np.tile(data[dim0], len(data[dim1])),
                dim1: np.repeat(data[dim1], len(data[dim0])),
                data.name: data.values.flatten(),
            }
        )

        legend_values = np.unique(data_frame[dim1].values)
        colors = get_color_palette(len(legend_values))
        for rank, legend_value in enumerate(legend_values):
            dat = data_frame[data_frame.moss_biomass == legend_value]
            plt.plot(
                dat[dim0].values,
                dat[data.name].values,
                lw=2,
                color=colors[rank],
                label=legend_value,
            )

        x_label = f"{data[dim0].attrs['long_name']} ({data[dim0].attrs['units']})"
        legend_label = f"{data[dim1].attrs['long_name']} ({data[dim1].attrs['units']})"
        varname = data.attrs["long_name"]
        y_label = f"{varname} ({data.attrs['units']})"

        plt.xlabel(x_label, fontsize=11)
        plt.ylabel(y_label, fontsize=11)
        plt.title(f"Simulated {varname}", fontsize=11)
        plt.legend(loc="upper left", title=legend_label)

        if save_fig:
            fig_name = os.path.join(plot_dir, f"moss_plot_{data.name}.png")
            plt.savefig(fig_name)
        plt.close()
