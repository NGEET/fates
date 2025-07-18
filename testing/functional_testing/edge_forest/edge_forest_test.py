"""
Concrete class for running the allometry functional tests for FATES.
"""

import os
import xarray as xr
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from utils import round_up
from utils_plotting import sample_colormap, blank_plot
from functional_class import FunctionalTest


class EdgeForestTest(FunctionalTest):
    """Quadratic test class"""

    name = "edge_forest"

    def __init__(self, test_dict):
        super().__init__(
            EdgeForestTest.name,
            test_dict["test_dir"],
            test_dict["test_exe"],
            test_dict["out_file"],
            test_dict["use_param_file"],
            test_dict["other_args"],
        )
        self.plot = True

    def plot_output(self, run_dir: str, save_figs: bool, plot_dir: str):
        """Plots all edge forest plots

        Args:
            run_dir (str): run directory
            out_file (str): output file name
            save_figs (bool): whether or not to save the figures
            plot_dir (str): plot directory to save the figures to
        """

        # read in edge forest data
        edge_forest_dat = xr.open_dataset(os.path.join(run_dir, self.out_file))

        # Plot all bins
        da = edge_forest_dat["frac_in_every_bin"]
        da_norm = edge_forest_dat["frac_in_every_bin_norm"]
        self.plot_edge_forest_frac_allbins(
            da,
            da_norm,
            "Fraction of forest in each edge bin",
            da.attrs["units"],
            save_figs,
            plot_dir,
        )

        # Plot individual bins
        plot_dict = {
            "frac_in_bin_gaussian": {
                "varname": "Fraction of forest in first bin with Gaussian fit",
                "units": "unitless",
            },
            "frac_in_bin_lognormal": {
                "varname": "Fraction of forest in first bin with lognormal fit",
                "units": "unitless",
            },
            "frac_in_bin_quadratic": {
                "varname": "Fraction of forest in first bin with quadratic fit",
                "units": "unitless",
            },
        }
        for plot, attributes in plot_dict.items():
            self.plot_edge_forest_frac_onebin(
                edge_forest_dat[plot],
                attributes["varname"],
                attributes["units"],
                save_figs,
                plot_dir,
            )

    @staticmethod
    def plot_edge_forest_frac_allbins(
        data: xr.Dataset, data_norm: xr.Dataset, varname: str, units: str, save_fig: bool, plot_dir: str = None
    ):
        """Plot the fraction of forest in all bins

        Args:
            data (xarray DataArray): the data array of the variable to plot
            data_norm (xarray DataArray): the data array of the normalized variable to plot
            var (str): variable name (for data structure)
            varname (str): variable name for plot labels
            units (str): variable units for plot labels
            save_fig (bool): whether or not to write out plot
            plot_dir (str): if saving figure, where to write to
        """

        x = data["frac_forest"]
        max_x = x.max()

        y = data
        y_norm = data_norm
        max_y = round_up(y.max())

        blank_plot(max_x, 0.0, max_y, 0.0, draw_horizontal_lines=True)

        bins = data["bin"].values

        for b, this_bin in enumerate(bins):
            plt.plot(
                x.values,
                y.sel(bin=this_bin).values,
                color=sample_colormap("jet_r", len(bins), b),
                label=str(this_bin),
            )

        for b, this_bin in enumerate(bins):
            plt.plot(
                x.values,
                y_norm.sel(bin=this_bin).values,
                color=sample_colormap("jet_r", len(bins), b),
                label="__nolegend__",
                linestyle="--",
            )

        plt.xlabel("Frac. forest in site", fontsize=11)
        plt.ylabel(f"{varname} ({units})", fontsize=11)
        plt.title(f"Simulated {varname} for input parameter file\n(dashed=adjusted)", fontsize=11)
        plt.legend(loc="best", title="Bin distance to\nnonforest (m)", alignment="left")

        if save_fig:
            fig_name = os.path.join(plot_dir, f"edge_forest_plot_{data.name}.png")
            plt.savefig(fig_name)

    @staticmethod
    def plot_edge_forest_frac_onebin(
        data: xr.Dataset, varname: str, units: str, save_fig: bool, plot_dir: str = None
    ):
        """Plot the fraction of forest in a given bin

        Args:
            data (xarray DataArray): the data array of the variable to plot
            var (str): variable name (for data structure)
            varname (str): variable name for plot labels
            units (str): variable units for plot labels
            save_fig (bool): whether or not to write out plot
            plot_dir (str): if saving figure, where to write to
        """

        # This is left over from AllometryTest, which had two dimensions
        data_frame = pd.DataFrame(
            {
                "frac_forest": np.tile(data.frac_forest, 1),
                data.name: data.values.flatten(),
            }
        )

        max_forest = data_frame["frac_forest"].max()
        max_var = round_up(data_frame[data.name].max())

        blank_plot(max_forest, 0.0, max_var, 0.0, draw_horizontal_lines=True)

        plt.plot(
            data_frame.frac_forest.values,
            data_frame[data.name].values,
            lw=2,
        )

        plt.xlabel("Frac. forest in site", fontsize=11)
        plt.ylabel(f"{varname} ({units})", fontsize=11)
        plt.title(f"Simulated {varname} for input parameter file", fontsize=11)

        if save_fig:
            fig_name = os.path.join(plot_dir, f"edge_forest_plot_{data.name}.png")
            plt.savefig(fig_name)
