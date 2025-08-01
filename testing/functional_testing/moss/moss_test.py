"""
Concrete class for running the moss functional tests for FATES.
"""

import os
import xarray as xr
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from utils import round_up, blank_plot, get_color_palette
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

        # read in moss data
        moss_dat = xr.open_dataset(os.path.join(run_dir, self.out_file))

        # Make plots from variables with outputs dimensioned: cumulative leaf area x moss biomass
        var_list = ["out_al"]
        for var in var_list:
            self.plot_moss_cla_x_mossbiomass(
                moss_dat[var],
                save_figs,
                plot_dir,
            )

    @staticmethod
    def plot_moss_cla_x_mossbiomass(
        data: xr.Dataset, save_fig: bool, plot_dir: str = None
    ):
        """Plot a variable with outputs dimensioned: cumulative leaf area x moss biomass

        Args:
            data (xarray DataArray): the data array of the variable to plot
            save_fig (bool): whether or not to write out plot
            plot_dir (str): if saving figure, where to write to
        """
        data_frame = pd.DataFrame(
            {
                "cla": np.tile(data.cla, len(data.moss_biomass)),
                "moss_biomass": np.repeat(data.moss_biomass, len(data.cla)),
                data.name: data.values.flatten(),
            }
        )

        max_cla = data_frame["cla"].max()
        max_var = round_up(data_frame[data.name].max())

        blank_plot(max_cla, 0.0, max_var, 0.0, draw_horizontal_lines=True)

        moss_biomasses = np.unique(data_frame.moss_biomass.values)
        colors = get_color_palette(len(moss_biomasses))
        for rank, moss_biomass in enumerate(moss_biomasses):
            dat = data_frame[data_frame.moss_biomass == moss_biomass]
            plt.plot(
                dat.cla.values,
                dat[data.name].values,
                lw=2,
                color=colors[rank],
                label=moss_biomass,
            )

        plt.xlabel("Cumulative leaf area in plot (m2)", fontsize=11)
        varname = data.attrs["long_name"]
        units = data.attrs["units"]
        plt.ylabel(f"{varname} ({units})", fontsize=11)
        plt.title(f"Simulated {varname}", fontsize=11)
        plt.legend(loc="upper left", title="Moss biomass (kg/m2 plot)")

        if save_fig:
            fig_name = os.path.join(plot_dir, f"allometry_plot_{data.name}.png")
            plt.savefig(fig_name)
