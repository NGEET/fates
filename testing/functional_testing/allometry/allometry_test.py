"""
Concrete class for running the allometry functional tests for FATES.
"""
import os
import xarray as xr
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from utils import round_up, get_color_palette, blank_plot
from functional_class import FunctionalTest


class AllometryTest(FunctionalTest):
    """Quadratic test class
    """

    name = "allometry"

    def __init__(self, test_dict):
        super().__init__(
            AllometryTest.name,
            test_dict["test_dir"],
            test_dict["test_exe"],
            test_dict["out_file"],
            test_dict["use_param_file"],
            test_dict["other_args"],
        )
        self.plot = True

    def plot_output(self, run_dir: str, save_figs: bool, plot_dir: str):
        """Plots all allometry plots

        Args:
            run_dir (str): run directory
            out_file (str): output file name
            save_figs (bool): whether or not to save the figures
            plot_dir (str): plot directory to save the figures to
        """

        # read in allometry data
        allometry_dat = xr.open_dataset(os.path.join(run_dir, self.out_file))

        plot_dict = {
            "height": {
                "varname": "height",
                "units": "m",
            },
            "bagw": {
                "varname": "aboveground biomass",
                "units": "kgC",
            },
            "blmax": {
                "varname": "maximum leaf biomass",
                "units": "kgC",
            },
            "crown_area": {
                "varname": "crown area",
                "units": "m$^2$",
            },
            "sapwood_area": {
                "varname": "sapwood area",
                "units": "m$^2$",
            },
            "bsap": {
                "varname": "sapwood biomass",
                "units": "kgC",
            },
            "bbgw": {
                "varname": "belowground biomass",
                "units": "kgC",
            },
            "fineroot_biomass": {
                "varname": "fineroot biomass",
                "units": "kgC",
            },
            "bstore": {
                "varname": "storage biomass",
                "units": "kgC",
            },
            "bdead": {
                "varname": "deadwood biomass",
                "units": "kgC",
            },
            "total_biomass_parts": {
                "varname": "total biomass (calculated from parts)",
                "units": "kgC",
            },
            "total_biomass_tissues": {
                "varname": "total biomass (calculated from tissues)",
                "units": "kgC",
            },
        }
        for plot, attributes in plot_dict.items():
            self.plot_allometry_var(
                allometry_dat[plot],
                attributes["varname"],
                attributes["units"],
                save_figs,
                plot_dir,
            )

        self.plot_total_biomass(allometry_dat, save_figs, plot_dir)

    @staticmethod
    def plot_allometry_var(
        data: xr.Dataset, varname: str, units: str, save_fig: bool, plot_dir: str = None
    ):
        """Plot an allometry variable

        Args:
            data (xarray DataArray): the data array of the variable to plot
            var (str): variable name (for data structure)
            varname (str): variable name for plot labels
            units (str): variable units for plot labels
            save_fig (bool): whether or not to write out plot
            plot_dir (str): if saving figure, where to write to
        """
        data_frame = pd.DataFrame(
            {
                "dbh": np.tile(data.dbh, len(data.pft)),
                "pft": np.repeat(data.pft, len(data.dbh)),
                data.name: data.values.flatten(),
            }
        )

        max_dbh = data_frame["dbh"].max()
        max_var = round_up(data_frame[data.name].max())

        blank_plot(max_dbh, 0.0, max_var, 0.0, draw_horizontal_lines=True)

        pfts = np.unique(data_frame.pft.values)
        colors = get_color_palette(len(pfts))
        for rank, pft in enumerate(pfts):
            dat = data_frame[data_frame.pft == pft]
            plt.plot(
                dat.dbh.values,
                dat[data.name].values,
                lw=2,
                color=colors[rank],
                label=pft,
            )

        plt.xlabel("DBH (cm)", fontsize=11)
        plt.ylabel(f"{varname} ({units})", fontsize=11)
        plt.title(f"Simulated {varname} for input parameter file", fontsize=11)
        plt.legend(loc="upper left", title="PFT")

        if save_fig:
            fig_name = os.path.join(plot_dir, f"allometry_plot_{data.name}.png")
            plt.savefig(fig_name)

    @staticmethod
    def plot_total_biomass(data: xr.Dataset, save_fig: bool, plot_dir: str):
        """Plot two calculations of total biomass against each other

        Args:
            data (xarray DataSet): the allometry dataset
        """
        data_frame = pd.DataFrame(
            {
                "dbh": np.tile(data.dbh, len(data.pft)),
                "pft": np.repeat(data.pft, len(data.dbh)),
                "total_biomass_parts": data.total_biomass_parts.values.flatten(),
                "total_biomass_tissues": data.total_biomass_tissues.values.flatten(),
            }
        )

        max_biomass = np.maximum(
            data_frame["total_biomass_parts"].max(),
            data_frame["total_biomass_tissues"].max(),
        )

        blank_plot(max_biomass, 0.0, max_biomass, 0.0, draw_horizontal_lines=False)

        pfts = np.unique(data_frame.pft.values)
        colors = get_color_palette(len(pfts))
        for rank, pft in enumerate(pfts):
            data = data_frame[data_frame.pft == pft]
            plt.scatter(
                data.total_biomass_parts.values,
                data.total_biomass_parts.values,
                color=colors[rank],
                label=pft,
            )

        plt.xlabel("Total biomass (kgC) from parts", fontsize=11)
        plt.ylabel("Total biomass (kgC) from tissues", fontsize=11)
        plt.title("Simulated total biomass for input parameter file", fontsize=11)
        plt.legend(loc="upper left", title="PFT")

        if save_fig:
            fig_name = os.path.join(
                plot_dir, "allometry_plot_total_biomass_compare.png"
            )
            plt.savefig(fig_name)
