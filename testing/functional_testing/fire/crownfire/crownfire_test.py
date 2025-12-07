"""
Concrete class for running the crownfire functional test for FATES.
"""
import os
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt

from functional_class import FunctionalTest
from utils_plotting import blank_plot

COLORS = ["#793922", "#99291F", "#CC9728", "#6B8939", "#2C778A", "#2C378A"]
MPERMIN_TO_KMPERHOUR = 0.06


class CrownFireTest(FunctionalTest):
    """CrownFire test class"""

    name = "crownfire"

    def __init__(self, test_dict):
        super().__init__(
            CrownFireTest.name,
            test_dict["test_dir"],
            test_dict["test_exe"],
            test_dict["out_file"],
            test_dict["use_param_file"],
            test_dict["other_args"],
        )
        self.plot = True

    def plot_output(self, run_dir: str, save_figs: bool, plot_dir: str):
        """Plot output associated with crown fire tests

        Args:
            run_dir (str): run directory
            out_file (str): output file
            save_figs (bool): whether or not to save the figures
            plot_dir (str): plot directory
        """

        # read in crown fire data
        cfire_dat = xr.open_dataset(os.path.join(run_dir, self.out_file))

        self.plot_passcrwn_fi(cfire_dat, save_figs, plot_dir)
        self.plot_lfmc(cfire_dat, save_figs, plot_dir)
        self.plot_rosact_fm10(cfire_dat, save_figs, plot_dir)
        self.plot_ci_fm10(cfire_dat, save_figs, plot_dir)

    @staticmethod
    def plot_passcrwn_fi(data: xr.Dataset, save_fig: bool, plot_dir: str = None):
        """Plot min fire intensity required to ignite crown fuel

        Args:
            data (xarray DataSet): the data set
            save_fig (bool): whether or not to write out plot
            plot_dir (str): if saving figure, where to write to
        """
        data_frame = pd.DataFrame(
            {
                "CBH": np.tile(data.CBH, len(data.CWC)),
                "CWC": np.repeat(data.CWC, len(data.CBH)),
                "passive_crown_fi": data.passive_crown_fi.values.flatten(),
            }
        )


        max_CBH = data_frame["CBH"].max()

        plt.figure(figsize=(7, 5))
        axis = plt.subplot(111)
        axis.spines["top"].set_visible(False)
        axis.spines["bottom"].set_visible(False)
        axis.spines["right"].set_visible(False)
        axis.spines["left"].set_visible(False)

        axis.get_xaxis().tick_bottom()
        axis.get_yaxis().tick_left()

        plt.xlim(0.0, max_CBH)

        plt.yticks(fontsize=10)
        plt.xticks(fontsize=10)
        plt.grid(True)

        CWC = np.unique(data_frame.CWC.values)

        for i, cwc in enumerate(CWC):
            dat = data_frame[data_frame.CWC == cwc]
            plt.plot(
                dat.CBH.values,
                dat["passive_crown_fi"].values,
                lw=2,
                color=COLORS[i],
                label=cwc,
            )

        plt.xlabel("Canopy base height (m)", fontsize=11)
        plt.ylabel("Threshold fire intensity to ignite crown fuels (kW m$^{-1}$)", fontsize=11)
        plt.legend(loc="upper left", title="Canopy water content (%)")

        if save_fig:
            fig_name = os.path.join(plot_dir, "passcrown_fi_plot.png")
            plt.savefig(fig_name)

    @staticmethod
    def plot_lfmc(data: xr.Dataset, save_fig: bool, plot_dir: str = None):
        """Plot live fuel moisture content

        Args:
            data (xarray DataSet): the data set
            save_fig (bool): whether or not to write out plot
            plot_dir (str): if saving figure, where to write to
        """
        data_frame = pd.DataFrame(
            {
                "smp": np.tile(data.smp, len(data.smp_alpha)),
                "smp_alpha": np.repeat(data.smp_alpha, len(data.smp)),
                "LFMC": data.LFMC.values.flatten(),
            }
        )

        min_smp = float(data_frame["smp"].min())
        max_lfmc = 150.0

        blank_plot(0.0, -10.0, 150.0, 0.0, draw_horizontal_lines=True)

        smp_alpha_vals = np.unique(data_frame.smp_alpha.values)

        for i, alpha in enumerate(smp_alpha_vals):
            dat = data_frame[data_frame["smp_alpha"] == alpha]
            dat = dat.sort_values("smp")
           
            plt.plot(
                dat.smp.values,
                dat["LFMC"].values,
                lw=2,
                color=COLORS[i],
                label=alpha,
            )

        plt.xlabel("Soil matric potential (MPa)", fontsize=11)
        plt.ylabel("Live fuel moisture (%)", fontsize=11)
        plt.legend(loc="upper right", title="Soil matric potential coeff")
        plt.xlim(0.0,-10.0)

        if save_fig:
            fig_name = os.path.join(plot_dir, "lfmc_plot.png")
            plt.savefig(fig_name)

    @staticmethod
    def plot_rosact_fm10(data: xr.Dataset, save_fig: bool, plot_dir: str = None):
        """Plot active crown fire spread rate

        Args:
            data (xarray DataSet): the data set
            save_fig (bool): whether or not to write out plot
            plot_dir (str): if saving figure, where to write to
        """
        data_frame = pd.DataFrame(
            {   
                "wind": np.tile(data.wind, len(data.dratio)),
                "drying_ratio": np.repeat(data.dratio, len(data.wind)),
                "ros_active_fm10": data.ROSACT_FM10.values.flatten(),
            }
        )

        data_frame["wind_kmhr"] = data_frame.wind * MPERMIN_TO_KMPERHOUR  # match usual wind speed unit in crown fire plots

        max_wind = data_frame["wind_kmhr"].max()
        max_ros_active = data_frame["ros_active_fm10"].max()

        blank_plot(max_wind, 0.0, max_ros_active, 0.0, draw_horizontal_lines=True)

        dr_vals = np.unique(data_frame.drying_ratio.values)
        for i, dr in enumerate(dr_vals):
            dat = data_frame[data_frame["drying_ratio"] == dr]

            plt.plot(
            dat.wind_kmhr.values,
            dat["ros_active_fm10"].values,
            lw=2,
            color=COLORS[i],
            label = dr,
            )

        plt.xlabel("Wind speed (km hr$^{-1}$)", fontsize=11)
        plt.ylabel("Active crown fire ROS (m min$^{-1}$)", fontsize=11)
        plt.legend(loc="upper right", title="Drying ratio")

        if save_fig:
            fig_name = os.path.join(plot_dir, "ros_active_plot.png")
            plt.savefig(fig_name)

    @staticmethod
    def plot_ci_fm10(data: xr.Dataset, save_fig: bool, plot_dir: str = None):
        """Plot crowning index

        Args:
            data (xarray DataSet): the data set
            save_fig (bool): whether or not to write out plot
            plot_dir (str): if saving figure, where to write to
        """
        data_frame = pd.DataFrame(
            {
                "CBD": np.tile(data.CBD, len(data.dratio)),
                "drying_ratio": np.repeat(data.dratio, len(data.CBD)),
                "CI_FM10": data.CI_FM10.values.flatten(),
            }
        )

        data_frame["CI_kmhr"] = data_frame.CI_FM10* MPERMIN_TO_KMPERHOUR

        max_CBD = data_frame["CBD"].max()
        max_CI = data_frame["CI_kmhr"].max()

        blank_plot(max_CBD, 0.0, max_CI, 0.0, draw_horizontal_lines=True)

        dr_vals = np.unique(data_frame.drying_ratio.values)

        for i, dr in enumerate(dr_vals):
            dat = data_frame[data_frame["drying_ratio"] == dr]

            plt.plot(
            dat.CBD.values,
            dat["CI_kmhr"].values,
            lw=2,
            color=COLORS[i],
            label = dr
            )

        plt.xlabel("Canopy bulk density (kg m$^{-3}$)", fontsize=11)
        plt.ylabel("Crowning index (km hr$^{-1}$)", fontsize=11)
        plt.legend(loc = "upper right", title="Drying ratio")

        if save_fig:
            fig_name = os.path.join(plot_dir, "ci_plot.png")
            plt.savefig(fig_name)