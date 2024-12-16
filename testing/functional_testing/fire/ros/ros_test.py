"""
Concrete class for running the ros functional test for FATES.
"""
import os
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
from functional_class import FunctionalTest
from utils import blank_plot

COLORS = ["#793922", "#6B8939", "#99291F", "#CC9728", "#2C778A"]
CM_TO_FT = 30.48
KJKG_TO_BTULB = 0.947817 / 2.20462


class ROSTest(FunctionalTest):
    """ROS test class"""

    name = "ros"

    def __init__(self, test_dict):
        super().__init__(
            ROSTest.name,
            test_dict["test_dir"],
            test_dict["test_exe"],
            test_dict["out_file"],
            test_dict["use_param_file"],
            test_dict["other_args"],
        )
        self.plot = True

    def plot_output(self, run_dir: str, save_figs: bool, plot_dir: str):
        """Plot output associated with fuel tests

        Args:
            run_dir (str): run directory
            out_file (str): output file
            save_figs (bool): whether or not to save the figures
            plot_dir (str): plot directory
        """

        # read in ros data
        ros_dat = xr.open_dataset(os.path.join(run_dir, self.out_file))

        self.plot_prop_flux(ros_dat, save_figs, plot_dir)
        self.plot_reaction_vel(ros_dat, save_figs, plot_dir)
        self.plot_qig(ros_dat, save_figs, plot_dir)
        self.plot_eps(ros_dat, save_figs, plot_dir)

    @staticmethod
    def plot_prop_flux(data: xr.Dataset, save_fig: bool, plot_dir: str = None):
        """Plot propagating flux

        Args:
            data (xarray DataSet): the data set
            save_fig (bool): whether or not to write out plot
            plot_dir (str): if saving figure, where to write to
        """
        data_frame = pd.DataFrame(
            {
                "SAV": np.tile(data.SAV, len(data.packing_ratio)),
                "packing_ratio": np.repeat(data.packing_ratio, len(data.SAV)),
                "prop_flux": data.prop_flux.values.flatten(),
            }
        )

        data_frame["SAV_ft"] = (
            data_frame.SAV * CM_TO_FT
        )  # covert to ft to compare with original Rothermel equations

        max_SAV = data_frame["SAV_ft"].max()
        max_prop_flux = 0.14

        plt.figure(figsize=(7, 5))
        axis = plt.subplot(111)
        axis.spines["top"].set_visible(False)
        axis.spines["bottom"].set_visible(False)
        axis.spines["right"].set_visible(False)
        axis.spines["left"].set_visible(False)

        axis.get_xaxis().tick_bottom()
        axis.get_yaxis().tick_left()

        plt.xlim(0.0, max_SAV)
        plt.ylim(0.0, max_prop_flux)

        plt.yticks(fontsize=10)
        plt.xticks(fontsize=10)
        plt.grid(True)

        packing_ratio = np.unique(data_frame.packing_ratio.values)

        for i, beta in enumerate(packing_ratio):
            dat = data_frame[data_frame.packing_ratio == beta]
            plt.plot(
                dat.SAV_ft.values,
                dat["prop_flux"].values,
                lw=2,
                color=COLORS[i],
                label=beta,
            )

        plt.xlabel("Surface-area-to-volume ratio (ft$^{-1}$)", fontsize=11)
        plt.ylabel("Propagating flux ratio", fontsize=11)
        plt.legend(loc="upper left", title="Packing ratio")

        if save_fig:
            fig_name = os.path.join(plot_dir, "prop_flux_plot.png")
            plt.savefig(fig_name)

    @staticmethod
    def plot_reaction_vel(data: xr.Dataset, save_fig: bool, plot_dir: str = None):
        """Plot reaction velocity

        Args:
            data (xarray DataSet): the data set
            save_fig (bool): whether or not to write out plot
            plot_dir (str): if saving figure, where to write to
        """
        data_frame = pd.DataFrame(
            {
                "beta_ratio": np.tile(data.beta_ratio, len(data.SAV_ind)),
                "SAV": np.repeat(data.SAV_ind, len(data.beta_ratio)),
                "reaction_vel": data.reaction_velocity.values.flatten(),
            }
        )

        data_frame["SAV_ft"] = data_frame.SAV * CM_TO_FT

        max_beta = data_frame["beta_ratio"].max()
        max_reaction_vel = 18

        blank_plot(max_beta, 0.0, max_reaction_vel, 0.0, draw_horizontal_lines=True)

        SAV_vals = np.unique(data_frame.SAV_ft.values)
        colors = COLORS
        colors.reverse()

        for i, sav in enumerate(SAV_vals):
            dat = data_frame[data_frame.SAV_ft == sav]
            plt.plot(
                dat.beta_ratio.values,
                dat["reaction_vel"].values,
                lw=2,
                color=colors[i],
                label=sav,
            )

        plt.xlabel("Relative packing ratio", fontsize=11)
        plt.ylabel("Reaction velocity (min$^{-1}$)", fontsize=11)
        plt.legend(loc="upper right", title="Surface-area-to-volume ratio (ft$^{-1}$)")

        if save_fig:
            fig_name = os.path.join(plot_dir, "reaction_vel_plot.png")
            plt.savefig(fig_name)

    @staticmethod
    def plot_qig(data: xr.Dataset, save_fig: bool, plot_dir: str = None):
        """Plot heat of preignition

        Args:
            data (xarray DataSet): the data set
            save_fig (bool): whether or not to write out plot
            plot_dir (str): if saving figure, where to write to
        """
        data_frame = pd.DataFrame(
            {
                "fuel_moisture": data.fuel_moisture,
                "q_ig": data.q_ig.values.flatten(),
            }
        )

        data_frame["fuel_moisture_perc"] = data_frame.fuel_moisture * 100.0
        data_frame["q_ig_btu"] = (
            data_frame.q_ig * KJKG_TO_BTULB
        )  # match Rothermel graph units

        max_moist = 200.0
        max_qig = 2500.0

        blank_plot(max_moist, 0.0, max_qig, 0.0, draw_horizontal_lines=True)

        plt.plot(
            data_frame.fuel_moisture_perc.values,
            data_frame["q_ig_btu"].values,
            lw=2,
            color="k",
        )

        plt.xlabel("Fuel moisture (%)", fontsize=11)
        plt.ylabel("Heat of Preignition (Btu lb$^{-1}$)", fontsize=11)

        if save_fig:
            fig_name = os.path.join(plot_dir, "qig_plot.png")
            plt.savefig(fig_name)

    @staticmethod
    def plot_eps(data: xr.Dataset, save_fig: bool, plot_dir: str = None):
        """Plot effective heating number

        Args:
            data (xarray DataSet): the data set
            save_fig (bool): whether or not to write out plot
            plot_dir (str): if saving figure, where to write to
        """
        data_frame = pd.DataFrame(
            {
                "SAV": data.SAV,
                "eps": data.eps.values.flatten(),
            }
        )

        data_frame["SAV_ft"] = data_frame.SAV * CM_TO_FT

        max_SAV = 3500.0
        max_eps = 1.0

        blank_plot(max_SAV, 0.0, max_eps, 0.0, draw_horizontal_lines=True)

        plt.plot(
            data_frame.SAV_ft.values,
            data_frame["eps"].values,
            lw=2,
            color="k",
        )

        plt.xlabel("Surface-area-to-volume ratio (ft$^{-1}$)", fontsize=11)
        plt.ylabel("Effective heating number", fontsize=11)

        if save_fig:
            fig_name = os.path.join(plot_dir, "eps_plot.png")
            plt.savefig(fig_name)
