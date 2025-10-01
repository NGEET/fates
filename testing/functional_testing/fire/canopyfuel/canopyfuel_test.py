
"""
Concrete class for running the canopyfuel functional test for FATES.
"""
import os
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt

from functional_class import FunctionalTest

class CanopyFuelTest(FunctionalTest):
    """Canopy fuel test class"""

    name = "canopyfuel"

    def __init__(self, test_dict):
        super().__init__(
            CanopyFuelTest.name,
            test_dict["test_dir"],
            test_dict["test_exe"],
            test_dict["out_file"],
            test_dict["use_param_file"],
            test_dict["other_args"],
        )
        self.plot = True

    def plot_output(self, run_dir: str, save_figs: bool, plot_dir: str):
        """Plot output associated with canopy fuel tests

        Args:
            run_dir (str): run directory
            out_file (str): output file
            save_figs (bool): whether or not to save the figures
            plot_dir (str): plot directory
        """

        # read in canopy fuel data
        cfuel_dat = xr.open_dataset(os.path.join(run_dir, self.out_file))

        self.plot_barchart(
            cfuel_dat,
            "CBD",
            "Canopy bulk density",
            "kg m$^{-3}$",
            save_figs,
            plot_dir,
            by_fuel_model=False,
            stacked= False,
        )
        self.plot_barchart(
            cfuel_dat,
            "CBH",
            "Canopy base height",
            "m",
            save_figs,
            plot_dir,
            by_fuel_model=False,
            stacked= False,
        )
        self.plot_barchart(
            cfuel_dat,
            "FI_init",
            "Fire intensity to initiate crown fire",
            "kW m$^{-1}$",
            save_figs,
            plot_dir,
            by_fuel_model=False,
            stacked= False,
        )
        self.plot_barchart(
            cfuel_dat,
            "FI",
            "Fire intensity of surface fire",
            "kW m$^{-1}$",
            save_figs,
            plot_dir,
            by_fuel_model=True,
            stacked= False,
        )
        self.plot_barchart(
            cfuel_dat,
            "ROS_active",
            "Active crown fire ROS",
            "m min$^{-1}$",
            save_figs,
            plot_dir,
            by_fuel_model=True,
            stacked= False,
        )
        self.plot_barchart(
            cfuel_dat,
            "ROS_min",
            "Critical rate of spread",
            "m min$^{-1}$",
            save_figs,
            plot_dir,
            by_fuel_model=False,
            stacked= False,
        )
        self.plot_barchart(
            cfuel_dat,
            "ROS_final",
            "Final rate of spread",
            "m min$^{-1}$",
            save_figs,
            plot_dir,
            by_fuel_model=True,
            stacked= False,
        )
        self.plot_barchart(
            cfuel_dat,
            "FI_final",
            "Final fire intensity",
            "kW m$^{-1}$",
            save_figs,
            plot_dir,
            by_fuel_model=True,
            stacked= False,
        )

    @staticmethod
    def plot_barchart(
        cfuel_dat: xr.Dataset,
        var: str,
        varname: str,
        units: str,
        save_figs: bool,
        plot_dir: str,
        by_fuel_model: bool = False,
        stacked: bool = False,
          ):
        """ Plot canopy fuel test outputs as bar plot
        Args:
        cfuel_dat (xr.Dataset): canopy fuel test data output
        var (str): variable to plot
        varname (str): variable name for y-axis
        units (str): unit expression
        save_figs (str): dir to save figure
        plot_dir (str): which dir to save figs
        by_fuel_model (bool, optional): whether or not the bar plot is grouped by fuel model, default to True

        """
        assert var in cfuel_dat, f"{var!r} not found in dataset variables."
        da = cfuel_dat[var]

        # figure out the dimensions we care about
        has_patch = "patch_type" in da.dims
        has_fm    = "fuel_model" in da.dims

        if not has_patch:
            raise ValueError(f"{var} is missing 'patch_type' dimension; got dims {da.dims}")
    
        patch_types = [str(p) for p in cfuel_dat["patch_type"].values]


        if by_fuel_model and has_fm:
            da2 = da.transpose("patch_type", "fuel_model")
            arr = da2.to_numpy()
            fm_labels = [str(f) for f in cfuel_dat["fuel_model"].values]

            fig, ax = plt.subplots(figsize=(8, 4.6))
            x = np.arange(len(patch_types))
            n_fm = arr.shape[1]
            width = 0.8 / n_fm

            if stacked:
                bottom = np.zeros(len(patch_types))
                for j in range(n_fm):
                    ax.bar(
                        x,
                        arr[:, j],
                        width=0.8, 
                        bottom=bottom,
                        label=fm_labels[j],
                    )
                    bottom += arr[:, j]
            else:
                for j in range(n_fm):
                    offset = (j - (n_fm - 1) / 2) * width
                    ax.bar(
                        x + offset,
                        arr[:, j],
                        width=width,
                        label=fm_labels[j],
                        )
            ax.set_xticks(x)
            ax.set_xticklabels(patch_types, rotation=90)
            ax.legend(title="fuel_model", loc="center left", bbox_to_anchor=(1, 0.5))
        else:
            if has_fm:
                da_plot = da.mean(dim="fuel_model")
            else:
                da_plot = da
            
            da_plot = da_plot.transpose("patch_type", ...)
            y = da_plot.to_numpy()
            fig, ax = plt.subplots(figsize=(7.5, 4.2))
            ax.bar(patch_types, y)
            ax.set_xticklabels(patch_types, rotation=90)

        ax.set_ylabel(f"{varname} ({units})", fontsize=11)
        ax.set_xlabel("Patch Type")
        ax.set_title(varname)
        plt.tight_layout()



        if save_figs:
            fig_name = os.path.join(plot_dir, f"{varname}_plot.png")
            plt.savefig(fig_name)