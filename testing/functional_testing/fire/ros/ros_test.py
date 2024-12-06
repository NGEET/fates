"""
Concrete class for running the fuel functional test for FATES.
"""
import os
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
from functional_class import FunctionalTest


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
        
        data_frame['SAV_ft'] = data_frame.SAV*30.48

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
        colors = ['#6B8939', '#99291F','#CC9728', '#2C778A']
        
        for i, beta in enumerate(packing_ratio):
            dat = data_frame[data_frame.packing_ratio == beta]
            plt.plot(
                dat.SAV_ft.values,
                dat['prop_flux'].values,
                lw=2,
                color=colors[i],
                label=beta,
            )

        plt.xlabel("Surface-area-to-volume ratio (ft$^{-1}$)", fontsize=11)
        plt.ylabel("Propagating flux ratio", fontsize=11)
        plt.legend(loc="upper left", title="Packing ratio")

        if save_fig:
            fig_name = os.path.join(plot_dir, "prop_flux_plot.png")
            plt.savefig(fig_name)
        