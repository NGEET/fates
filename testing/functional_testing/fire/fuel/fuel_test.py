"""
Concrete class for running the fuel functional test for FATES.
"""
import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from functional_class_with_drivers import FunctionalTestWithDrivers


class FuelTest(FunctionalTestWithDrivers):
    """Fuel test class"""

    name = "fuel"

    def __init__(self, test_dict):
        super().__init__(
            test_dict["datm_file"],
            FuelTest.name,
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

        # define all timeseries vars
        # concise: define all time-series plots in one list
        ts_specs = [
        ("NI", "Nesterov Index", None, "Nesterov_plot.png"),
        ("fuel_moisture_dead", "Dead fuel Moisture", "fuel_model", "fuel_moisture_dead_plot.png"),
        ("fuel_moisture_live", "Live fuel Moisture", "fuel_model", "fuel_moisture_live_plot.png"),
        ]

        fuel_dat = xr.open_dataset(os.path.join(run_dir, self.out_file))

        for var, ylabel, hue, filename in ts_specs:
            if var in fuel_dat:
                self.plot_timeseries(fuel_dat, var, ylabel, hue, save_figs, plot_dir, filename)

        self.plot_barchart(
            fuel_dat,
            "fuel_loading",
            "Fuel loading",
            "kgC m$^{-2}$",
            save_figs,
            plot_dir,
        )
        self.plot_barchart(
            fuel_dat,
            "weighted_loading_dead",
            "Weighted dead fuel loading",
            "kgC m$^{-2}$",
            save_figs,
            plot_dir,
            by_litter_type = False,
        )
        self.plot_barchart(
            fuel_dat,
            "weighted_loading_live",
            "Weighted live fuel loading",
            "kgC m$^{-2}$",
            save_figs,
            plot_dir,
            by_litter_type=False,
        )
        self.plot_barchart(
            fuel_dat,
            "frac_loading",
            "Fractional fuel loading",
            "0-1",
            save_figs,
            plot_dir,
        )
        self.plot_barchart(
            fuel_dat,
            "bulk_density",
            "Fuel bulk density",
            "kg m$^{-3}$",
            save_figs,
            plot_dir,
            by_litter_type=False,
        )
        self.plot_barchart(
            fuel_dat,
            "SAV",
            "Fuel surface area to volume ratio",
            "cm$^{-1}$",
            save_figs,
            plot_dir,
            by_litter_type=False,
        )

    @staticmethod
    def plot_barchart(
        fuel_dat: xr.Dataset,
        var: str,
        varname: str,
        units: str,
        save_figs: bool,
        plot_dir: bool,
        by_litter_type: bool = True,
    ):
        """Plots fuel data output as a bar chart

        Args:
            fuel_dat (xr.Dataset): fuel data output
            var (str): variable to plot
            varname (str): variable name for x axis
            units (str): units description
            save_figs (bool): whether or not to save figure
            plot_dir (bool): where to save figure
            by_litter_type (bool, optional): whether the bar chart is by litter type. Defaults to True.
        """

        litter_classes = [
            "twigs",
            "small branches",
            "large branches",
            "trunks",
            "dead leaves",
            "live grass",
        ]
        colors = [
            "darksalmon",
            "peru",
            "saddlebrown",
            "black",
            "moccasin",
            "yellowgreen",
        ]
        fuel_models = [str(f) for f in fuel_dat.fuel_model.values]

        if by_litter_type:
            data_dict = {
                lc: fuel_dat.isel(litter_class=i)[var].values
                for i, lc in enumerate(litter_classes)
            }
        else:
            data_dict = fuel_dat[var].values

        _, ax = plt.subplots()
        if by_litter_type:
            bottom = np.zeros(len(fuel_models))
            for i, (litter_class, dat) in enumerate(data_dict.items()):
                ax.bar(
                    fuel_models,
                    dat,
                    0.5,
                    label=litter_class,
                    bottom=bottom,
                    color=colors[i],
                )
                bottom += dat
            plt.legend(loc="center left", bbox_to_anchor=(1, 0.5))
        else:
            ax.bar(fuel_models, data_dict, color="darkcyan")

        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.75, box.height])
        plt.ylabel(f"{varname} ({units})", fontsize=11)
        plt.xticks(rotation=90)
        plt.xlabel("Fuel Model")

        if save_figs:
            fig_name = os.path.join(plot_dir, f"{varname}_plot.png")
            plt.savefig(fig_name)

    @staticmethod
    def plot_timeseries(
        fuel_dat: xr.Dataset, 
        var: str,
        ylabel: str, 
        hue: str | None,
        save_figs: bool, 
        plot_dir: str,
        filename: str,
        ):
        """Plot output as time series

        Args:
            fuel_dat (Xarray Dataset): output fuel data
            var (str): variable to plot
            ylabel (str): label for y axis 
            hue (str): grouping variable (optional)
            save_figs (bool): whether or not to save the figures
            plot_dir (str): plot directory
            filename (str): file name for saving figures
        """
        if var not in fuel_dat:
            return
        
        fig, ax = plt.subplots()
        if hue is None:
            fuel_dat[var].plot(ax=ax)
        else:
            fuel_dat[var].plot(hue = hue, ax=ax)

        ax.set_xlabel("Time", fontsize=11)
        ax.set_ylabel(ylabel, fontsize=11)
        plt.tight_layout()
    
        if save_figs:
            fig_name = os.path.join(plot_dir, filename)
            plt.savefig(fig_name, dpi=150, bbox_inches="tight")
        plt.close(fig)

