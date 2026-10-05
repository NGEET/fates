"""
Concrete class for running a LeafLevelPhoto functional tests for FATES.
"""
import os
import xarray as xr
import matplotlib.pyplot as plt
from framework.functional_test import FunctionalTest


class LeafLevelPhoto(FunctionalTest):
    """LeafLevelPhoto test class
    """
    name = "leaf_level_photo"

    # matches FatesConstantsMod's t_water_freeze_k_1atm, used only to convert
    # the temp sweep's K output to degC for plotting
    _T_FREEZE_K = 273.15

    # matches FatesTestEnvironmentMod's sea_level_press [Pa]
    _CAN_PRESS_PA = 101325.0

    # (swept dimension, x-axis label, plot-title fragment) for each of the
    # sweeps
    _SWEEPS = [
        ("par", "Absorbed PAR ($\\mu$mol m$^{-2}$ s$^{-1}$)", "PAR"),
        ("co2", "CO$_2$ (ppm)", "CO$_2$"),
        ("vpd", "Vapor pressure deficit (kPa)", "VPD"),
        ("temp", "Leaf temperature ($^{\\circ}$C)", "leaf temperature"),
        ("soilfrac", "Soil water content (fraction of saturation)", "soil water content"),
    ]

    def plot_output(self, run_dir: str, save_figs: bool, plot_dir: str):
        """Plots - update this to plot your output

        Args:
            run_dir (str): run directory
            out_file (str): output file name
            save_figs (bool): whether or not to save the figures
            plot_dir (str): plot directory to save the figures to
        """
        data = xr.open_dataset(os.path.join(run_dir, self.out_file))
        for dim, xlabel, title in self._SWEEPS:
            self.plot_sweep(data, dim, xlabel, title, save_figs, plot_dir)


    @staticmethod
    def _style_axis(axis):
        """Applies the shared minimalist axis styling used across these plots

        Args:
            axis (matplotlib.axes.Axes): axis to style
        """
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
        axis.tick_params(bottom=False, left=False)
        axis.set_axisbelow(True)
        axis.grid(axis="y", lw=0.5, alpha=0.3, linestyle="--")


    @classmethod
    def plot_sweep(cls, data: xr.Dataset, dim: str, xlabel: str, title: str,
                    save_fig: bool, plot_dir: str = None):
        """Plots gross/net photosynthesis, stomatal conductance, and
        intracellular CO2 against one swept variable, for the single PFT
        the test was run with

        Args:
            data (xarray Dataset): the leaf-level photosynthesis dataset
            dim (str): swept dimension name (e.g. "par" for
                anet_bypar/agross_bypar/gs_bypar/ci_bypar)
            xlabel (str): x-axis label
            title (str): title fragment, e.g. "PAR" -> "... vs. PAR"
            save_fig (bool): whether or not to write out the figure
            plot_dir (str): if saving figure, where to write to
        """
        x = data[dim].values
        if dim == "temp":
            x = x - cls._T_FREEZE_K
        elif dim == "vpd":
            x = x / 1000.0  # Pa -> kPa
        elif dim == "co2":
            x = x / cls._CAN_PRESS_PA * 1.0e6  # Pa -> ppm

        fig, axes = plt.subplots(2, 2, figsize=(11, 8))
        panels = [
            (axes[0, 0], f"agross_by{dim}", "Gross photosynthesis\n($\\mu$molC m$^{-2}$ s$^{-1}$)"),
            (axes[0, 1], f"anet_by{dim}", "Net photosynthesis\n($\\mu$molC m$^{-2}$ s$^{-1}$)"),
            (axes[1, 0], f"gs_by{dim}", "Stomatal conductance\n($\\mu$mol H$_2$O m$^{-2}$ s$^{-1}$)"),
            (axes[1, 1], f"ci_by{dim}", "Intracellular CO$_2$\n(ppm)"),
        ]

        for axis, varname, ylabel in panels:
            values = data[varname].values  # (n_sweep,)
            if varname.startswith("ci_by"):
                values = values / cls._CAN_PRESS_PA * 1.0e6  # Pa -> ppm
            axis.plot(x, values, lw=1.2, color="tab:blue")
            cls._style_axis(axis)
            axis.set_xlabel(xlabel, fontsize=10)
            axis.set_ylabel(ylabel, fontsize=10)

        fig.suptitle(f"Leaf-level photosynthesis vs. {title}", fontsize=12)

        if save_fig:
            fig.savefig(os.path.join(plot_dir, f"leaf_level_photo_{dim}.png"))
