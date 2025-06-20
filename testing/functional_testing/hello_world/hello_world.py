"""
Concrete class for running a HelloWorld functional tests for FATES.
"""
import os
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from functional_class import FunctionalTest


class HelloWorld(FunctionalTest):
    """Quadratic test class
    """

    name = "allometry"

    def __init__(self, test_dict):
        super().__init__(
            HelloWorld.name,
            test_dict["test_dir"],
            test_dict["test_exe"],
            test_dict["out_file"],
            test_dict["use_param_file"],
            test_dict["other_args"],
        )
        self.plot = True

    def plot_output(self, run_dir: str, save_figs: bool, plot_dir: str):
        """Plots - update this to plot your output

        Args:
            run_dir (str): run directory
            out_file (str): output file name
            save_figs (bool): whether or not to save the figures
            plot_dir (str): plot directory to save the figures to
        """

        pass

