"""
Concrete class for running a GreatCircle functional tests for FATES.
"""
import os
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from framework.functional_test import FunctionalTest


class GreatCircle(FunctionalTest):
    """GreatCircle test class
    """
    name = "great_circle"

    def plot_output(self, run_dir: str, save_figs: bool, plot_dir: str):
        """Plots - update this to plot your output

        Args:
            run_dir (str): run directory
            out_file (str): output file name
            save_figs (bool): whether or not to save the figures
            plot_dir (str): plot directory to save the figures to
        """
        pass

