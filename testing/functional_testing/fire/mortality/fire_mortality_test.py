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

class FireMortTest(FunctionalTest):
    """Fire mortality test class"""

    name = "fire_mortality"

    def __init__(self, test_dict):
        super().__init__(
            FireMortTest.name,
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
