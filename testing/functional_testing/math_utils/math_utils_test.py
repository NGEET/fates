"""
Concrete class for running the quadtratic functional tests for FATES.
"""
import os
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from utils import get_color_palette
from functional_class import FunctionalTest


class QuadraticTest(FunctionalTest):
    """Quadratic test class
    """

    name = "quadratic"

    def __init__(self, test_dict):
        super().__init__(
            QuadraticTest.name,
            test_dict["test_dir"],
            test_dict["test_exe"],
            test_dict["out_file"],
            test_dict["use_param_file"],
            test_dict["other_args"],
        )
        self.plot = True

    def plot_output(self, run_dir: str, save_figs: bool, plot_dir: str):
        """Reads in and plots quadratic formula test output

        Args:
            run_dir (str): run directory
            out_file (str): output file
            save_figs (bool): whether or not to save the figures
            plot_dir (str): plot directory
        """

        # read in quadratic data
        quadratic_dat = xr.open_dataset(os.path.join(run_dir, self.out_file))

        # plot output
        self.plot_quad_and_roots(
            quadratic_dat.a.values,
            quadratic_dat.b.values,
            quadratic_dat.c.values,
            quadratic_dat.root1.values,
            quadratic_dat.root2.values,
        )
        if save_figs:
            fig_name = os.path.join(plot_dir, "quadratic_test.png")
            plt.savefig(fig_name)

    @staticmethod
    def plot_quad_and_roots(a_coeff, b_coeff, c_coeff, root1, root2):
        """Plots a set of quadratic formulas (ax**2 + bx + c) and their two roots

        Args:
            a_coeff (float array): set of a coefficients
            b_coeff (float array): set of b coefficients
            c_coeff (float array): set of b coefficients
            root1 (float array): set of first real roots
            root2 (float array): set of second real roots
        """
        num_equations = len(a_coeff)

        plt.figure(figsize=(7, 5))
        x_vals = np.linspace(-10.0, 10.0, num=20)

        colors = get_color_palette(num_equations)
        for i in range(num_equations):
            y_vals = a_coeff[i] * x_vals**2 + b_coeff[i] * x_vals + c_coeff[i]
            plt.plot(x_vals, y_vals, lw=2, color=colors[i])
            plt.scatter(root1[i], root2[i], color=colors[i], s=50)
            plt.axhline(y=0.0, color="k", linestyle="dotted")
