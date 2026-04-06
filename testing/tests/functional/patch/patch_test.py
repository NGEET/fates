"""
Concrete class for running the allometry functional tests for FATES.
"""
from framework.functional_test import FunctionalTest


class PatchTest(FunctionalTest):
    """Patch test class
    """
    name = 'patch'

    def plot_output(self, run_dir: str, save_figs: bool, plot_dir: str):
        """Plots all allometry plots

        Args:
            run_dir (str): run directory
            out_file (str): output file name
            save_figs (bool): whether or not to save the figures
            plot_dir (str): plot directory to save the figures to
        """
        pass
