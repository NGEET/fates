from abc import ABC, abstractmethod
from utils import str_to_bool, str_to_list

class Heuristic(ABC):
    """Class for running FATES heuristics"""

    def __init__(self, name, test_dir, test_exe, out_file, use_param_file, other_args):
        super().__init__(name, test_dir, test_exe, 'functional')
        self.out_file = out_file
        self.use_param_file = str_to_bool(use_param_file)
        self.other_args = str_to_list(other_args)
    
    @abstractmethod
    def plot_output(self, run_dir:str, out_file:str, save_figs:bool, plot_dir:str):
        """plot output for this test

        Args:
            run_dir (str): path to run directory
            out_file (str): name of output file to read in
            save_figs (bool): whether or not to save figs to png
            plot_dir (str): where to save figs
        """