from abc import ABC, abstractmethod
from utils import str_to_bool, str_to_list

class FunctionalTest(ABC):
    """Class for running FATES functional tests"""

    def __init__(self, name:str, test_dir:str, test_exe:str, out_file:str,
                 use_param_file:str, other_args:str):
        self.name = name
        self.test_dir = test_dir
        self.test_exe = test_exe
        self.out_file = out_file
        self.use_param_file = str_to_bool(use_param_file)
        self.other_args = str_to_list(other_args)
        self.plot = False
        
    @abstractmethod
    def plot_output(self, run_dir:str, save_figs:bool, plot_dir:str):
        pass
    