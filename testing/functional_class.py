import os
import urllib.request
from abc import ABC, abstractmethod
from utils import str_to_bool, str_to_list

class FunctionalTest(ABC):
    """Class for running FATES functional tests"""

    def __init__(self, name:str, test_dir:str, test_exe:str, out_file:str,
                 use_param_file:str, datm_file:str, datm_file_url:str, other_args:str):
        self.name = name
        self.test_dir = test_dir
        self.test_exe = test_exe
        self.out_file = out_file
        self.use_param_file = str_to_bool(use_param_file)
        self.datm_file = None
        self.other_args = str_to_list(other_args)
        self.plot = False

        # Check that datm exists and save its absolute path
        if datm_file:
            self.datm_file = os.path.abspath(datm_file)
            if not os.path.exists(self.datm_file):
                if not datm_file_url:
                    raise FileNotFoundError(f"datm_file not found: '{self.datm_file}'")
                datm_file_dir = os.path.dirname(self.datm_file)
                if not os.path.isdir(datm_file_dir):
                    os.makedirs(datm_file_dir)
                print(f"Downloading datm_file from {datm_file_url}")
                urllib.request.urlretrieve(datm_file_url, self.datm_file)

    @abstractmethod
    def plot_output(self, run_dir:str, save_figs:bool, plot_dir:str):
        pass
    