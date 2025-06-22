import os
import urllib.request
from functional_class import FunctionalTest


class FunctionalTestWithDrivers(FunctionalTest):
    """Class for running FATES functional tests with driver files"""

    def __init__(self, datm_file: str, datm_file_url: str, *args):

        # Check that datm exists and save its absolute path
        self.datm_file = os.path.abspath(datm_file)
        if not os.path.exists(self.datm_file):
            if not datm_file_url:
                raise FileNotFoundError(f"datm_file not found: '{self.datm_file}'")
            datm_file_dir = os.path.dirname(self.datm_file)
            if not os.path.isdir(datm_file_dir):
                os.makedirs(datm_file_dir)
            print(f"Downloading datm_file from {datm_file_url}")
            urllib.request.urlretrieve(datm_file_url, self.datm_file)

        super().__init__(*args)
