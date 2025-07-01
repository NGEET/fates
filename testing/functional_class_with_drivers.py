import os
from functional_class import FunctionalTest


class FunctionalTestWithDrivers(FunctionalTest):
    """Class for running FATES functional tests with driver files"""

    def __init__(self, datm_file: str, *args):

        # Check that datm exists and save its absolute path
        self.datm_file = os.path.abspath(datm_file)
        if not os.path.exists(self.datm_file):
            raise FileNotFoundError(f"datm_file not found: '{self.datm_file}'")

        super().__init__(*args)
