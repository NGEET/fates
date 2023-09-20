import os 
import sys
import xarray as xr
import matplotlib.pyplot as plt

#
_CTSM_PYTHON = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "..", "python")
)
sys.path.insert(1, _CTSM_PYTHON)

from ctsm import add_cime_to_path
from CIME.utils import run_cmd_no_fail

from shutil import copy, rmtree

def main(build_dir):
  
  # build the test
   build_test(build_dir)
  
  # run the test
   #run_test("biogeophys/test/solar_rad_test", "FATES_rad_exe")
   # 
   # plot_output()

def build_test(build_dir):
  
  ## super hacky right now!!
  os.chdir("..")

  run_command = ["../../cime/scripts/fortran_unit_testing/run_tests.py",
                 "--build-dir", build_dir]
  run_cmd_no_fail(" ".join(run_command), combine_output=True)
  
  copy(os.path.join(build_dir, "__command_line_test__/__command_line_test__/fates_biogeophys_test/solar_rad_test/FATES_rad_exe"),
           "biogeophys/test/solar_rad_test/")
  
def run_test(dir, prog):
  
  # run the program
  os.chdir(dir)
  run_command = [f"./{prog}"]
  output = run_cmd_no_fail(" ".join(run_command), combine_output=True)
  print(output)
  
def plot_output():
  dat = xr.open_dataset('/Users/afoster/Documents/ncar/CTSM/src/fates/biogeophys/test/solar_rad_test/radiation_out.nc')
  plt.figure()
  dat.alb_dir.plot()
  plt.show()
  
if __name__ == "__main__":
  main("build.temp")