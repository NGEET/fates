import os 
import sys

#
_CTSM_PYTHON = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "..", "python")
)
sys.path.insert(1, _CTSM_PYTHON)

from ctsm import add_cime_to_path
from CIME.utils import run_cmd_no_fail

from shutil import rmtree, copy

def main(build_dir):
  os.chdir("..")
  if os.path.isdir(build_dir):
    rmtree(build_dir)

  ## super hacky right now!!
  run_command = ["../../cime/scripts/fortran_unit_testing/run_tests.py",
                 "--build-dir", "build.temp"]
  run_cmd_no_fail(" ".join(run_command), combine_output=True)
  
  copy(os.path.join(build_dir, "__command_line_test__/__command_line_test__/fates_biogeophys_test/solar_rad_test/FATES_rad_exe"),
           "biogeophys/test/solar_rad_test/")
  
  os.chdir("biogeophys/test/solar_rad_test")
  run_command = ["./FATES_rad_exe"]
  output = run_cmd_no_fail(" ".join(run_command), combine_output=True)
  print(output)
  
  
if __name__ == "__main__":
  main("build.temp")