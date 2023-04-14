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

def main():
  os.chdir("..")
  if os.path.isdir("build"):
    rmtree("build")

  ## super hacky right now!!
  run_command = ["../../cime/scripts/fortran_unit_testing/run_tests.py",
                 "--build-dir", "build"]
  run_cmd_no_fail(" ".join(run_command), combine_output=True)
  
  copy("build/__command_line_test__/__command_line_test__/unit_testing/test/radiation_test/FATES_rad_test",
           "unit_testing/test/radiation_test/")
  
  os.chdir("unit_testing/test/radiation_test")
  run_command = ["./FATES_rad_test"]
  output = run_cmd_no_fail(" ".join(run_command), combine_output=True)
  print(output)
  
  
if __name__ == "__main__":
  main()