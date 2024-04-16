#!/usr/bin/env python3

import os
import sys
from build_fortran_tests import build_unit_tests

_FATES_PYTHON = os.path.join(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(1, _FATES_PYTHON)

from utils import add_cime_lib_to_path
add_cime_lib_to_path()

from CIME.utils import run_cmd_no_fail


if __name__ == "__main__":
    
    ## Arguments
    clean = True
    build = True
    build_dir = "../_build"
    name = "fates_unit_tests"
    make_j = 8
    cmake_directory = os.path.abspath("../")
    param_file = None
    param_cdl_path = "../parameter_files/fates_params_default.cdl"
    
    ## Constants
    test_dir = "fates_allom_test"
    test_exe = "FATES_allom_exe"
    
    if build:
        build_unit_tests(build_dir, name, cmake_directory, make_j, clean=clean)
    
    build_dir_path = os.path.abspath(build_dir)    
    exe_path = os.path.join(build_dir_path, test_dir, test_exe)
    
    if param_file is None:
      file_basename = os.path.basename(param_cdl_path).split(".")[-2]
      file_nc_name = f"{file_basename}.nc"
      file_gen_command = [
        "ncgen -o",
        os.path.join(file_nc_name),
        param_cdl_path
      ]
      run_cmd_no_fail(" ".join(file_gen_command), combine_output=True)
          
    out = run_cmd_no_fail(exe_path, combine_output=True)
    print(out)