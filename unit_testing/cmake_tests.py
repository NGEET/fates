#!/usr/bin/env python3

import os
import sys
from shutil import rmtree

def main():
    cmake_stage("radiation test", 
                "/Users/afoster/Documents/ncar/CTSM/src/fates/unit_testing", 
                "test_pfunit_path")

def cmake_stage(
    name,
    test_spec_dir,
    pfunit_path,
    cmake_args=None,
    clean=False,
):
    """Run cmake in the current working directory.

    Arguments:
    name - Name for output messages.
    test_spec_dir - Test specification directory to run CMake on.
    """
    if clean:
        if os.path.isfile("CMakeCache.txt"):
            os.remove("CMakeCache.txt")
        if os.path.isdir("CMakeFiles"):
            rmtree("CMakeFiles")

    if not os.path.isfile("CMakeCache.txt"):

        print("Running cmake for " + name + ".")

        cmake_command = [
            "cmake",
            "-C Macros.cmake",
            test_spec_dir,
           # "-DCIMEROOT=" + _CIMEROOT,
           # "-DSRC_ROOT=" + get_src_root(),
           # "-DCIME_CMAKE_MODULE_DIRECTORY="
           # + os.path.abspath(os.path.join(_CIMEROOT, "CIME", "non_py", "src", "CMake")),
           # "-DCMAKE_BUILD_TYPE=" + build_type,
           # "-DPFUNIT_MPIRUN='" + mpirun_command + "'",
           # "-DPFUNIT_PATH=" + pfunit_path,
        ]

        if cmake_args is not None:
            cmake_command.extend(cmake_args.split(" "))

        #run_cmd_no_fail(" ".join(cmake_command), combine_output=True)
        print(" ".join(cmake_command))

if __name__ == "__main__":
    main()