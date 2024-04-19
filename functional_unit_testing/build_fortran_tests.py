import os
import sys
import shutil

_FATES_PYTHON = os.path.join(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(1, _FATES_PYTHON)

from utils import add_cime_lib_to_path
add_cime_lib_to_path()

from CIME.utils import get_src_root, run_cmd_no_fail, expect, stringify_bool
from CIME.build import CmakeTmpBuildDir
from CIME.XML.machines import Machines
from CIME.BuildTools.configure import configure, FakeCase
from CIME.XML.env_mach_specific import EnvMachSpecific

_CIMEROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../../cime")

def run_cmake(name, test_dir, pfunit_path, netcdf_c_path, netcdf_f_path, cmake_args):
    """Run cmake for the fortran unit tests
    Arguments:
    name (str) - name for output messages
    test_dir (str) - directory to run Cmake in
    pfunit_path (str) - path to pfunit
    netcdf_c_path (str) - path to netcdf
    netcdf_f_path (str) - path to netcdff
    clean (bool) - clean the build
    """
    if not os.path.isfile("CMakeCache.txt"):
        print(f"Running cmake for {name}.")

        # directory with cmake modules
        cmake_module_dir = os.path.abspath(os.path.join(_CIMEROOT, "CIME", "non_py", "src", "CMake"))

        # directory with genf90
        genf90_dir = os.path.join(_CIMEROOT, "CIME", "non_py", "externals", "genf90")

        cmake_command = [
          "cmake",
          "-C Macros.cmake",
          test_dir,
          f"-DCIMEROOT={_CIMEROOT}",
          f"-DSRC_ROOT={get_src_root()}",
          f"-DCIME_CMAKE_MODULE_DIRECTORY={cmake_module_dir}",
          "-DCMAKE_BUILD_TYPE=CESM_DEBUG",
          f"-DCMAKE_PREFIX_PATH={pfunit_path}",
          "-DUSE_MPI_SERIAL=ON",
          "-DENABLE_GENF90=ON",
          f"-DCMAKE_PROGRAM_PATH={genf90_dir}"
        ]

        if netcdf_c_path is not None:
            cmake_command.append(f"-DNETCDF_C_PATH={netcdf_c_path}")

        if netcdf_f_path is not None:
            cmake_command.append(f"-DNETCDF_F_PATH={netcdf_f_path}")

        cmake_command.extend(cmake_args.split(" "))

        run_cmd_no_fail(" ".join(cmake_command), combine_output=True)

def find_library(caseroot, cmake_args, lib_string):
    """Find the library installation we'll be using, and print its path

    Args:
        caseroot (str): Directory with pfunit macros
        cmake_args (str): The cmake args used to invoke cmake (so that we get the correct makefile vars)
    """
    with CmakeTmpBuildDir(macroloc=caseroot) as cmaketmp:
        all_vars = cmaketmp.get_makefile_vars(cmake_args=cmake_args)

        all_vars_list = all_vars.splitlines()
        for all_var in all_vars_list:
            if ":=" in all_var:
                expect(all_var.count(":=") == 1, f"Bad makefile: {all_var}")
                varname, value = [item.strip() for item in all_var.split(":=")]
                if varname == lib_string:
                    return value

        expect(False, f"{lib_string} not found for this machine and compiler")

        return None

def prep_build_dir(build_dir, clean):
    """Create (if necessary) build directory and clean contents (if asked to)

    Args:
        build_dir (str): build directory name
        clean (bool): whether or not to clean contents
    """

    # create the build directory
    build_dir_path = os.path.abspath(build_dir)
    if not os.path.isdir(build_dir_path):
      os.mkdir(build_dir_path)

    # change into that directory
    os.chdir(build_dir_path)

    # clean up any files if we want to
    if clean:
        clean_cmake_files()

    return build_dir_path

def clean_cmake_files():
    """Deletes all files related to build

    """
    if os.path.isfile("CMakeCache.txt"):
        os.remove("CMakeCache.txt")
    if os.path.isdir("CMakeFiles"):
        shutil.rmtree("CMakeFiles")

    cwd_contents = os.listdir(os.getcwd())

    # Clear contents to do with cmake cache
    for file in cwd_contents:
      if (
          file in ("Macros.cmake", "env_mach_specific.xml")
          or file.startswith("Depends")
          or file.startswith(".env_mach_specific")
          ):
        os.remove(file)

def get_extra_cmake_args(build_dir, mpilib):
    """Makes a fake case to grab the required cmake arguments
    Args:
        build_dir (str): build directory name
        mpilib (str): MPI library name
    """
    # get the machine objects file
    machobj = Machines()

    # get compiler
    compiler = machobj.get_default_compiler()

    # get operating system
    os_ = machobj.get_value("OS")

    # Create the environment, and the Macros.cmake file
    #
    #
    configure(
        machobj,
        build_dir,
        ["CMake"],
        compiler,
        mpilib,
        True,
        "nuopc",
        os_,
        unit_testing=True,
    )
    machspecific = EnvMachSpecific(build_dir, unit_testing=True)

    # make a fake case
    fake_case = FakeCase(compiler, mpilib, True, "nuopc", threading=False)

    cmake_args = (
        "{}-DOS={} -DMACH={} -DCOMPILER={} -DDEBUG={} -DMPILIB={} -Dcompile_threaded={} -DCASEROOT={}".format(
            "",
            os_,
            machobj.get_machine_name(),
            compiler,
            stringify_bool(True),
            mpilib,
            stringify_bool(False),
            build_dir
        )
    )

    return cmake_args

def run_make(name, make_j, clean=False, verbose=False):
    """Run make in current working directory

    Args:
        name (str): Name for output messages
        make_j (int): number of processes to use for make
        clean (bool, optional): whether or not to clean Defaults to False.
        verbose (bool, optional): verbose error logging for make Defaults to False.
    """

    print(f"Running make for {name}.")

    if clean:
        run_cmd_no_fail("make clean")

    make_command = ["make", "-j", str(make_j)]

    if verbose:
        make_command.append("VERBOSE=1")

    run_cmd_no_fail(" ".join(make_command), combine_output=True)

def build_exists(build_dir, test_dir, test_exe=None):
    """Checks to see if the build directory and associated executables exist.

        Args:
          build_dir (str): build directory
          test_dir (str): test directory
          test_exe (str): test executable
    """

    build_path = os.path.abspath(build_dir)
    if not os.path.isdir(build_path):
        return False

    if not os.path.isdir(os.path.join(build_path, test_dir)):
        return False

    if test_exe is not None:
        if not os.path.isfile(os.path.join(build_path, test_dir, test_exe)):
            return False

    return True

def build_unit_tests(build_dir, name, cmake_directory, make_j, clean=False):
    """Build the unit test executables

    Args:
        build_dir (str): build directory
        name (str): name for set of tests
        cmake_directory (str): directory where the make CMakeLists.txt file is
        make_j (int): number of processes to use for make
        clean (bool, optional): whether or not to clean the build first. Defaults to False.
    """
    # create the build directory
    full_build_path = prep_build_dir(build_dir, clean=clean)

    # get cmake args and the pfunit and netcdf paths
    cmake_args = get_extra_cmake_args(full_build_path, "mpi-serial")
    pfunit_path = find_library(full_build_path, cmake_args, "PFUNIT_PATH")

    if not "NETCDF" in os.environ:
        netcdf_c_path = find_library(full_build_path, cmake_args, "NETCDF_C_PATH")
        netcdf_f_path = find_library(full_build_path, cmake_args, "NETCDF_FORTRAN_PATH")
    else:
        netcdf_c_path = None
        netcdf_f_path = None

    # change into the build dir
    os.chdir(full_build_path)

    # run cmake and make
    run_cmake(name, cmake_directory, pfunit_path, netcdf_c_path, netcdf_f_path, cmake_args)
    run_make(name, make_j, clean=clean)

