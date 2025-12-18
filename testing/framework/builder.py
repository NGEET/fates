"""
Builds/compiles any tests within the FATES repository
"""

import os
import shutil
import logging
from pathlib import Path
from dataclasses import dataclass, field
from typing import Optional

# local imports
from framework.utils.path import add_cime_lib_to_path, path_to_cime

# initialize CIME path
add_cime_lib_to_path()

try:
    from CIME.utils import get_src_root, run_cmd_no_fail
    from CIME.build import CmakeTmpBuildDir
    from CIME.XML.machines import Machines
    from CIME.BuildTools.configure import configure, FakeCase
    from CIME.XML.env_mach_specific import EnvMachSpecific
except ImportError as e:
    # fail fast if environment isn't set up
    raise ImportError(
        f"CIME dependencies missing. Ensure environment is configured. Error: {e}"
    ) from e

# setup logging
logging.basicConfig(
    level=logging.WARNING, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

# constants
_CIMEROOT = path_to_cime()
_MPI_LIBRARY = "mpi-serial"


@dataclass
class BuildConfig:
    """Encapsulates all configuration needed for the build."""

    build_dir: Path
    cmake_dir: Path
    make_j: int
    clean: bool = False
    verbose: bool = False
    compiler: str = field(init=False)
    machine_os: str = field(init=False)
    mpilib: str = _MPI_LIBRARY

    # paths resolved during setup
    pfunit_path: Optional[str] = None
    netcdf_c_path: Optional[str] = None
    netcdf_f_path: Optional[str] = None
    cmake_args: str = ""

    def __post_init__(self):
        """Ensure paths are Path objects"""
        self.build_dir = Path(self.build_dir).resolve()
        self.cmake_dir = Path(self.cmake_dir).resolve()


class TestBuilder:
    """Orchestrates the build process."""

    def __init__(self, config: BuildConfig):
        self.config = config

        if self.config.verbose:
            logging.getLogger().setLevel(logging.INFO)
        else:
            logging.getLogger().setLevel(logging.WARNING)

    def build(self):
        """Main entry point for building tests."""

        self._prep_build_dir()

        # gether environment info and cmake args
        self._configure_environment()

        # execute build steps
        self._run_cmake()
        self._run_make()

    def _prep_build_dir(self):
        """Creates directory and cleans if requested."""
        if not self.config.build_dir.exists():
            self.config.build_dir.mkdir(parents=True, exist_ok=True)

        if self.config.clean:
            logger.info("Cleaning build directory: %s", self.config.build_dir)
            self._clean_cmake_files()

    def _clean_cmake_files(self):
        """Deletes all files related to cmake build."""
        to_remove = ["CMakeCache.txt", "Macros.cmake", "env_mach_specific.xml"]
        to_remove_dirs = ["CMakeFiles"]

        for f in to_remove:
            p = Path(f)
            if p.is_file():
                p.unlink()

        for d in to_remove_dirs:
            p = Path(d)
            if p.is_dir():
                shutil.rmtree(p)

        for p in Path(".").glob("Depends*"):
            p.unlink()
        for p in Path(".").glob(".env_mach_specific*"):
            p.unlink()

    def _configure_environment(self):
        """Generates cmake arguments and finds library paths
        NOTE: this section interacts heavily with CIME global state/objects
        """

        machobj = Machines()
        self.config.compiler = machobj.get_default_compiler()
        self.config.machine_os = machobj.get_value("OS")

        # configure CIME (creates Macros.cmake)
        configure(
            machobj,
            str(self.config.build_dir),
            ["CMake"],
            self.config.compiler,
            self.config.mpilib,
            True,
            "nuopc",
            self.config.machine_os,
            unit_testing=True,
        )

        # load environment specific to this machine
        machspecific = EnvMachSpecific(str(self.config.build_dir), unit_testing=True)
        fake_case = FakeCase(
            self.config.compiler, self.config.mpilib, True, "nuopc", threading=False
        )
        machspecific.load_env(fake_case)

        self._generate_cmake_args(machobj)
        self._find_libraries()

    def _generate_cmake_args(self, machobj):
        args_list = [
            f"-DOS={self.config.machine_os}",
            f"-DMACH={machobj.get_machine_name()}",
            f"-DCOMPILER={self.config.compiler}",
            f"-DDEBUG={'ON'}",
            f"-DMPILIB={self.config.mpilib}",
            f"-Dcompile_threaded={'OFF'}",
            f"-DCASEROOT={self.config.build_dir}",
        ]
        self.config.cmake_args = " ".join(args_list)

    def _find_libraries(self):
        """Locates PFUNIT and NETCDF paths."""
        self.config.pfunit_path = self._query_makefile_var("PFUNIT_PATH")

        if "NETCDF" not in os.environ:
            self.config.netcdf_c_path = self._query_makefile_var("NETCDF_C_PATH")
            self.config.netcdf_f_path = self._query_makefile_var("NETCDF_FORTRAN_PATH")

    def _query_makefile_var(self, var_name: str) -> Optional[str]:
        """Helper to query variables from CIME makefile generation

        Args:
            var_name (str): variable to query

        Returns:
            Optional[str]: library path
        """
        # using the CmakeTmpBuildDir context manager provided by CIME
        with CmakeTmpBuildDir(macroloc=str(self.config.build_dir)) as cmaketmp:
            all_vars = cmaketmp.get_makefile_vars(cmake_args=self.config.cmake_args)

            # parsing logic
            for line in all_vars.splitlines():
                if ":=" in line:
                    parts = line.split(":=")
                    if len(parts) == 2:
                        key, val = parts[0].strip(), parts[1].strip()
                        if key == var_name:
                            return val
        logger.warning("Library path variable %s not found.", var_name)
        return None

    def _run_cmake(self):
        """Run cmake"""
        cmake_module_dir = _CIMEROOT / "CIME/non_py/src/CMake"
        genf90_dir = _CIMEROOT / "CIME/non_py/externals/genf90"

        cmd = [
            "cmake",
            "-C",
            "Macros.cmake",
            str(self.config.cmake_dir),
            f"-DCIMEROOT={_CIMEROOT}",
            f"-DSRC_ROOT={get_src_root()}",
            f"-DCIME_CMAKE_MODULE_DIRECTORY={cmake_module_dir}",
            "-DCMAKE_BUILD_TYPE=CESM_DEBUG",
            f"-DCMAKE_PREFIX_PATH={self.config.pfunit_path}",
            "-DUSE_MPI_SERIAL=ON",
            "-DENABLE_GENF90=ON",
            f"-DCMAKE_PROGRAM_PATH={genf90_dir}",
        ]
        if self.config.netcdf_c_path:
            cmd.append(f"-DNETCDF_C_PATH={self.config.netcdf_c_path}")
        if self.config.netcdf_f_path:
            cmd.append(f"-DNETCDF_F_PATH={self.config.netcdf_f_path}")

        # append extra args
        cmd.extend(self.config.cmake_args.split(" "))

        logger.info("Running cmake...")
        run_cmd_no_fail(
            " ".join(cmd), from_dir=self.config.build_dir, combine_output=True
        )

    def _run_make(self):
        """Run make"""
        cmd = ["make", "-j", str(self.config.make_j)]
        if self.config.verbose:
            cmd.append("VERBOSE=1")

        logger.info("Running make...")
        run_cmd_no_fail(
            " ".join(cmd), from_dir=self.config.build_dir, combine_output=True
        )


def build_tests(
    build_dir: Path,
    cmake_dir: Path,
    make_j: int,
    clean: bool = False,
    verbose: bool = False,
):
    """Wrapper function

    Args:
        build_dir (str): _description_
        cmake_directory (str): _description_
        make_j (int): _description_
        clean (bool, optional): _description_. Defaults to False.
        verbose (bool, optional): _description_. Defaults to False.
    """
    config = BuildConfig(
        build_dir=build_dir,
        cmake_dir=cmake_dir,
        make_j=make_j,
        clean=clean,
        verbose=verbose,
    )
    builder = TestBuilder(config)
    builder.build()


def build_exists(build_dir: str, test_dir: str, test_exe: Optional[str] = None) -> bool:
    """Checks to see if build artifacts exist.

    Args:
      build_dir (str): build directory
      test_dir (str): test directory
      test_exe (str, optional): test executable. Defaults to None.
    Returns:
      bool: whether or not build directory and associated executables exist
    """

    build_path = Path(build_dir).resolve()
    test_path = build_path / test_dir

    if not test_path.is_dir():
        return False

    if test_exe and not (test_path / test_exe).is_file():
        return False

    return True
