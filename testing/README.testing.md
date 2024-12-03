# FATES Testing

These scripts set up, build, and run FATES functional tests and unit tests.

By "functional" test, we mean a standalone Fortran program that runs pieces of the FATES 
production code, potentially outputs results (i.e. to a netcdf file), and potentially 
plots or runs some other test on the output in python. These tests do not necessarily
have a pass/fail status, but are meant to be more hands-on for the user.

Unit tests do have a pass/fail outcome, and are written as such. We accommodate both Ctest
and pfunit tests. 

## How to run 

To run the testing scripts, `run_functional_tests.py` or `run_unit_tests.py`, you will
need a few python packages. You can create a conda environment with these packages
using the `testing.yml` file: `conda env create --file=testing.yml`

Though these tests do not need host land model code, you will need the `cime` and `shr`
repositories, as well as a machine configuration file. If you are already set up to run
FATES on a machine (e.g. derecho at NCAR), you should be able to run these scripts out
of the box on that machine.

Additionally, these tests require netcdf and netcdff, as well as a fortran compiler (e.g. gnu),
esmf, and pfunit. See `cime_setup.md` for tips on how to do this.

Once you are set up, you should be able to just run the scripts. For the functional test
script, you can point to your own parameter file (cdl or nc; `./run_functional_tests -f my_param_file.nc`). 
If none is supplied the script will use the default cdl file in `parameter_files`.

You can run an individual set of tests by passing the script a comma-separated list of 
test names. See the `functional_tests.cfg` or `unit_tests.cfg` for the test names. If you
do not supply a list, the script will run all tests.

## How to create new tests

First, determine if you are going to create a functional or unit test. Remember,
unit tests must have a pass/fail outcome. These are best for testing edgecases, error
handling, or where we know exactly what the result should be from a method.

First, add your test to either the `functional_tests.cfg` or `unit_tests.cfg` config file,
depending on the test you want to create.

### Config file information

The `test_dir` is where cmake will place relevant libraries created from your test.
The convention is to call this directory `fates_{testname}_ftest` for functional tests
and `fates_{testname}_utest` for unit tests.

The `test_exe` (only applicable for functional tests) is the executable the script will
create based on your test program. The convention is to call it `FATES_{testname}_exe`.

The `out_file` (only applicable for functional tests) is the output file name that your test
may or may not create. Set the value to `None` if your test does not create one.

Set `use_param_file` to `True` if your test uses the FATES parameter file, and `False` 
otherwise. This is only applicable for functional tests.

Add any other arguments your test needs in the `other_args` list.
This is only applicable for functional tests.

### Cmake setup

Under the `testing/functional_testing` or `testing/unit_testing` directory, depending
on your test type, create a new directory for your test, e.g. "my_new_test".

In the file `testing/CMakeLists.txt` add your test directory to the list of tests, e.g.:

`add_subdirectory(functional_testing/my_new_test fates_new_test_ftest)`

The first argument must match your directory name, and the second must match the
`test_dir` value you set up in the config file above.

Inside your new testing directory create a new `CMakeLists.txt` file to tell cmake
how to build and compile your program. It may be easiest to copy an existing one
from another similar test and then update the relevant information.

Importantly, the sources must be set to tell the program which file(s) contain your
test program. Additionally, for functional tests, the executable name must match what 
you set up in the config file above.

### Fortran program

Write your Fortran tests. For functional tests, this should be an actual Fortran `program`.
For unit tests this will be either a ctest or a pfunit test. See existing tests for examples.

For functional tests, if you output an output file, the name must match what you set up
in the config file above.

### Python setup - functional tests

For functional tests, you will need to add your test as a new concrete class based on
the abstract FunctionalTest class. See examples in the `functional_testing` directory. 
Most of the work involves creating a plotting function for your test.

You will then need to add this class as an import statment at the top of the 
`run_functional_tests.py` script.
