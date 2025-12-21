# FATES Testing Framework

This directory contains the infratstructure to set up, build, and execute FATES functional
and unit tests.

## Test Definitions

* **Functional Tests**: Standalone Fortran programs that exercse specific modules of the
FATES production code. They often output NetCDF data and include a Python class for
automated plotting and science validation. These are "hands-on" tests for developers.

* **Unit Tests**: Discrete tests with a clear Pass/Fail status, utilizing **pfUnit** or
**CTest**. These are best for testing edge cases, error handling, and logical branches
within individual subroutines.

## Environment Setup

### Python Requirements

The framework requires Python 3.12+ with `xarray`, `pandas`, `scipy`, `netCDF`, `numpy`,
and `matplotlib`. You can create a compatible environment using the provided file:

```bash
conda env create --file=environment.yml
conda activate fates_testing
```

### System Requirements

While these tests do not require a host land model (e.g., CTSM/ELM), they require:

* **CIME and shr**: These must be present in your source tree.
* **Libraries**: NetCDF (C and Fortran), ESMF, and pfUnit.
* **Machine Tags**: Valid machine and compiler configurations (e.g., NCAR Derecho/Izumi).
If FATES is already running on your machine, these scripts should work out of the box.

See `docs/cime_setup.md` for a step-by-step guide to setting up the environment on your
own machine (Mac/Linux).

## Execution

### Running Functional Tests

```bash
# Run all functional tests using default parameters
./run_functional_tests.py

# Run specific tests with a custom parameter file
./run_functional_tests.py --test-list allometry,fire --param-file my_params.json
```

### Running Unit Tests

```bash
# Run all unit tests
./run_unit_tests.py
```

*Note: Use the `--help` flag on either script to see full options for build directories,
skipping runs, and saving figures.*

## Creating New Tests

To maintain consistency, always use the boilerplate generator to start a new test.

### Step 1. Generate Scaffolding

Run the generator script. This updates the CMake system, adds config entries, and creates
your test directory.

```bash
# Example: Create a new functional test named 'hydro_stress'
./generate_empty_test.py functional --test-name hydro_stress
```

### Step 2: Implement the Fortran Logic

Navigate to `tests/[unit|functional]/[test_name]` and edit the generated `.F90` or `.pf`
file.

* **Functional**: Write a standard Fortran program.
* **Unit**: Write a pfUnit-compatible test module.

### Step 3: Configure Metadata

Edit `config/functional.cfg` or `config/unit.cfg`. The generator provides defaults, but
you may need to update:

* `out_file`: The name of the NetCDF file your Fortran code generates.
* `use_param_file`: Set to `True` to pass the FATES JSON parameter file to your binary.
* `datm_file`: (Optional) Path to a driver file, relative to the FATES root.

### Step 4: Python Analysis (Functional Only)

The generator creates a Python file (`[test_name]_test.py`) in your test directory. Edit
the `plot_output` method to define how your test results should be visualized or validated.

## Directory Structure

* `framework/`: Core logic for loading, building, and executing tests.
* `config/`: Registry of all active tests and their arguments.
* `tests/`: The source code for individual test cases.
* `templates/`: Boilerplate files used by the generator.

## Troubleshooting

* **Missing Drivers**: If a functional tests requires a `datm_file`, ensure the path in
the `.cfg` is correct relative to the FATES root. The script will perform a pre-flight
check and fail if the file is missing.

* **Build Failures**: Ensure your environment modules (compiler, netcdf, esmf) match the
configuration in your CIME machine tags.

* **Library Paths**: If the executable fails to start, ensure your `LD_LIBRARY_PATH`
(or `DYLD_LIBRARY_PATH` on macOS) includes the paths to the NetCDF and ESMF shared
libraries.

