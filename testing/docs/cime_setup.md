# CIME and Environment Setup Guide

This guide covers setting up the specialized Fortran environment required to build and
run FATES standalone tests on a local machine (Mac/Linux).

## System Dependencies (Homebrew)

FATES requires a specific stack of scientific libraries. If you are on a Mac,
use ([Homebrew](https://brew.sh/)) to install these core dependencies:

```bash
# Core Build Tools
brew install cmake gcc mpich git subversion

# NetCDF Stack (Critical for FATES)
brew install netcdf nco ncview

# Scientific Libraries
brew install lapack
```

## Setting up pFunit

Standalone Unit Tests require pFUnit. Important: It must be built with the same
compiler you intend to use for FATES.

**1. Clone and Build:**

```bash
git clone https://github.com/Goddard-Fortran-Ecosystem/pFUnit.git
cd pFUnit
mkdir build && cd build

# Set your install prefix to a local directory
export PFUNIT_INSTALL=$HOME/software/pfunit
cmake -DSKIP_OPENMP=YES -DCMAKE_INSTALL_PREFIX=$PFUNIT_INSTALL ..
make -j 8
make tests
make install
```

**2. Verify:** Ensure `$PFUNIT_INSTALL/bin/pfunit-config` exists.

## CIME Configuration

CIME acts as the "manager" that tells the FATES build system where your libraries are.
We use a local `.cime` directory to store these "Machine Tags."

### Initialze the .cime directory

If you don't have a `.cime` folder in your home directory, create one:

```bash
cd ~
git clone https://github.com/billsacks/mac_cime_configuration.git .cime
cd .cime
```

### Configure `config_machines.xml`

Find your computer's hostname by typing `hostname -s` in your terminal. Edit
`config_machines.xml` and update the following for your machine entry:

* `<MACH>`: Set this to your hostname.
* `<OS>`: `Darwin` (for Mac) or `Linux`.
* `<COMPILERS>`: `gnu`.
* `<ESMFMKFILE>`: The path to your `ESMF.mod` file (usually found in your ESMF installation directory). See below.

### Create your Compiler File

Rename the template to match your machine:

```bash
cp cmake_macros/gnu_green.cmake cmake_macros/gnu_$(hostname -s).cmake
```

Edit this new `.cmake` file and update the paths to match your Homebrew/Manual installs.
You can find these paths using these commands:

| Variable                | Command to find value                        |
| ----------------------- | -------------------------------------------- |
| **NETCDF_C_PATH**       | `nc-config --prefix`                         |
| **NETCDF_FORTRAN_PATH** | `nf-config --prefix`                         |
| **LDFLAGS**             | `nc-config --libs and nf-config --flibs`     |
| **PFUNIT_PATH**         | The `$PFUNIT_INSTALL` path used in Section 2 |

### Environment Variables

To ensure the FATES building scripts can communicate with CIME, you must define your
machine name and model type in your shell profile (e.g., `~/.zshrc` or `~/.bashrc`).

Add the following lines to the bottom of your profile:

```bash
# Tells CIME which machine entry to use from config_machines.xml
export CIME_MACHINE='hostname' 

# Tells CIME to use the CESM build structure logic
export CIME_MODEL='cesm'
```

*Note: Replace 'hostname` with the hostname you defined in `config_machines.xml`.*

After saving, remember to source your profile: `source ~/.zshrc`.

## ESMF Build

**1. Download:** [ESMF Releases](https://github.com/esmf-org/esmf/releases)

**2. Environment Variables:**

```bash
export ESMF_DIR=$(pwd)  # where you will install ESMF
export ESMF_INSTALL_PREFIX=$ESMF_DIR/install_dir
export ESMF_COMM=mpich
export ESMF_COMPILER=gfortranclang  # For Mac
```

**3. Build:**

```bash
gmake -j4 lib
gmake install
```

**4. Link:** Update the `<ESMFMKFILE>` path in your .cime/config_machines.xml to point 
to `$ESMF_DIR/install_dir/lib/esmf.mk`.
