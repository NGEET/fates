# Instructions for setting up CIME on your personal computer

## Mac and Linux Users

### Downloads and Installs

1. *For Mac Users Only*: Install Apple Developer Tools, if you haven't already
2. Install homebrew ([link here](https://brew.sh/))

### PFunit

Clone pfunit from GitHub:

```bash
git clone https://github.com/Goddard-Fortran-Ecosystem/pFUnit.git
```

Install and build:

```bash
export PFUNIT=/Users/afoster/pfunit_builds/pfunit4.11.1-gfortran14.2.0

mkdir build-dir
cd build-dir
cmake -DSKIP_OPENMP=YES -DCMAKE_INSTALL_PREFIX=$PFUNIT/install ..
make -j 8

make tests
make install
```

#### Homebrew Installs

```bash
brew install subversion git bash-completion

brew install kdiff3

brew install gcc

brew install netcdf 

brew install nco

brew install ncview

brew install mpich

brew install cmake

brew install markdown

brew install sloccount

brew install pyqt --with-python3

brew install lapack
```

For compilers to find `lapack` you may need to set:

```bash
export LDFLAGS="-L/usr/local/opt/lapack/lib"
export CPPFLAGS="-I/usr/local/opt/lapack/include"
```

For compilers to find `libomp` you may need to set:

```bash
export LDFLAGS="-L/usr/local/opt/libomp/lib"
export CPPFLAGS="-I/usr/local/opt/libomp/include"
```

```bash
brew install highlight

brew install git-when-merged

brew install the_silver_searcher
```

*For Mac users*: `brew cask install mactex`

### ESMF Installation

Download ESMF: ([link here](https://github.com/esmf-org/esmf/releases))

Set an environment variable `ESMF_DIR` to where you want to install ESMF, for example: `export ESMF_DIR=/Users/afoster/esmf/esmf-8.4.0`

Next set up some other environment variables:

```bash
export ESMF_INSTALL_PREFIX=$ESMF_DIR/install_dir
export ESMF_COMM=mpich
export ESMF_COMPILER=gfortranclang
```

Inside the download, run:

```bash
gmake -j4 lib
gmake install
```

### CIME Setup

You'll need to set up a `.cime` directory with some specific files in it. I use Bill Sack's setup:

```bash
cd
git clone https://github.com/billsacks/mac_cime_configuration.git .cime
git checkout -b "my_configuration"
```

You'll need to modify some of the files and folders.

1. In the top-level `config_machines.xml`, change the `MACH` values to match your machine name (can obtain this via the terminal by typing `hostname`)
2. Rename the `green` directory to your computer's hostname.
3. In the `{hostname}/config_machines.xml`, update relevant fields like `MACH`, `DESC`, `OS`, and `SUPPORTED_BY`
4. In `{hostname}/config_machines.xml`, update the `ESMFMKFILE` based on the path to what you just made above.
5. Rename `gnu_green.cmake` to `gnu_{hostname}.cmake`.
6. Inside `gnu_{hostname}.cmake` update `NETCDF_C_PATH`, `NETCDF_FORTRAN_PATH`, `APPEND LDFLAGS` (both), and `PFUNIT_PATH` to your netcdf and pfunit paths (see below).

#### Libraries

The `NETCDF_C_PATH` should be the output of `nc-config --prefix`.
The `NETCDF_FORTRAN_PATH` should be the output of `nf-config --prefix`
Then update the `APPEND LDFLAGS` section as the output from:

```bash
nc-config --libs
```

```bash
nf-config --flibs 
```

Pfunit should be in a directiory called `installed` in the build directory.
