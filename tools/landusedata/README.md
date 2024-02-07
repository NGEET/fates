# FATES LUH2 data tool README

## Purpose

This tool takes the raw Land Use Harmonization (https://luh.umd.edu/), or LUH2, data files as
input and prepares them for use with FATES.  The tool concatenates the various raw data sets into
a single file and provides the ability to regrid the source data resolution to a target
resolution that the user designates.  The output data is then usable by FATES, mediated through
a host land model (currently either CTSM or E3SM).

For more information on how FATES utilizes this information see https://github.com/NGEET/fates/pull/1040.

## Installation

This tool requires the usage of conda with python3.  See https://docs.conda.io/en/latest/miniconda.html#installing
for information on installing conda on your system.  To install the conda environment necessary to run the tool
execute the following commands:

conda env create -f conda-luh2.yml

This will create a conda environment named "luh2".  To activate this environment run:

conda activate luh2

For more information on creating conda environments see
https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file

Note that it is planned that a subset of host land model (hlm) and hlm supported machines will incoporate this tool into the surface dataset workflow.
As such, if you are working on one of these machines, the output from this tool may be precomputed and available for the grid resolution of interest.

## Usage

After activating the "luh2" environment the tool can be run from the command line with the following minimum required inputs:

python luh2.py -l <raw-luh2-datafile> -s <luh2-static-datafile> -r <regrid-targetfile> -w <regridder-output> -o <outputfile>

The description of the minimum required input arguments is as follows:
- raw-luh2-datafile: this is one of three raw luh2 datafiles, either states, transitions, or management.  This is the data to be regridded and used by FATES.
- luh2-static-datafile: supplementary 0.25 deg resolution static data used in the construction of the raw luh2 datafiles.  This is utilized to help set the gridcell mask for the output file.
- regrid-targetfile: host land model surface data file intended to be used in conjunction with the fates run at a specific grid resolution.  This is used as the regridder target resolution.
- regridder-output: the path and filename to write out the regridding weights file or to use an existing regridding weights file.
- outputfile: the path and filename to which the output is written

The tool is intended to be run three times, sequentially, to concatenate the raw states, transitions, and management data into a single file.  After the first run of
the tool, a merge option should also be included in the argument list pointing to the most recent output file.  This will ensure that the previous regridding run
will be merged into the current run as well as reusing the previously output regridding weights file (to help reduce duplicate computation).
The luh2.sh file in this directory provides an example shell script in using the python tool in this sequential manner.  The python tool itself provides additional
help by passing the `--help` option argument to the command line call.

## Description of directory contents

- luh2.py: main luh2 python script
- luh2mod.py: python module source file for the functions called in luh2.py
- luh2.sh: example bash shell script file demonstrating how to call luh2.py
- conda-luh2.yml: conda enviroment yaml file which defines the minimum set of package dependencies for luh2.py
