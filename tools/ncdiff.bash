#!/usr/bin/env bash

## script that compares the differences between two netcdf files.
## two arguments are the paths to two files to compare

tempfile1=$(mktemp)
tempfile2=$(mktemp)

ncdump $1 >> ${tempfile1}
ncdump $2 >> ${tempfile2}

diff  ${tempfile1} ${tempfile2}

rm ${temp_file1}
rm ${temp_file2}
