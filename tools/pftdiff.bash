#!/usr/bin/env bash

## script to compare two PFTs in a FATES parameter file.  takes three arguments: 
## first argument is the parameter file name
## second argument is the first pft
## third argument is the second pft

tempfile1=$(mktemp)
tempfile2=$(mktemp)
tempfile3=$(mktemp)
tempfile4=$(mktemp)

toolsdir=$(dirname "$0")

$toolsdir/FatesPFTIndexSwapper.py --pft-indices=$2 --fin=$1 --fout=${tempfile1} &>/dev/null
$toolsdir/FatesPFTIndexSwapper.py --pft-indices=$3 --fin=$1 --fout=${tempfile2} &>/dev/null

ncdump ${tempfile1} >> ${tempfile3}
ncdump ${tempfile2} >> ${tempfile4}

diff  ${tempfile3} ${tempfile4}

rm ${tempfile1}
rm ${tempfile2}
rm ${tempfile3}
rm ${tempfile4}
