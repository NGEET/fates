#!/usr/bin/env python

### This script takes a ED2 style inventory init file and converts it to a file compatible with FATES.
# It accepts the following flags:
# --type : patch or cohort
# --fin : input filename
# --fout : output file name

import argparse
import pandas as pd
import sys

def main():
    parser = argparse.ArgumentParser(description='Parse command line arguments to this script.')
    #
    parser.add_argument('--type', dest='fatestype', type=str, help="patch or cohort. Required.", required=True)
    parser.add_argument('--fin', dest='fnamein', type=str, help="Input filename.  Required.", required=True)
    parser.add_argument('--fout', dest='fnameout', type=str, help="Output filename.  Required.", required=True) 

    args = parser.parse_args()

    # open the input data
    dsin = pd.read_csv(args.fnamein, delim_whitespace=True)

    # if patch file delete unnecessary patch columns
    if args.fatestype == 'patch' :
        keep_col = ['time', 'patch', 'trk', 'age', 'area']
        newds = dsin[keep_col]
        

    # if cohort file delete unnecessary cohort columns
    elif args.fatestype == 'cohort' :
        keep_col = ['time', 'patch', 'dbh', 'pft', 'nplant']
        newds = dsin[keep_col]
        
    else :
        print("type must be one of patch or cohort")


    newds.to_csv(args.fnameout , index=False, sep=' ')    
# ========================================================================================================
# This is the actual call to main

if __name__ == "__main__":
    main()
    
