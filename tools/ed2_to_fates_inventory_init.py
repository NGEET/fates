#!/usr/bin/env python

### This script takes a ED2 style inventory init file and converts it to a file compatible with FATES.
# It accepts the following flags:
# --fin : input filename

import argparse
import pandas as pd
import sys

def main():
    parser = argparse.ArgumentParser(description='Parse command line arguments to this script.')
    #
    parser.add_argument('--fin', dest='fnamein', type=str, help="Input filename.  Required.", required=True)
  
    args = parser.parse_args()

    # is it a pss or css file?
    filetype  = args.fnamein.split('.')[1]

    # make the new file name
    base_filename = args.fnamein.split('.')[0]
    output_filename = f"{base_filename}_{'fates'}.{filetype}"

    # open the input data
    dsin = pd.read_csv(args.fnamein, delim_whitespace=True)

    # if patch file delete unnecessary patch columns
    if filetype == 'pss' :
        keep_col = ['time', 'patch', 'trk', 'age', 'area']
        newds = dsin[keep_col]
        

    # if cohort file delete unnecessary cohort columns
    elif filetype == 'css' :
        keep_col = ['time', 'patch', 'dbh', 'pft', 'nplant']
        newds = dsin[keep_col]
        
    else :
        print("type must be one of patch or cohort")


    newds.to_csv(output_filename, index=False, sep=' ')    
# ========================================================================================================
# This is the actual call to main

if __name__ == "__main__":
    main()
    
