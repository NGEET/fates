#!/usr/bin/env python

# =======================================================================================
# =======================================================================================

import os
import argparse
import shutil
import tempfile
import sys
import datetime
import time
import numpy as np
import code  # For development: code.interact(local=dict(globals(), **locals()))
import json

# ========================================================================================
# ========================================================================================
#                                        Main
# ========================================================================================
# ========================================================================================

def main():

    parser = argparse.ArgumentParser(description='Parse command line arguments to this script.')
    #
    parser.add_argument('--fin', '--input', dest='inputfname', type=str, help="Input filename.  Required.", required=True)
    parser.add_argument('--fout','--output', dest='outputfname', type=str, help="Output filename. Required if not --overwrite.")
    parser.add_argument('--overwrite', dest='overwrite', help="Use this flag to overwrite the input file",action="store_true")
    parser.add_argument('--pft-indices',dest='pft_indices',nargs='+',help="Integer list of pfts to keep in output file.",required=True)
    parser.add_argument('--silent',dest='silent',help="Suppress output messages",action="store_true")
                        
    args = parser.parse_args()

    with open(args.inputfname, 'r') as file:
        data = json.load(file)

        maxpft = data['dimensions']['fates_pft']

        pft_indices = []
        for index in args.pft_indices:
            for item in index.strip('=[]').split(','):
                if(int(item) > maxpft):
                    print(f'\nYou provided an index that is larger than')
                    print(f'the size of the variable\'s dataset.')
                    print(f'index: {int(item)}, total indices: {int(max_index)}')
                    print(f'You arguments for --indices were {args.indices}.')
                    print(f'Exiting\n')
                    exit(2)
                else:
                    pft_indices.append(int(item)-1)
        
        # Find all variables with fates_pft in the dims...
        for key, struct in data['variables'].items():
            if('fates_pft' in struct['dims']):
                if(not args.silent):
                    print(f'Changing pft dimensions for: {key}')
                np_array = np.array(struct['data'])
                if(len(struct['dims'])>1):
                    # fates_pft is always the second dimension
                    subset_array = np_array[:, pft_indices]
                else:
                    subset_array = np_array[pft_indices]
                data['variables'][key]['data'] = subset_array.tolist()

        maxpft = data['dimensions']['fates_pft'] = len(pft_indices)

    if args.outputfname is not None:
        output_filename = args.outputfname
        if(args.overwrite is True):
            print(f'\nYou have specified both an output file, and to overwrite.')
            print(f'You can\'t do both, you must pick one.')
            print(f'see: --output and --overwrite')
            print(f'Exiting\n')
            exit(2)
    else:
        if(args.overwrite is True):
            output_filename = args.inputfname
            
    
    with open(output_filename, 'w') as outfile:
    
        # 2. Use json.dump(data, file_object) to write the Python object
        #    directly to the file as JSON format.
        # indent=4 is highly recommended for readability (pretty-printing).
        json.dump(data, outfile, indent=2)     

        

# =======================================================================================
# This is the actual call to main

if __name__ == "__main__":
    main()

