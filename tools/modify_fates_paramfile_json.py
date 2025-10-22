#!/usr/bin/env python

#
# Original: C. Koven, 2018
# Refactored for json: R. Knox 2025

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
    parser.add_argument('--var','--variable', dest='varname', type=str, help="What variable to modify? Required.",required=True)
    parser.add_argument('--q','--queryvar', dest='queryvar', help="Report variable info only, then exit.",action="store_true")
    parser.add_argument('--indices', nargs='+', help='List of variable indices (1-len) to change (flattens 2d arrays, second dim is the inner).')
    parser.add_argument('--pft-name.')
    parser.add_argument('--values', nargs='+', help='List of values to change, must coincide with indices')
    parser.add_argument('--overwrite', dest='overwrite', help="Use this flag to overwrite the input file",action="store_true")
    parser.add_argument('--silent',dest='silent',help="Suppress output messages",action="store_true")

    args = parser.parse_args()
    
    with open(args.inputfname, 'r') as file:
        data = json.load(file)

        try:
            dim_names = data['variables'][args.varname]['dims']
        except KeyError:
            print(f'\n\nThe variable name "{args.varname}" is not in the dataset. Exiting.\n')
            exit(2)
            
        dim_sizes = list(np.shape(data['variables'][args.varname]['data']))
        max_index = np.prod(dim_sizes)
        vardata   = data['variables'][args.varname]['data']
        
        if(len(dim_sizes)>1):
            dtype = type(data['variables'][args.varname]['data'][0][0])
        else:
            dtype = type(data['variables'][args.varname]['data'][0])
        
        if(args.queryvar):
            print(f'Reporting: Variable: {args.varname}:')
            print(f'           dimension names: {dim_names}')
            print(f'           dimension sizes: {dim_sizes}')
            print(f'           Current Values:')
            if(len(dim_sizes)>1):
                for i in range(dim_sizes[0]):
                    print(f'           {data["variables"][args.varname]["data"][i][:]}')
            else:
                print(f'           {data["variables"][args.varname]["data"][:]}')
            print(f'           Index Pattern:')
            if(len(dim_sizes)>1):
                for i in range(dim_sizes[0]):
                    i0 = i*dim_sizes[1]+1
                    i1 = i0+dim_sizes[1]
                    #code.interact(local=dict(globals(), **locals()))
                    print(f'           {list(range(i0,i1))}')
            else:
                #code.interact(local=dict(globals(), **locals()))
                print(f'           {list(range(1,dim_sizes[0]+1))}')
            exit(0)

        indices = []
        for index in args.indices:
            for item in index.strip('=[]').split(','):
                if(int(item) > max_index):
                    print(f'\nYou provided an index that is larger than')
                    print(f'the size of the variable\'s dataset.')
                    print(f'index: {int(item)}, total indices: {int(max_index)}')
                    print(f'You arguments for --indices were {args.indices}.')
                    print(f'Exiting\n')
                    exit(2)
                else:
                    indices.append(int(item))
       
        values = []
        for val in args.values:
            for item in val.strip('=[]').split(','):
                if dtype is float:
                    try:
                        values.append(float(item))
                    except ValueError:
                        print(f'\nYou provided an incompatible value for this variable.')
                        print(f'The variable "{args.varname}" is type {dtype}.')
                        print(f'You arguments for --values were {args.values}.')
                        print(f'Exiting\n')
                        exit(2)
                else:
                    values.append(str(item))

        if(not args.silent):
            print(f'\n')
            
        # i is 1-size, not 0 through size
        for i,idx in enumerate(indices):
            if(len(dim_sizes)>1):
                k = np.remainder(idx-1,dim_sizes[1])
                j = int(np.floor(float(idx-1)/float(dim_sizes[1])))
                if(not args.silent):
                    old_val = data['variables'][args.varname]['data'][j][k]
                    print(f'Changing {args.varname}[{j},{k}] from {old_val} to {values[i]}')
                data['variables'][args.varname]['data'][j][k] = values[i]
            else:
                j = idx-1
                if(not args.silent):
                    old_val = data['variables'][args.varname]['data'][j]
                    print(f'Changing {args.varname}[{j}] from {old_val} to {values[i]}')
                data['variables'][args.varname]['data'][j] = values[i]

        if(not args.silent):
            print(f'\n')
            
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

