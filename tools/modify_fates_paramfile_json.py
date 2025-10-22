#!/usr/bin/env python

#### this script modifies a FATES parameter file. It accepts the following flags
# --var or --variable: variable.
# --pft or --PFT: PFT number. If this is missing, script will assume that its a global variable that is being modified.
# --input or --fin: input filename.
# --output or --fout: output filename.  If missing, will assume its directly modifying the input file, and will prompt unless -O is specified
# --O or --overwrite: overwrite output file without asking.
# --value or --val: value to put in variable
# --s or --silent: don't write anything on successful execution.
####
#
# Written by C. Koven, 2018
#

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


def comma_separated_list(input_string):
    """Splits a string by commas and returns a list of stripped strings."""
    # Use list comprehension to split and strip whitespace from each item
    return [item.strip() for item in input_string.split(',')]

# ========================================================================================
# ========================================================================================
#                                        Main
# ========================================================================================
# ========================================================================================

def main():

    parser = argparse.ArgumentParser(description='Parse command line arguments to this script.')
    #
    parser.add_argument('--fin', '--input', dest='inputfname', type=str, help="Input filename.  Required.", required=True)
    parser.add_argument('--fout','--output', dest='outputfname', type=str, help="Output filename.  Required.", required=True)
    parser.add_argument('--var','--variable', dest='varname', type=str, help="What variable to modify? Required.", required=True)
    parser.add_argument('--q','--queryvar', dest='queryvar', help="Report variable info only, then exit.",action="store_true")
    parser.add_argument('--indices', nargs='+', help='List of variable indices (1-len) to change (flattens 2d arrays, second dim is the inner).')
    parser.add_argument('--values', nargs='+', help='List of values to change, must coincide with indices')
    parser.add_argument('--overwrite', dest='overwrite', help="Use this flag to overwrite the input file",action="store_true")
     
    args = parser.parse_args()
    
    # work with the file in some random temporary place so that if something goes wrong,
    # then nothing happens to original file and it doesn't make a persistent output file
    tempdir = tempfile.mkdtemp()
    tempfilename = os.path.join(tempdir, 'temp_fates_param_file.json')

    
    

    with open(args.inputfname, 'r') as file:
        # json.load() reads the file object and converts the JSON data 
        # (like arrays, objects, strings) into corresponding Python objects (lists, dictionaries, strings).
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
            print(f'           dimension names: {dim_sizes}')
            print(f'           data: {vardata}')
            exit(0)

        indices = []
        for index in args.indices:
            for item in index.split(','):
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
            for item in val.split(','):
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

        # i is 1-size, not 0 through size
        for i,idx in enumerate(indices):
            if(len(dim_sizes)>1):
                k = np.remainder(idx-1,dim_sizes[1])
                j = int(np.floor(float(idx-1)/float(dim_sizes[1])))
                data['variables'][args.varname]['data'][j][k] = values[i]
            else:
                j = idx-1
                data['variables'][args.varname]['data'][j] = values[i]

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


        
                
    # data is a nested dictionary
    # >>> data['dimensions']['fates_NCWD']
    #     4
    #>>> type(data['dimensions']['fates_NCWD'])
    #    <class 'int'>
    #>>> data['variables']['fates_q10_froz']
    #    {'dims': ['scalar'], 'long_name': 'Q10 for frozen-soil respiration rates', 'units': 'unitless', 'data': [1.5]}


# =======================================================================================
# This is the actual call to main

if __name__ == "__main__":
    main()

