#!/usr/bin/env python

#
# Original: C. Koven, 2018
# Refactored for json: R. Knox 2025
#
# Script Arguments:
#
# --fin <string>                  : Input filename.  Required
# --fout <string>                 : Output filename. Required if not --overwrite
# --param  <string>               : Name of the parameter to modify or query. Required
# --queryparam <string>             : Report parameter info only, then exit. Not required
# --indices  <int,int,etc>        : List of the parameter indices (1-len) to change (flattens 2d arrays, second dim is the inner).')
# --pft-names <string,string,etc> : List of PFT names to identify indices to change in param
# --values    <value,value,etc>   : List of values to change, must coincide with indices, Required if not queryparam
# --overwrite                     : Use this flag to overwrite the input file
# --silent                        : Suppress output messages
# --listparams                      : List out all the parameters in --fin and exit
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
import write_json

# ========================================================================================
# ========================================================================================
#                                        Main
# ========================================================================================
# ========================================================================================

def main():


    time_str = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    change_log = time_str+': modify_fates_paramfile_json.py'+' '.join(sys.argv[1:])
    
    parser = argparse.ArgumentParser(description='Parse command line arguments to this script.')
    #
    parser.add_argument('--fin', '--input', dest='inputfname', \
                        type=str, help="Input filename.  Required.", required=True)
    parser.add_argument('--fout','--output', dest='outputfname', \
                        type=str, help="Output filename. Required if not --overwrite.")
    parser.add_argument('--param', dest='paramname', type=str, \
                        help="Name of the parameter that you want to modify. Required if not listing.")
    parser.add_argument('--q','--queryparam', dest='queryparam', \
                        help="Report parameter info only, then exit.",action="store_true")
    parser.add_argument('--indices',dest='indices', nargs='+', \
                        help='List of parameter indices (1-len) to change (flattens 2d arrays, second dim is the inner).')
    parser.add_argument('--pft-names',dest='pftnames', nargs='+', \
                        help='List of PFT names to change.')
    parser.add_argument('--values', nargs='+', \
                        help='List of values to change, must coincide with indices')
    parser.add_argument('--overwrite', dest='overwrite', \
                        help="Use this flag to overwrite the input file",action="store_true")
    parser.add_argument('--silent',dest='silent', \
                        help="Suppress output messages",action="store_true")
    parser.add_argument('--listparams',dest='listparams', \
                        help="List out all the parameters and exit",action="store_true")
    
    args = parser.parse_args()
    
    with open(args.inputfname, 'r') as file:
        data = json.load(file)

    # Trivial use case: Just list out parameters and exit
    # -------------------------------------------------------------------------
    
    if(args.listparams):
        print('Listing out parameters names:')
        for key in data['parameters'].keys():
            print(key.strip('"'))
            exit(0)
    else:
        if(args.paramname is None):
            print('You must provide a parameter name if not using --listparams')
            exit(2)
        if(not args.queryparam):
            if(args.values is None):
                print('You must provide values (--values) to overwrite the desired parameter with')
                exit(2)

    # All use cases below here assume that you have a target parameter in mind
    # So lets get some information about that parameter
    # -------------------------------------------------------------------------
    
    try:
        dim_names = data['parameters'][args.paramname]['dims']
    except KeyError:
        print(f'\n\nThe parameter name "{args.paramname}" is not in the dataset. Exiting.\n')

    dim_sizes = list(np.shape(data['parameters'][args.paramname]['data']))
    max_index = np.prod(dim_sizes)

    if(args.queryparam):
        print(f'')
        print(f'Reporting: Parameter: {args.paramname}:')
        print(f'           dimension names: {dim_names}')
        print(f'           dimension sizes: {dim_sizes}')
        print(f'           Current Values:')
        if(len(dim_sizes)>1):
            for i in range(dim_sizes[0]):
                print(f'           {data["parameters"][args.paramname]["data"][i][:]}')
        else:
            print(f'           {data["parameters"][args.paramname]["data"][:]}')
        print(f'           Index Pattern:')
        if(len(dim_sizes)>1):
            for i in range(dim_sizes[0]):
                i0 = i*dim_sizes[1]+1
                i1 = i0+dim_sizes[1]
                print(f'           {list(range(i0,i1))}')
        else:
            print(f'           {list(range(1,dim_sizes[0]+1))}')
        exit(0)

    if ( (args.indices is not None) and (args.pftnames is not None)):
        print(f'You must choose either --indices OR --pft-names when changing values.')
        print(f'You cant choose both. Exiting.')
        exit(2)

    # Everything below here assumes that you are going to actually make modification
    # and not just query things. First step then is to modify the changelog of
    # the file
    # ------------------------------------------------------------------------------

    data['attributes']['history']
    old_hist = data['attributes']['history']
    new_hist = old_hist+'  '+change_log+'.'
    data['attributes']['history'] = new_hist

    
    # Identify the type of the parameter (float,string or integer)
    # -------------------------------------------------------------------------------
    
    dtype_str = data['parameters'][args.paramname]['dtype'].strip()
    if(dtype_str=='float'):
        dtype = type(1.0)
    elif(dtype_str=='integer'):
        dtype = type(1)
    elif(dtype_str=='string'):
        dtype = type('string')
    else:
        print(f'Could not identify the datatype of parameter {args.paramname}')
        print(f'"dtype" should be either "float","integer" or "string"')
        print(f'found: {dtype_str}.  Exiting')
        exit(2)

    # Identify the index positions of the values that you wish to change.
    # There are two mutually exclusive ways to do this, specifying either
    # --indices or --pft-names. In both scenarios we want to generate
    # a list of indices associated with the parameter vector, however
    # in this step we assume the index starts at 1, we will remove
    # 1 position later to be python compliant
    # ----------------------------------------------------------------------

    if (args.indices is not None):
        indices = []
        for index in args.indices:
            for item in index.strip('=[]').split(','):
                if(int(item) > max_index):
                    print(f'\nYou provided an index that is larger than')
                    print(f'the size of the parameter\'s dataset.')
                    print(f'index: {int(item)}, total indices: {int(max_index)}')
                    print(f'You arguments for --indices were {args.indices}.')
                    print(f'Exiting\n')
                    exit(2)
                else:
                    indices.append(int(item))


    if(args.pftnames is not None):

        #Two requirements, must not be a 2d parameter, and must have pft dimension
        if(len(dim_sizes)>1):
            print(f"\nYou specified parameter index changes by pft (--pft-names)")
            print(f' but this parameter "{args.paramname} has two dimensions.')
            print(f' This is incompatible, with pft-name, instead use --idices')
            print(f' and make sure to use help to understand how indices map into 2d')
            print(f' arrays. Exiting')
            exit(2)
        if(not('fates_pft' in data['parameters'][args.paramname]['dims'])):
            print(f"\nYou specified parameter index changes by pft (--pft-names)")
            print(f' but this parameter "{args.paramname} does not have pft as a dimension.')
            print(f' Exiting')
            exit(2)

        indices = []
        cleaned_pft_names = [s.strip() for s in data['parameters']['fates_pftname']['data']]
        #code.interact(local=dict(globals(), **locals()))
        for name in args.pftnames:
            for item in name.strip('=[]').split(','):
                if not(item in data['parameters']['fates_pftname']['data']):
                    print(f"\nYou provided a pft name to filter by, but")
                    print(f" that name was not found in the parameter ")
                    print(f" array fates_pftname")
                    print(f" fates_pftname: {data['parameters']['fates_pftname'].keys()}")
                    print(f" you specified: {item}")
                    print(f" Exiting\n")
                else:
                    # indices uses a 1-numpft convention.
                    indices.append(cleaned_pft_names.index(item)+1)

    # Save the new values in a list, to some simple string
    # manipulations to clean the input values, and use the known
    # data type to store the values correctly (int/float/string)
    # ---------------------------------------------------------------------
                    
    values = []
    for val in args.values:
        for item in val.strip('=[]').split(','):
            if dtype is float:
                try:
                    values.append(float(item))
                except ValueError:
                    print(f'\nYou provided an incompatible value for this parameter.')
                    print(f'The parameter "{args.paramname}" is type {dtype}.')
                    print(f'You arguments for --values were {args.values}.')
                    print(f'Exiting\n')
                    exit(2)
            elif dtype is int:
                values.append(int(item))

            else:
                values.append(str(item))

    if(not args.silent):
        print(f'\n')

    # Apply the changes to the parameter's data vector
    # -----------------------------------------------------------------
    
    for i,idx in enumerate(indices):
        if(len(dim_sizes)>1):
            k = np.remainder(idx-1,dim_sizes[1])
            j = int(np.floor(float(idx-1)/float(dim_sizes[1])))
            if(not args.silent):
                old_val = data['parameters'][args.paramname]['data'][j][k]
                print(f'Changing {args.paramname}[{j},{k}] from {old_val} to {values[i]}')
            data['parameters'][args.paramname]['data'][j][k] = values[i]
        else:
            j = idx-1
            if(not args.silent):
                old_val = data['parameters'][args.paramname]['data'][j]
                print(f'Changing {args.paramname}[{j}] from {old_val} to {values[i]}')
            data['parameters'][args.paramname]['data'][j] = values[i]

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
        #json.dump(data, outfile, indent=2)
        if(not args.silent):
            print(f'Writing to: {output_filename}')

        write_json.traverse_data(outfile,data)

        if(not args.silent):
            print(f'Writing complete')
    
        
# =======================================================================================
# This is the actual call to main

if __name__ == "__main__":
    main()

