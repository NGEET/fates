#!/usr/bin/env python

#### this script sorts a FATES parameter file. It accepts the following flags
# --input or --fin: input filename.
# --output or --fout: output filename.  If missing, will assume its directly modifying the input file, and will prompt unless -O is specified

import sys
import os
import argparse
import datetime
import json
import code  # For development: code.interact(local=dict(globals(), **locals()))

# 1. Get the directory of the currently executing script
# os.path.dirname(__file__) gets the directory of main_script.py
script_dir = os.path.dirname(os.path.abspath(__file__))

# 2. Insert this directory at the beginning of the Python search path
# This tells Python to look in the same folder as main_script.py for imports
if script_dir not in sys.path:
    sys.path.insert(0, script_dir)

import write_json

# program sorts the variables based on the provided list, and pulls them one at a time
# from an existing file and adds them to a new file in the sorted order.
# input/output based on code here: https://gist.github.com/guziy/8543562

def main():
    parser = argparse.ArgumentParser(description='Parse command line arguments to this script.')
    #
    parser.add_argument('--fin', '--input', dest='fnamein', type=str, help="Input filename.  Required.", required=True)
    parser.add_argument('--fout','--output', dest='fnameout', type=str, help="Output filename.  Required.", required=True)
    parser.add_argument('--O','--overwrite', dest='overwrite', help="If present, automatically overwrite the output file.", action="store_true")
    parser.add_argument('--debug', dest='debug', help="If present, output more diagnostics", action="store_true")
    parser.add_argument('--silent', dest='silent', help="If present, prevents printing messages", action="store_true")
    #
    args = parser.parse_args()

    with open(args.fnamein, 'r') as file:
        data = json.load(file)

    # Note that this is a recursive dictionary
    # data['dimensions'] is a dictionary with string:int key pairs
    # data['parameters'][<key>]['dims'] generates a list of keys

    sortorder_list = [
        ['fates_history_age_bins'],
        ['fates_history_coage_bins'],
        ['fates_history_height_bins'],
        ['fates_history_size_bins'],
        ['fates_history_damage_bins'],
        ['fates_hydr_organs'],
        ['fates_prt_organs'],
        ['fates_plant_organs'],
        ['fates_pft', 'fates_string_length'],
        ['fates_hydr_organs', 'fates_string_length'],
        ['fates_prt_organs', 'fates_string_length'],
        ['fates_plant_organs', 'fates_string_length'],
        ['fates_litterclass', 'fates_string_length'],
        ['fates_landuseclass', 'fates_string_length'],
        ['fates_pft'],
        ['fates_hydr_organs', 'fates_pft'],
        ['fates_leafage_class', 'fates_pft'],
        ['fates_prt_organs', 'fates_pft'],
        ['fates_plant_organs', 'fates_pft'],
        ['fates_hlm_pftno', 'fates_pft'],
        ['fates_litterclass'],
        ['fates_NCWD'],
        ['fates_landuseclass'],
        ['fates_landuseclass', 'fates_pft'],
        ['scalar']]

    # Check to make sure that the dimensions defined in
    # all parameters can be found in the sorting list
    
    for param_name,param in data['parameters'].items():
        found_dim = False
        for dim in sortorder_list:
            if (param['dims'] == dim):
                found_dim = True
        if(not found_dim):
            print(f'Could not find the dimensions for parameter: {param_name}')
            print(f' in the set of defined dimensions to control order.')
            print(f' Exiting.')
            exit(2)


    # Construct a new dictionary, based on the original
    # With everything but the parameter data
    
    new_data = {}
    new_data['attributes'] = data['attributes']
    new_data['parameters'] = {}
    new_data['dimensions'] = {}

    # Add a log to the history
    time_str = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    change_log = time_str+': sort_parameters.py'+' '.join(sys.argv[1:])
    old_hist = data['attributes']['history']
    new_hist = old_hist+'  '+change_log+'.'
    new_data['attributes']['history'] = new_hist
    
    # Lets check to make sure all these dimensions
    # are actually used, if they are not used by any parameters
    # then we cull it from the new list
    for dim_name,dim_size in data['dimensions'].items():
        found_dim = False
        for param_name,param in data['parameters'].items():
            if (dim_name in param['dims']):
                found_dim = True
        new_data['dimensions'][dim_name] = dim_size

    # Now lets loop through our list of dimension groups
    # And for each of these groups, sort in alphabetical key order
    for dimgroup in sortorder_list:
        sublist = {}
        for param_name,param in data['parameters'].items():
            if (param['dims'] == dimgroup):
                sublist[param_name] = param
        for param_name in sorted(sublist):
            new_data['parameters'][param_name] = sublist[param_name]

    # Final Check, the old and new dictionaries should have the
    # same number of parameters in them...
    if(len(new_data['parameters']) != len(data['parameters'])):
        print(f'The new dictionary does not have the same number')
        print(f' of parameters as the original dictionary.')
        print(f' Exiting')
        exit(2)
    
            
    if args.fnameout is not None:
        output_filename = args.fnameout
        if(args.overwrite is True):
            print(f'\nYou have specified both an output file, and to overwrite.')
            print(f'You can\'t do both, you must pick one.')
            print(f'see: --output and --overwrite')
            print(f'Exiting\n')
            exit(2)
    else:
        if(args.overwrite is True):
            output_filename = args.fnamein

    with open(output_filename, 'w') as outfile:
    
        if(not args.silent):
            print(f'Writing to: {output_filename}')

        write_json.traverse_data(outfile,new_data)

        if(not args.silent):
            print(f'Writing complete')

# Call main
if __name__ == "__main__":
    main()
