#!/usr/bin/env python

#### this script modifies the default FATES parameter file to generate
#    a file used in testing E3SM
#    Parser code was based off of modify_fates_paramfile.py

import os
import argparse
import code  # For development: code.interact(local=dict(globals(), **locals()))
import json
import datetime
import sys

# 1. Get the directory of the currently executing script
# os.path.dirname(__file__) gets the directory of main_script.py
script_dir = os.path.dirname(os.path.abspath(__file__))

# 2. Insert this directory at the beginning of the Python search path
# This tells Python to look in the same folder as main_script.py for imports
if script_dir not in sys.path:
    sys.path.insert(0, script_dir)

import write_json


def main():


    time_str = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    change_log = time_str+': modify_fates_paramfile_json.py'+' '.join(sys.argv[1:])
    
    parser = argparse.ArgumentParser(description='Parse command line arguments to this script.')

    parser.add_argument('--fin', '--input', dest='inputfname', \
                        type=str, help="Input filename.  Required.", required=True)
    parser.add_argument('--silent',dest='silent', \
                        help="Suppress output messages",action="store_true")

    args = parser.parse_args()
    
    with open(args.inputfname, 'r') as file:
        batch_data = json.load(file)

    
    base_file = batch_data['base_file'].strip()
    new_file  = batch_data['new_file'].strip()

    # This should be a list of integers
    pft_trim_list = ",".join(str(i) for i in batch_data['pft_trim_list'])
    
    swapcmd="../tools/pft_index_swapper.py --pft-indices="+pft_trim_list+" --fin="+base_file+" --fout="+new_file
    os.system(swapcmd)
        
    with open(new_file, 'r') as file:
        base_data = json.load(file)

    if "pft_parameters" in batch_data['parameters']:
        for pfts,pft_dic in batch_data['parameters']['pft_parameters'].items():
            ipfts = [int(x) for x in pfts.split(',')]
            for param_name,data_list in pft_dic.items():
                # 2D list
                if(isinstance(data_list[0], list)):
                    dim2 = len(data_list)
                    for index,i in enumerate(ipfts):
                        ipft = i-1 #because python uses 0 and our file uses 1 indexing
                        for id2 in range(dim2):
                            base_data['parameters'][param_name]['data'][id2][ipft] = data_list[id2][index]
                # 1D list
                else:
                    for index,i in enumerate(ipfts):
                        ipft = i-1 #because python uses 0 and our file uses 1 indexing
                        base_data['parameters'][param_name]['data'][ipft] = data_list[index]
                        
    if "non_pft_parameters" in batch_data['parameters']:
        for param_name,data_list in batch_data['parameters']['non_pft_parameters'].items():
            base_data['parameters'][param_name]['data'] = data_list

    with open(new_file, 'w') as outfile:
        write_json.traverse_data(outfile,base_data)
    
        
if __name__ == "__main__":
    main()
