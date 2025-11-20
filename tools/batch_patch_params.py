#!/usr/bin/env python

#### this script modifies the default FATES parameter file to generate
#    a file used in testing E3SM
#    Parser code was based off of modify_fates_paramfile.py

import os
import argparse
import code  # For development: code.interact(local=dict(globals(), **locals()))
import json
import write_json
import datetime
import sys

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
    pft_trim_list = ", ".join(str(i) for i in batch_data['pft_trim_list'])
    
    #code.interact(local=dict(globals(), **locals()))
    
    swapcmd="../tools/pft_index_swapper.py --pft-indices="+pft_trim_list+" --fin="+base_file+" --fout="+new_file
    os.system(swapcmd)
        
    with open(new_file, 'r') as file:
        base_data = json.load(file)

    for pft,pft_dic in batch_data['parameters']['pft_parameters'].items():
        ipft = int(pft)-1
        for param_name,data_list in pft_dic.items():
            if(len(data_list)==1):
                base_data['parameters'][param_name]['data'] = data_list
            else:
                for i,val in enumerate(data_list):
                    #code.interact(local=dict(globals(), **locals()))
                    base_data['parameters'][param_name]['data'][i][ipft] = val

    for param_name,data_list in batch_data['parameters']['non_pft_parameters'].items():
        base_data['parameters'][param_name]['data'] = data_list

    with open(new_file, 'w') as outfile:
        write_json.traverse_data(outfile,base_data)
    
        
if __name__ == "__main__":
    main()
