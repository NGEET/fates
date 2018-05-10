#!/usr/bin/env python

# =======================================================================================
#
# This python script will open an input FATES parameter file, and given a list of PFT
# indices supplied by the user, will create a new parameter file with PFTs entries cloned
# from the original file as-per the list of indices supplied by the user.
#
# First Added, Ryan Knox: Thu Jan 11 13:36:14 PST 2018
# =======================================================================================

import numpy as np
import sys
import getopt
import code  # For development: code.interact(local=locals())
from datetime import datetime
from scipy.io import netcdf
import matplotlib.pyplot as plt

# =======================================================================================
# Parameters
# =======================================================================================

pft_dim_name = 'fates_pft'


class timetype:   
    
    # This is time, like the thing that always goes forward and cant be seen
    # or touched, insert creative riddle here

    def __init__(self,ntimes):
    
        self.year  = -9*np.ones((ntimes))
        self.month = -9*np.ones((ntimes))
        # This is a floating point decimal day
        self.day   = -9.0*np.ones((ntimes))

        # This is a decimal datenumber
        self.datenum = -9.0*np.ones((ntimes))


def usage():
     print('')
     print('=======================================================================')
     print('')
     print(' python FatesPFTIndexSwapper.py -h --pft-indices <integer position> ')
     print('                                --fin <netcdf-file-in> ')
     print('                                --fout <netcdf-file-out>')
     print('')
     print('')
     print(' -h --help ')
     print('     print this help message')
     print('')
     print('')
     print(' --pft-indices <integer positions ie 1,1,2,3,5,7>')
     print('     This is a comma delimited list of integer positions of the PFTs')
     print('     to be copied into the new file. Note that first pft position')
     print('     is treated as 1 (not C or python like), and any order or multiples')
     print('     of indices can be chosen')
     print('')
     print('')
     print(' --fin <netcdf-file-in>')
     print('     This is the full path to the netcdf file you are basing off of')
     print('')
     print('')
     print(' --fout <netcdf-file-out>')
     print('     This is the full path to the netcdf file you are writing to.')
     print('')
     print('')
     print('=======================================================================')


def interp_args(argv):

    argv.pop(0)  # The script itself is the first argument, forget it

    # Name of the conversion file

    input_fname = "none"
    output_fname = "none"
    donor_pft_indices = -9
    donot_pft_indices_str = ''
    try:
        opts, args = getopt.getopt(argv, 'h',["fin=","fout=","pft-indices="])
    except getopt.GetoptError as err:
        print('Argument error, see usage')
        usage()
        sys.exit(2)

    for o, a in opts:
        if o in ("-h", "--help"):
            usage()
            sys.exit(0)
        elif o in ("--fin"):
            input_fname = a
        elif o in ("--fout"):
            output_fname = a
        elif o in ("--pft-indices"):
            donor_pft_indices_str = a.strip()
        else:
            assert False, "unhandled option"


    if (input_fname == "none"):
        print("You must specify an input file:\n\n")
        usage()
        sys.exit(2)

    if (output_fname == "none"):
        print("You must specify an output file:\n\n")
        usage()
        sys.exit(2)    

    if (donor_pft_indices_str == ''):
        print("You must specify at least one donor pft index!\n\n")
        usage()
        sys.exit(2)
    else:
        donor_pft_indices = []
        for strpft in donor_pft_indices_str.split(','):
            donor_pft_indices.append(int(strpft))        


    return (input_fname,output_fname,donor_pft_indices)


# ========================================================================================
# ========================================================================================
#                                        Main
# ========================================================================================
# ========================================================================================

def main(argv):

    # Interpret the arguments to the script
    [input_fname,output_fname,donor_pft_indices] = interp_args(argv)

    num_pft_out = len(donor_pft_indices)

    # Open the netcdf files
    fp_out = netcdf.netcdf_file(output_fname, 'w')
    
    fp_in  = netcdf.netcdf_file(input_fname, 'r')

    for key, value in sorted(fp_in.dimensions.iteritems()):
        if(key==pft_dim_name):
            fp_out.createDimension(key,int(num_pft_out))
            print('Creating Dimension: {}={}'.format(key,num_pft_out))
        else:
            fp_out.createDimension(key,int(value))
            print('Creating Dimension: {}={}'.format(key,value))

    for key, value in sorted(fp_in.variables.iteritems()):
        print('Creating Variable: ',key)
        #   code.interact(local=locals())
        
        
        in_var  = fp_in.variables.get(key)
        

        # Idenfity if this variable has pft dimension
        pft_dim_found = -1
        pft_dim_len   = len(fp_in.variables.get(key).dimensions)

        for idim, name in enumerate(fp_in.variables.get(key).dimensions):
            # Manipulate data 
            if(name==pft_dim_name):
                pft_dim_found = idim
            
        # Copy over the input data
        # Tedious, but I have to permute through all combinations of dimension position
        if( pft_dim_len == 0 ):
            out_var = fp_out.createVariable(key,'f',(fp_in.variables.get(key).dimensions))
            out_var.assignValue(float(fp_in.variables.get(key).data))
        elif(pft_dim_found==-1):
            out_var = fp_out.createVariable(key,'f',(fp_in.variables.get(key).dimensions))
            out_var[:] = in_var[:]
        elif( (pft_dim_found==0) & (pft_dim_len==1) ):           # 1D fates_pft
            out_var = fp_out.createVariable(key,'f',(fp_in.variables.get(key).dimensions))
            tmp_out  = np.zeros([num_pft_out])
            for id,ipft in enumerate(donor_pft_indices):
                tmp_out[id] = fp_in.variables.get(key).data[ipft-1]
            out_var[:] = tmp_out


        elif( (pft_dim_found==1) & (pft_dim_len==2) ):           # 2D hdyro_organ - fate_pft
            out_var = fp_out.createVariable(key,'f',(fp_in.variables.get(key).dimensions))
            dim2_len = fp_in.dimensions.get(fp_in.variables.get(key).dimensions[0])
            tmp_out  = np.zeros([dim2_len,num_pft_out])
            for id,ipft in enumerate(donor_pft_indices):
                for idim in range(0,dim2_len):
                    tmp_out[idim,id] = fp_in.variables.get(key).data[idim,ipft-1]
            out_var[:] = tmp_out

        elif( (pft_dim_found==0) & (pft_dim_len==2) ):          # fates_pft - string_length
            out_var = fp_out.createVariable(key,'c',(fp_in.variables.get(key).dimensions))
            dim2_len = fp_in.dimensions.get(fp_in.variables.get(key).dimensions[1])
            out_var[:] = np.empty([num_pft_out,dim2_len], dtype="S{}".format(dim2_len))
            for id,ipft in enumerate(donor_pft_indices):
                out_var[id] = fp_in.variables.get(key).data[ipft-1]

        else:
            print('This variable has a dimensioning that we have not considered yet.')
            print('Please add this condition to the logic above this statement.')
            print('Aborting')
            for idim, name in enumerate(fp_in.variables.get(key).dimensions):
               print("idim: {}, name: {}".format(idim,name))
            exit(2)

        out_var.units     = in_var.units
        out_var.long_name = in_var.long_name

    fp_out.history = "This file was made from FatesPFTIndexSwapper.py \n Input File = {} \n Indices = {}"\
          .format(input_fname,donor_pft_indices)

    #var_out.mode = var.mode
    #fp.flush()

    fp_in.close()
    fp_out.close()

    print('Cloneing complete!')
    exit(0)




# =======================================================================================
# This is the actual call to main
   
if __name__ == "__main__":
    main(sys.argv)








