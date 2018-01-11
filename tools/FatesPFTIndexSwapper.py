# =======================================================================================
# =======================================================================================

import numpy as np
import sys
import getopt
import code  # For development: code.interact(local=locals())
from datetime import datetime
from matplotlib.dates import date2num, num2date
import csv
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
     print(' python CloneFATESPFTFile.py -h --num-pfts=>n> ')
     print('                                --pft-index=<integer position> ')
     print('                                --fin=<netcdf-file-in> ')
     print('                                --fout=<netcdf-file-out>')
     print('')
     print('')
     print(' -h --help ')
     print('     print this help message')
     print('')
     print('')
     print(' --pft-index=<integer position>')
     print('     This is the PFT index of the base file that you want copied into the new file')
     print('')
     print('')
     print(' --num-pfts=<n>')
     print('     This is the desired number of pfts you want in the output file.')
     print('')
     print('')
     print(' --fin=<netcdf-file-in>')
     print('     This is the full path to the netcdf file you are basing off of')
     print('')
     print('')
     print(' --fout=<netcdf-file-out>')
     print('     This is the full path to the netcdf file you are writing to.')
     print('')
     print('')
     print('=======================================================================')


def interp_args(argv):

    argv.pop(0)  # The script itself is the first argument, forget it

    # Name of the conversion file

    input_fname = "none"
    output_fname = "none"
    donor_pft_index = -9
    num_pft_out = -9
    try:
        opts, args = getopt.getopt(argv, 'h',["fin=","fout=","pft-index=","num-pfts="])

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
        elif o in ("--pft-index"):
            donor_pft_index = int(a.strip())
        elif o in ("--num-pfts"):
            num_pft_out = int(a.strip())
        else:
            assert False, "unhandled option"


    if (input_fname == "none"):
        print("You must specify an input file:")
        usage()
        sys.exit(2)

    if (output_fname == "none"):
        print("You must specify an output file:")
        usage()
        sys.exit(2)    

    if (donor_pft_index == -9):
        print("You must specify the donor pft index, > 0:")
        usage()
        sys.exit(2)

    if (num_pft_out == -9):
        print("You must specify the number of output pfts")
        usage()
        sys.exit(2)


    return (input_fname,output_fname,donor_pft_index,num_pft_out)


# ========================================================================================
# ========================================================================================
#                                        Main
# ========================================================================================
# ========================================================================================

def main(argv):

    # Interpret the arguments to the script
    [input_fname,output_fname,donor_pft_index,num_pft_out] = interp_args(argv)
    

    # Open the netcdf files
    fp_out = netcdf.netcdf_file(output_fname, 'w')
    
    fp_in  = netcdf.netcdf_file(input_fname, 'r')

#    code.interact(local=locals())

    for key, value in sorted(fp_in.dimensions.iteritems()):
        print('Creating Dimension: ',value)
        if(key==pft_dim_name):
            fp_out.createDimension(key,int(num_pft_out))
        else:
            fp_out.createDimension(key,int(value))

    for key, value in sorted(fp_in.variables.iteritems()):
        print('Creating Variable: ',key)
        #   code.interact(local=locals())
        
        out_var = fp_out.createVariable(key,'f',(fp_in.variables.get(key).dimensions))
        in_var  = fp_in.variables.get(key)
        out_var.units     = in_var.units
        out_var.long_name = in_var.long_name

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
            out_var.assignValue(float(fp_in.variables.get(key).data))
        elif(pft_dim_found==-1):
            out_var[:] = in_var[:]
        elif( (pft_dim_found==0) & (pft_dim_len==1) ):           # 1D fates_pft
            tmp_out = fp_in.variables.get(key).data[donor_pft_index-1] * np.ones([num_pft_out])
            out_var[:] = tmp_out
        elif( (pft_dim_found==1) & (pft_dim_len==2) ):           # 2D hdyro_organ - fate_pft
            dim2_len = fp_in.dimensions.get(fp_in.variables.get(key).dimensions[0])
            tmp_out  = np.ones([dim2_len,num_pft_out])
            for idim in range(0,dim2_len):
                tmp_out[idim,:] = tmp_out[idim,:] * fp_in.variables.get(key).data[idim,donor_pft_index-1]
            out_var[:] = tmp_out
        else:
            print('This variable has a dimensioning that we have not considered yet.')
            print('Please add this condition to the logic above this statement.')
            print('Aborting')
            exit(2)

    fp_out.history = "This file was made from CloneHLMPFTFile.py"


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








