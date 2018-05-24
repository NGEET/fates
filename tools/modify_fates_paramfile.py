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
from scipy.io import netcdf as nc
import argparse
import shutil
import tempfile
import sys
import datetime
import time

# ========================================================================================
# ========================================================================================
#                                        Main
# ========================================================================================
# ========================================================================================

def main():
    parser = argparse.ArgumentParser(description='Parse command line arguments to this script.')
    #
    parser.add_argument('--var','--variable', dest='varname', type=str, help="What variable to modify? Required.", required=True)
    parser.add_argument('--pft','--PFT', dest='pftnum', type=int, help="PFT number to modify. If this is missing and --allPFTs is not specified, will assume a global variable.")
    parser.add_argument('--allPFTs', '--allpfts', dest='allpfts', help="apply to all PFT indices. Cannot use at same time as --pft argument.", action="store_true")
    parser.add_argument('--fin', '--input', dest='inputfname', type=str, help="Input filename.  Required.", required=True)
    parser.add_argument('--fout','--output', dest='outputfname', type=str, help="Output filename.  Required.", required=True)
    parser.add_argument('--val', '--value', dest='val', type=float, help="New value of PFT variable.  Required.", required=True)
    parser.add_argument('--O','--overwrite', dest='overwrite', help="If present, automatically overwrite the output file.", action="store_true")
    parser.add_argument('--silent', '--s', dest='silent', help="prevent writing of output.", action="store_true")
    parser.add_argument('--nohist', dest='nohist', help="prevent recording of the edit in the history attribute of the output file", action="store_true")
    #
    args = parser.parse_args()
    #
    # work with the file in some random temporary place so that if something goes wrong, then nothing happens to original file and it doesn't make a persistent output file
    tempdir = tempfile.mkdtemp()
    tempfilename = os.path.join(tempdir, 'temp_fates_param_file.nc')
    #
    try:
        shutil.copyfile(args.inputfname, tempfilename)
        #
        ncfile = nc.netcdf_file(tempfilename, 'a')
        #
        var = ncfile.variables[args.varname]
        #
        ### check to make sure that, if a PFT is specified, the variable has a PFT dimension, and if not, then it doesn't. and also that shape is reasonable.
        ndim_file = len(var.dimensions)
        ispftvar = False
        # for purposes of current state of this script, assume 1D 
        if ndim_file > 1:
            raise ValueError('variable dimensionality is too high for this script')
        if ndim_file < 1:
            raise ValueError('variable dimensionality is too low for this script. FATES assumes even scalars have a 1-length dimension')
        for i in range(ndim_file):
            if var.dimensions[i] == 'fates_pft':
                ispftvar = True
                npft_file = var.shape[i]
                pftdim = 0
            elif var.dimensions[i] == 'fates_scalar':
                npft_file = None
                pftdim = None
            else:
                raise ValueError('variable is not on either the PFT or scalar dimension')
        if (args.pftnum == None and ispftvar) and not args.allpfts:
            raise ValueError('pft value is missing but variable has pft dimension.')
        if (args.pftnum != None) and args.allpfts:
            raise ValueError("can't specify both a PFT number and the argument allPFTs.")
        if args.pftnum != None and not ispftvar:
            raise ValueError('pft value is present but variable does not have pft dimension.')
        if args.pftnum != None and ispftvar:
            if args.pftnum > npft_file:
                raise ValueError('PFT specified ('+str(args.pftnum)+') is larger than the number of PFTs in the file ('+str(npft_file)+').')
            if pftdim == 0:
                if not args.silent:
                    print('replacing prior value of variable '+args.varname+', for PFT '+str(args.pftnum)+', which was '+str(var[args.pftnum-1])+', with new value of '+str(args.val))
                var[args.pftnum-1] = args.val
        elif args.allpfts and ispftvar:
            if pftdim == 0:
                if not args.silent:
                    print('replacing prior values of variable '+args.varname+', for all PFTs, which were '+str(var[:])+', with new value of '+str(args.val))
                var[:] = args.val            
        elif args.pftnum == None and not ispftvar:
            if not args.silent:
                print('replacing prior value of variable '+args.varname+', which was '+str(var[:])+', with new value of '+str(args.val))
            var[:] = args.val
        else:
            raise ValueError('Nothing happened somehow.')
        #
        if not args.nohist:
            # write to the netcdf file history attribute what you just did.
            actionstring = 'modify_fates_paramfile.py '+' '.join(sys.argv[1:])
            timestampstring = datetime.datetime.fromtimestamp(time.time()).strftime('%a %b %d %Y, %H:%M:%S')
            #
            oldhiststr = ncfile.history
            newhiststr = oldhiststr + "\n "+timestampstring + ': ' + actionstring
            ncfile.history = newhiststr
        #
        ncfile.close()
        #
        #
        # now move file from temporary location to final location
        #
        # check to see if output file exists
        if os.path.isfile(args.outputfname):
            if args.overwrite:
                if not args.silent:
                    print('replacing file: '+args.outputfname)
                os.remove(args.outputfname)
            else:
                raise ValueError('Output file already exists and overwrite flag not specified for filename: '+args.outputfname)
        #
        shutil.move(tempfilename, args.outputfname)
        shutil.rmtree(tempdir, ignore_errors=True)
    except:
        shutil.rmtree(tempdir, ignore_errors=True)
        raise


# =======================================================================================
# This is the actual call to main
   
if __name__ == "__main__":
    main()

