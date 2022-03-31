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
import numpy as np

# ========================================================================================
# ========================================================================================
#                                        Main
# ========================================================================================
# ========================================================================================

def main():
    parser = argparse.ArgumentParser(description='Parse command line arguments to this script.')
    #
    parser.add_argument('--var','--variable', dest='varname', type=str, help="What variable to modify? Required.", required=True)
    parser.add_argument('--pft','--PFT', dest='pftnum', type=int, help="PFT number to modify. If this and --pftname are missing and --allPFTs is not specified, will assume a global variable.")
    parser.add_argument('--pftname', '--PFTname', dest='pftname', type=str, help="Name of PFT to modify.  alternate argument to --pft or --allpfts")
    parser.add_argument('--allPFTs', '--allpfts', dest='allpfts', help="apply to all PFT indices. Cannot use at same time as --pft argument.", action="store_true")
    parser.add_argument('--fin', '--input', dest='inputfname', type=str, help="Input filename.  Required.", required=True)
    parser.add_argument('--fout','--output', dest='outputfname', type=str, help="Output filename.  Required.", required=True)
    parser.add_argument('--val', '--value', dest='val', type=str, help="New value of PFT variable.  Must be interpretable as either a real number or a comma-separated list of real numbers. Required.", required=True)
    parser.add_argument('--O','--overwrite', dest='overwrite', help="If present, automatically overwrite the output file.", action="store_true")
    parser.add_argument('--silent', '--s', dest='silent', help="prevent writing of output.", action="store_true")
    parser.add_argument('--nohist', dest='nohist', help="prevent recording of the edit in the history attribute of the output file", action="store_true")
    parser.add_argument('--changeshape', dest='changeshape', help="allow script to change shape of specified variable, and all other variables with the relevant dimension, if necessary", action="store_true")
    #
    args = parser.parse_args()
    #
    # work with the file in some random temporary place so that if something goes wrong, then nothing happens to original file and it doesn't make a persistent output file
    tempdir = tempfile.mkdtemp()
    tempfilename = os.path.join(tempdir, 'temp_fates_param_file.nc')
    ncfile_old = None
    rename_pft = False
    #
    try:
        outputval = float(args.val)
        if args.changeshape:
            raise Exception
    except:
        try:
            #print('output variable not interpretable as real. trying array')
            outputval = np.fromstring(args.val, sep=',', dtype=np.float32)
            if len(outputval) == 0:
                raise RuntimeError('output variable needs to have size greater than zero')
        except:
            if args.varname != 'fates_pftname':
                raise RuntimeError('output variable not interpretable as real or array')
            else:
                rename_pft = True
    #
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
        if ndim_file > 2:
            raise ValueError('variable dimensionality is too high for this script')
        for i in range(ndim_file):
            if var.dimensions[i] == 'fates_pft':
                ispftvar = True
                npft_file = var.shape[i]
                pftdim = i
                otherdimpresent = False
            elif var.dimensions[i] in ['fates_history_age_bins','fates_history_size_bins','fates_history_coage_bins','fates_history_height_bins','fates_NCWD','fates_litterclass','fates_leafage_class','fates_prt_organs','fates_hydr_organs','fates_variants']:
                otherdimpresent = True
                otherdimname = var.dimensions[i]
                otherdimlength = var.shape[i]
            elif var.dimensions[i] == 'fates_string_length' and rename_pft:
                otherdimpresent = True
                otherdimname = var.dimensions[i]
                otherdimlength = var.shape[i]
            else:
                raise ValueError('variable is not on either the PFT or scalar dimension')
        #
        if args.changeshape:
            ### if we are allowing the script to change the shape of the variable, then we need to figure out if that's really a thing that needs to happen.
            ### first identify what dimension we would change the shape of if we had to.
            length_specified = len(outputval)
            if length_specified != otherdimlength:
                ### ok, we find ourselves in the situation where we need to rewrite the netcdf from scratch with its revised shape.
                #
                # first lets chech to make sure the dimension we are changing can be changed without breaking things.
                plastic_dimensions_list = ['fates_history_age_bins','fates_history_size_bins','fates_history_coage_bins','fates_history_height_bins','fates_leafage_class']
                if otherdimname not in plastic_dimensions_list:
                    raise ValueError('asking to change the shape of a dimension, '+otherdimname+', that will probably break things')
                else:
                    print('WARNING: we need to change the dimension of '+otherdimname)
                ### close the file that's open and start over.
                ncfile.close()
                os.remove(tempfilename)
                ncfile = nc.netcdf_file(tempfilename, 'w')
                ncfile_old = nc.netcdf_file(args.inputfname, 'r')
                #
                try:
                    ncfile.history = ncfile_old.history
                except:
                    print('no history')
                #
                ### copy over and, when needed, modify the dimensions
                for name, dimlength in ncfile_old.dimensions.items():
                    #print(name, dimlength)
                    if name != otherdimname:
                        ncfile.createDimension(name, dimlength)
                    else:
                        ncfile.createDimension(name, length_specified)
                        #print(name, length_specified)
                #
                ### copy over and, when needed, modify the variables
                for name, variable in ncfile_old.variables.items():
                    variabledims = variable.dimensions
                    #print(name, variabledims)
                    x = ncfile.createVariable(name, variable.data.dtype, variable.dimensions)
                    try:
                        x.units = variable.units
                    except:
                        print('no units')
                    try:
                        x.long_name = variable.long_name
                    except:
                        print('no long name')
                    #
                    if len(variable.dimensions) > 0:
                        if not otherdimname in variable.dimensions:
                            x[:] = variable[:]
                        else:
                            if len(variable.dimensions) == 1:
                                if length_specified > otherdimlength:
                                    print('WARNING: Variable '+name+' has a dimension that has been reshaped. New length is longer than old, so its been filled in with zeros.')
                                    x[0:otherdimlength] = variable[0:otherdimlength]
                                    x[otherdimlength:length_specified] = 0
                                else:
                                    print('WARNING: Variable '+name+' has a dimension that has been reshaped. New length is shorter than old, so its been truncated.')
                                    x[0:length_specified] = variable[0:length_specified]
                            elif len(variable.dimensions) == 2:
                                if length_specified > otherdimlength:
                                    print('WARNING: Variable '+name+' has a dimension that has been reshaped. New length is longer than old, so its been filled in with zeros.')
                                    x[0:otherdimlength,:] = variable[0:otherdimlength,:]
                                    x[otherdimlength:length_specified,:] = 0
                                else:
                                    print('WARNING: Variable '+name+' has a dimension that has been reshaped. New length is shorter than old, so its been truncated.')
                                    x[0:length_specified,:] = variable[0:length_specified,:]
                    else:
                        x.assignValue(float(variable.data))
                #
                var = ncfile.variables[args.varname]
            else:
                # declare as none for now
                ncfile_old = None
        #
        if (args.pftnum == None and args.pftname == None and ispftvar) and not args.allpfts:
            raise ValueError('pft value is missing but variable has pft dimension.')
        if (args.pftnum != None or args.pftname != None) and args.allpfts:
            raise ValueError("can't specify both a PFT number and the argument allPFTs.")
        if (args.pftnum != None or args.pftname != None) and not ispftvar:
            raise ValueError('pft value is present but variable does not have pft dimension.')
        if (args.pftnum != None and args.pftname != None):
            raise ValueError('can only specify pft number or name, not both.')
        if (args.pftnum == None or args.pftname != None) and not args.allpfts and ispftvar:
            ## now we need to figure out what the number of the pft that has been given a name argument
            pftnamelist = []
            npftnames = ncfile.variables['fates_pftname'].shape[0]
            for i in range(npftnames):
                pftname_bytelist = list(ncfile.variables['fates_pftname'][i,:])
                pftname_stringlist = [i.decode('utf-8') for i in pftname_bytelist]
                pftnamelist.append(''.join(pftname_stringlist).strip())
            n_times_pft_listed = pftnamelist.count(args.pftname.strip())
            if n_times_pft_listed != 1:
                raise ValueError('can only index by PFT name if the chosen PFT name occurs once and only once.')
            pftnum = pftnamelist.index(args.pftname.strip())
            args.pftnum=pftnum +1
        if args.pftnum != None and ispftvar:
            if not rename_pft:
                if args.pftnum > npft_file:
                    raise ValueError('PFT specified ('+str(args.pftnum)+') is larger than the number of PFTs in the file ('+str(npft_file)+').')
                if pftdim == 0:
                    if not args.silent:
                        print('replacing prior value of variable '+args.varname+', for PFT '+str(args.pftnum)+', which was '+str(var[args.pftnum-1])+', with new value of '+str(outputval))
                    var[args.pftnum-1] = outputval
                if pftdim == 1:
                    if not args.silent:
                        print('replacing prior value of variable '+args.varname+', for PFT '+str(args.pftnum)+', which was '+str(var[:,args.pftnum-1])+', with new value of '+str(outputval))
                    var[:,args.pftnum-1] = outputval
            else:
                pftname_in_bytelist = list(ncfile.variables['fates_pftname'][args.pftnum-1,:])
                pftname_in_stringlist = [i.decode('utf-8') for i in pftname_in_bytelist]
                print('replacing prior value of pft name for PFT '+str(args.pftnum)+', which was "'+''.join(pftname_in_stringlist).strip()+'", with new value of "'+args.val+'"')
                var[args.pftnum-1] = args.val.ljust(otherdimlength)
        elif args.allpfts and ispftvar:
            if pftdim == 0:
                if not args.silent:
                    print('replacing prior values of variable '+args.varname+', for all PFTs, which were '+str(var[:])+', with new value of '+str(outputval))
                var[:] = outputval
            if pftdim == 1:
                if not args.silent:
                    print('replacing prior values of variable '+args.varname+', for all PFTs, which were '+str(var[:])+', with new value of '+str(outputval))
                var[:] = outputval
        elif args.pftnum == None and not ispftvar and ndim_file > 0:
            if not otherdimpresent:
                if not args.silent:
                    print('replacing prior value of variable '+args.varname+', which was '+str(var[:])+', with new value of '+str(outputval))
                var[:] = outputval
            else:
                #print(var.shape)
                #print(outputval.shape)
                if not args.silent:
                    print('replacing prior value of variable '+args.varname+', which was '+str(var[:])+', with new value of '+str(outputval))
                var[:] = outputval
        elif ndim_file < 1:
            if not args.silent:
                print('replacing prior value of scalar variable '+args.varname+', which was '+str(var.data)+', with new value of '+str(outputval))
            var.assignValue(outputval)
        else:
            raise ValueError('Nothing happened somehow.')
        #
        if not args.nohist:
            # write to the netcdf file history attribute what you just did.
            actionstring = 'modify_fates_paramfile.py '+' '.join(sys.argv[1:])
            timestampstring = datetime.datetime.fromtimestamp(time.time()).strftime('%a %b %d %Y, %H:%M:%S')
            #
            oldhiststr = ncfile.history.decode('utf-8')
            newhiststr = oldhiststr + "\n "+timestampstring + ': ' + actionstring
            ncfile.history = newhiststr
        #
        ncfile.close()
        if type(ncfile_old) != type(None):
            ncfile_old.close()
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

