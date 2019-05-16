#!/usr/bin/env python

#### this script sorts a FATES parameter file. It accepts the following flags
# --input or --fin: input filename.
# --output or --fout: output filename.  If missing, will assume its directly modifying the input file, and will prompt unless -O is specified

import netCDF4 as nc
import sys
import os
import argparse

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
    #
    args = parser.parse_args()
    #
    # open the input dataset
    dsin = nc.Dataset(args.fnamein, 'r')
    #
    # make empty lists to hold the variable names in. the first of these is a list of sub-lists,
    # one for each type of variable (based on dimensionality).
    # the second is the master list that will contain all variables.
    varnames_list = [[],[],[],[],[],[],[],[],[],[]]
    varnames_list_sorted = []
    #
    # sort the variables by dimensionality, but mix the PFT x other dimension in with the regular PFT-indexed variables
    dimtype_sortorder_dict = {
    (u'fates_history_age_bins',):0,
    (u'fates_history_height_bins',):1,
    (u'fates_history_size_bins',):2,
    (u'fates_pft', u'fates_string_length'):3,
    (u'fates_prt_organs', u'fates_string_length'):4,
    (u'fates_pft',):5,
    (u'fates_variants', u'fates_pft'):5,
    (u'fates_hydr_organs', u'fates_pft'):5,
    (u'fates_leafage_class', u'fates_pft'):5,
    (u'fates_prt_organs', u'fates_pft'):5,
    (u'fates_litterclass',):6,
    (u'fates_NCWD',):7,
    ():8}
    #
    # go through each of the variables and assign it to one of the sub-lists based on its dimensionality
    for v_name, varin in dsin.variables.iteritems():
        sortorder = dimtype_sortorder_dict[varin.dimensions]
        # if a KeyError, it means that the parameter has a dimension which isn't in dimtype_sortorder_dict. need to add it.
        varnames_list[sortorder].append(v_name)
    #
    # go through each of the lists and sort the variable names alphabetically,
    # and put them into a master list of all variables.
    for i in range(len(varnames_list)):
        varnames_list[i] = sorted(varnames_list[i], key=lambda L: (L.lower(), L))
        varnames_list_sorted.extend(varnames_list[i])
    #
    # write list of variables in ourput order
    if args.debug:
        print(varnames_list_sorted)
    #
    # open the output filename, deleting it if it exists already.
    if os.path.isfile(args.fnameout):
        if args.fnameout == args.fnamein:
            raise ValueError('Error: output file name is the same as the input file name.')
        elif args.overwrite:
            print('replacing file: '+args.fnameout)
            os.remove(args.fnameout)
        else:
            raise ValueError('Output file already exists and overwrite flag not specified for filename: '+args.fnameout)
    #
    dsout = nc.Dataset(args.fnameout,  "w")
    #
    #Copy dimensions
    for dname, the_dim in dsin.dimensions.iteritems():
        print dname, the_dim.size
        dsout.createDimension(dname, the_dim.size )
    #
    print
    #
    try:
        dsout.history = dsin.history
    except:
        print('no history!')
    #
    #
    # go through each variable in the order of the sorted master list, and copy the variable
    # as well as all metadata to the new file.
    for i in range(len(varnames_list_sorted)):
        v_name = varnames_list_sorted[i]
        varin = dsin.variables[v_name]
        outVar = dsout.createVariable(v_name, varin.datatype, varin.dimensions)
        print v_name
        #
        outVar.setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})
        outVar[:] = varin[:]
        #
        # copy global attributes
        dsout.setncatts({k: dsin.getncattr(k) for k in dsin.ncattrs()})#
    #
    # close the output file
    dsin.close()
    dsout.close()

# =======================================================================================
# This is the actual call to main

if __name__ == "__main__":
    main()
