from netCDF4 import Dataset
import sys
import os

# program sorts the variables based on the provided list, and pulls them one at a time 
# from an existing file and adds them to a new file in the sorted order.
# input/output based on code here: https://gist.github.com/guziy/8543562

def clobber(filename):
    try:
        print('replacing file: '+filename)
        os.remove(filename)
    except:
        print('file does not exist: '+filename)

### modify the paths below to point to the new and old file names
fnamein = 'fates_params_default.nc'
fnameout = 'fates_params_default_sorted.nc'

# open the input dataset
dsin = Dataset(fnamein)

# make empty lists to hold the variable names in. the first of these is a list of sub-lists, 
# one for each type of variable (based on dimensionality).
# the second is the master list that will contain all variables.
varnames_list = [[],[],[],[],[],[],[],[],[]]
varnames_list_sorted = []

# sort the variables by dimensionality, but mix the PFT x other dimension in with the regular PFT-indexed variables
dimtype_sortorder_dict = {(u'fates_history_height_bins',):0,
 (u'fates_history_size_bins',):1,
 (u'fates_history_age_bins',):2,
 (u'fates_scalar',):3,
 (u'fates_pft', u'fates_string_length'):4,
 (u'fates_pft',):5,
 (u'fates_variants', u'fates_pft'):5,
 (u'fates_hydr_organs', u'fates_pft'):5,
 (u'fates_litterclass',):6,
 (u'fates_NCWD',):7,
 ():8}

# go through each of the variables and assign it to one of the sub-lists based on its dimensionality
for v_name, varin in dsin.variables.iteritems():
    sortorder = dimtype_sortorder_dict[varin.dimensions]
    # if a KeyError, it means that the parameter has a dimension which isn't in dimtype_sortorder_dict. need to add it.
    varnames_list[sortorder].append(v_name)

# go through each of the lists and sort the variable names alphabetically,
# and put them into a master list of all variables.
for i in range(len(varnames_list)):
    varnames_list[i] = sorted(varnames_list[i], key=lambda L: (L.lower(), L))
    varnames_list_sorted.extend(varnames_list[i])

# open the output filename, deleting it if it exists already.
clobber(fnameout)
dsout = Dataset(fnameout,  "w", format="NETCDF3_CLASSIC")

#Copy dimensions
for dname, the_dim in dsin.dimensions.iteritems():
    print dname, len(the_dim)
    dsout.createDimension(dname, len(the_dim) if not the_dim.isunlimited() else None)

print

# go through each variable in the order of the sorted master list, and copy the variable
# as well as all metadata to the new file.
for i in range(len(varnames_list_sorted)):
    v_name = varnames_list_sorted[i]
    varin = dsin.variables[v_name]
    outVar = dsout.createVariable(v_name, varin.datatype, varin.dimensions)
    print v_name
    
    # Copy variable attributes
    outVar.setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})

    # copy data from input file to output file
    outVar[:] = varin[:]

# copy global attributes
dsout.setncatts({k: dsin.getncattr(k) for k in dsin.ncattrs()})

# close the output file
dsout.close()
dsin.close()
