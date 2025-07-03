#!/usr/bin/env python

# =======================================================================================
# This script modifies any FATES parameter file to update it to a new API spec.
#    It can change variable names
#    It can add new variables
#    It can add attributes
#    It can update attributes
#    It can add new dimensions
# =======================================================================================

import os
import argparse
import code  # For development: code.interact(local=dict(globals(), **locals()))
import xml.etree.ElementTree as et
import numpy as np

# Newer versions of scipy have dropped the netcdf module and
# netcdf functions are part of the io parent module
try:
    from scipy import io as nc

except ImportError:
    from scipy.io import netcdf as nc
	 
# =======================================================================================

def load_xml(xmlfile): 

    # This routine parses the XML tree

    xmlroot = et.parse(xmlfile).getroot()
    print("\nOpened: {}\n".format(xmlfile))

    base_cdl = xmlroot.find('base_file').text
    new_cdl = xmlroot.find('new_file').text

    pft_list = xmlroot.find('pft_list').text.replace(" ","")
    
    modroot = xmlroot.find('mods') 

    return(base_cdl,new_cdl,pft_list,modroot)

# =======================================================================================

def str2fvec(numstr):

    # Convert a list of strings into floating point numbers
    
    numvec = [float(i) for i in numstr.split(',')]
    return(numvec)

# =======================================================================================

def str2ivec(numstr):

    # Convert a list of strings into integer numbers
    
    intvec = [int(i) for i in numstr.split(',')]
    return(intvec)

# =======================================================================================

def createvar(ncfile,paramname,dimnames,units,longname,dcode,sel_values):

    # Create a new netcdf variable inside an existing netcdf dataset (append)
    ncvar = ncfile.createVariable(paramname,dcode,dimnames)
    ncvar.units = units
    ncvar.long_name = longname
    if( not dimnames):
        ncvar.assignValue(sel_values)
    else:
        ncvar[:] = sel_values
    ncfile.flush()
    
    return(ncfile,ncvar)

# =======================================================================================

def selectvalues(ncfile,dimnames,ipft_list,values,dtype):

    # Reduce a list of values so that onlythe chosen pft values are left. This
    # only works on float arrays currently.  We need to pass in a file
    # so that we can get the dimension sizes associated with the dimension names.

    if(len(ipft_list) != ncfile.dimensions['fates_pft']):
        print('you list of pfts in the xml file must be')
        print('the same size as the fates_pft dimension')
        print('in your destination file. exiting')
        print('len(ipft_list) = {}'.format(len(ipft_list)))
        print('fates_pft dim = {}'.format(ncfile.dimensions['fates_pft']))
        exit(2)

    # shift the pft list to a base of 0 instead of 1
        
    pft_dim = -1
    dim_size = [0 for i in range(0,len(dimnames))]

    for idim, name in enumerate(dimnames):
        dim_size[idim] = ncfile.dimensions[name]
        if(name=='fates_pft'):
            pft_dim = idim
            pft_dim_size = ncfile.dimensions['fates_pft']
    
    if(len(dimnames) == 1):

        if(pft_dim==0):
            dim1_list = ipft_list
        else:
            dim1_list = range(0,dim_size[0])
        
        sel_values = np.zeros([len(dim1_list)])
        for i,ipft in enumerate(dim1_list):
            sel_values[i] = values[ipft]
            
    elif(len(dimnames) == 2 ):

        if(dtype=="c"):

            if(pft_dim>0):
                print("problem with pft_dim: {},{}".format(dimnames[0],dimnames[1]));exit(2)
            
            if(pft_dim==0):
                dim1_list = ipft_list
            else:
                dim1_list = range(0,dim_size[0])

            sel_values = np.chararray([len(dim1_list),dim_size[1]])
            sel_values[:] = ""
            for i,ipft in enumerate(dim1_list):
                for j,cc in enumerate(values[ipft]):
                    sel_values[i][j] = cc

                
        elif(dtype=="d"):

            if(pft_dim==0):
                print("problem with pft_dim: {},{}".format(dimnames[0],dimnames[1]));exit(2)
                
            if(pft_dim==1):
                dim1_list = ipft_list
            else:
                dim1_list = range(0,dim_size[1])

            sel_values = np.zeros([dim_size[0],len(dim1_list)])
            for i,ipft in enumerate(dim1_list):
                for i2 in range(0,dim_size[0]):
                    id = i2*len(dim1_list)+ipft
                    sel_values[i2,i] = values[id]

    else:
                    
        # Scalar
        #code.interact(local=dict(globals(), **locals()))
        sel_values = float(values[0])

    return(sel_values)

# =======================================================================================

def removevar(base_nc,varname):

    # Remove a variable from a dataset. This is actually the hardest thing to do!
    # The trick here, is to copy the whole file, minus the variable of interest
    # into a temp file. Then completely remove the old file, and 

    fp_base = nc.netcdf_file(base_nc, 'r',mmap=False)

    new_nc = os.popen('mktemp').read().rstrip('\n')
    fp_new  = nc.netcdf_file(new_nc, 'w',mmap=False)

    found = False
    for key, value in sorted(fp_base.dimensions.items()):
        if( key == varname ):
            found = True
        else:
            fp_new.createDimension(key,int(value))

    for key, value in fp_base.variables.items():

        if(key == varname):
            found = True
        else:
            datatype = value.typecode()
            new_var = fp_new.createVariable(key,datatype,value.dimensions)
            if(value.data.size == 1):
                new_var.assignValue(float(value.data))
            else:
                new_var[:] = value[:].copy()
                
            new_var.units = value.units
            new_var.long_name = value.long_name
            #try:
            #    new_var.use_case = value.use_case
            #except:
            #    new_var.use_case = "undefined"
    
    fp_new.history = fp_base.history

    if(not found):
        print("was not able to find variable: ()".format(varname))
        exit(2)

    fp_new.flush()
    fp_base.close()
    fp_new.close()
    
    mvcmd = "(rm -f "+base_nc+";mv "+new_nc+" "+base_nc+")"
    os.system(mvcmd)
    
    
# =======================================================================================

def main():

    # Parse arguments
    parser = argparse.ArgumentParser(description='Parse command line arguments to this script.')
    parser.add_argument('--f', dest='xmlfile', type=str, help="XML control file  Required.", required=True)
    args = parser.parse_args()


    # Load the xml file, which contains the base cdl, the output cdl,
    #  and the parameters to be modified
    [base_cdl,new_cdl,pft_list,modroot] = load_xml(args.xmlfile)

    ipft_list = str2ivec(pft_list)
    
    for i,ipft in enumerate(ipft_list):
        ipft_list[i] = ipft_list[i]-1
    
    # Convert the base cdl file into a temp nc binary
    base_nc = os.popen('mktemp').read().rstrip('\n')
    gencmd = "ncgen -o "+base_nc+" "+base_cdl
    os.system(gencmd)
    modlist = []
    for mod in modroot:
        if(not('type' in mod.attrib.keys())):
            print("mod tag must have attribute type")
            print("exiting")
            exit(2)

        if(mod.attrib['type'].strip() == 'dimension_add'):

            try:
                dimname = mod.find('di').text.strip()
            except:
                print("{}, no dimension (di), exiting".format(mod.attrib['type']));exit(2)

            try:
                values = str2ivec(mod.find('val').text.strip())
            except:
                print("no values (val), exiting");exit(2)
                
            if(len(values)>1):
                print("The dimension size should be a scalar")
                exit(2)

            ncfile = nc.netcdf_file(base_nc,"a",mmap=False)
            ncfile.createDimension(dimname, values[0])
            ncfile.flush()
            ncfile.close()

            print("dimension: {}, size: {}, added".format(dimname,values[0]))
            
        elif(mod.attrib['type'].strip() == 'dimension_del'):

            try:
                dimname = mod.find('di').text.strip()
            except:
                print('define the dimension name to delete using a <na> tag')
                exit(2)

            # Find which parameters use this dimension
            ncfile = nc.netcdf_file(base_nc,"r",mmap=False)
            found = False
            for key, value in sorted(ncfile.dimensions.items()):
                if(key==dimname):
                    found=True

            if(not found):
                print("could not find {} for deletion".format(dimname));exit(2)
            
            for key, value in sorted(ncfile.variables.items()):
                hasdim = any([dim == dimname for dim in list(value.dimensions) ])
                
                if (hasdim):
                    print("parameter: {}, removed (to accomodate dimension removal)".format(key))
                    removevar(base_nc,key)
                    
            ncfile.close()
            
            removevar(base_nc,dimname)
            print("dimension: {}, removed".format(dimname))
            
        elif(mod.attrib['type'].strip() == 'variable_add'):
            print("Adding Parameter")
            try:
                paramname = mod.find('na').text.strip()
            except:
                print("no name (na), exiting");exit(2)

            #try:
            #    dtype = mod.find('dt').text.strip()
            #except:
            #    print("no data type (dt), exiting");exit(2)
                
            try:
                # print("trying dimnames: {}".format(paramname))
                dimnames = tuple(mod.find('di').text.replace(" ","").split(','))
            except:
                print("no data type (di), exiting");exit(2)
                
            try:
                # print("trying units: {}".format(paramname))
                units = mod.find('un').text.strip()
            except:
                print("no units (un), exiting");exit(2)
                
            try:
                # print("trying ln: {}".format(paramname))
                longname = mod.find('ln').text.strip()
            except:
                print("no long-name (ln), exiting");exit(2)

            ncfile = nc.netcdf_file(base_nc,"a",mmap=False)

            try:
                # print("trying val: {}".format(paramname))
                valstring = mod.find('val').text.strip()
                values = str2fvec(valstring)
            except Exception as emsg:
                # print("type: {}".format(type(valstring)))
                if(isinstance(valstring,type(None))):
                    print("Warning: no values (val). Setting undefined (i.e. '_'): {}\n".format(paramname))
                    sel_values = ncfile.variables['fates_dev_arbitrary_pft'].data
                    dcode = "d"
                elif(isinstance(valstring,str)):
                    dcode = "c"
                    values = valstring.split(',')
                    for i,val in enumerate(values):
                        values[i] = val.strip()
                        print("value: {},{}".format(i,values[i]))

                    sel_values = selectvalues(ncfile,list(dimnames),ipft_list,values,dcode)
                else:
                    print("exception, unknown values (val), exiting: {}".format(emsg));exit(2)
                    #print("no values (val), exiting");exit(2)
            else:
            #code.interact(local=dict(globals(), **locals()))
            
                if(dimnames[0]=='scalar' or dimnames[0]=='none' or dimnames[0]==''):
                    dimnames = ()
                elif(isinstance(values[0],float)):
                    dcode = "d"
                else:
                    print("Unknown value type: {} {}".format(type(values[0]),paramname));exit(2)

                sel_values = selectvalues(ncfile,list(dimnames),ipft_list,values,dcode)

            [ncfile,ncvar] = createvar(ncfile,paramname,dimnames,units,longname,dcode,sel_values)
            ncfile.flush()
            ncfile.close()

            print("parameter: {}, added".format(paramname))

            
        elif(mod.attrib['type'] == 'variable_del'):
            print("Deleting Parameter")
            try:
                paramname = mod.find('na').text.strip()
            except:
                print('must define the parameter name to delete, using name attribute')
                exit(2)
            removevar(base_nc,paramname)
            print("parameter: {}, removed".format(paramname))

            
        elif(mod.attrib['type'] == 'variable_change'):  

            print("Changing Parameter")
            
            try:
                paramname_o = mod.attrib['name'].strip()
            except:
                print("to change a parameter, the field must have a name attribute")
                exit(2)
                
            ncfile = nc.netcdf_file(base_nc,"a",mmap=False)
            ncvar_o = ncfile.variables[paramname_o]
            # dims_o  = ncvar_o.dimensions
            dtype_o = ncvar_o.typecode()
            units_o = ncvar_o.units.decode("utf-8")
            longname_o = ncvar_o.long_name.decode("utf-8")

            # Check for a parameter name change
            try:
                newparamname = mod.find('na').text.strip()
            except:
                newparamname = None

            # Check for a dimensionality change
            try:
                dimnames = tuple(mod.find('di').text.replace(" ","").split(','))
            except:
                dimnames = None

            # Change the parameter's name and/or the dimensionality
            if(not isinstance(newparamname,type(None)) or not isinstance(dimnames,type(None))):

                print("Changing Name")
                
                # Initialize the parameter name to pass to the create variable function
                # If this is None, a dimension update is happening and this will be updated
                # below
                paramname = newparamname

                # If no dimension change, use the previous dimensions
                if not dimnames:
                    dims_o = ncvar_o.dimensions
                else:
                    # If dimension to change to is scalar, set dims_o to an empty list
                    if(dimnames[0]=='scalar' or dimnames[0]=='none' or dimnames[0]==''):
                        dims_o = ()
                    else:
                        dims_o = dimnames

                # If there is no dimension change grab the original data
                if (isinstance(dimnames,type(None))):
                    if not dims_o:
                        ncvardata = float(ncvar_o.data)
                        # [ncfile,ncvar] = createvar(ncfile,paramname,dims_o,units_o,longname_o,dtype_o,float(ncvar_o.data))
                    else:
                        ncvardata = ncvar_o[:].copy()
                        # [ncfile,ncvar] = createvar(ncfile,paramname,dims_o,units_o,longname_o,dtype_o,ncvar_o[:].copy())
                # If there is a dimension change create zero filled data
                else:
                    # If there isn't a parameter name change, grab the old name
                    if (isinstance(paramname,type(None))):
                        paramname = paramname_o

                    # If changing to scalar
                    if not dims_o:
                        ncvardata = 0.
                    # read the dimensions and create a
                    else:
                        dimsize = []
                        for idim,name in enumerate(dimnames):
                            dimsize.append(ncfile.dimensions[name])
                        ncvardata = np.zeros(dimsize)

                # Create the new variable
                [ncfile,ncvar] = createvar(ncfile,paramname,dims_o,units_o,longname_o,dtype_o,ncvardata)

            else:
                ncvar = ncvar_o
                dims_o = ncvar_o.dimensions
                
            # Change the metadata:
            try:
                units = mod.find('un').text.strip()
            except:
                units = None
            if(not isinstance(units,type(None))):
                ncvar.units = units 
                
            try:
                longname = mod.find('ln').text.strip()
            except:
                longname = None
            if(not isinstance(longname,type(None))):
                ncvar.long_name = longname
                
            try:
                values = str2fvec(mod.find('val').text.strip())
            except:
                values = None
                
            if(not isinstance(values,type(None))):
                sel_values = selectvalues(ncfile,list(dims_o),ipft_list,values,dtype_o)
                    
                # Scalars have their own thing
                if(ncvar.data.size == 1):
                    ncvar.assignValue(sel_values)
                else:
                    ncvar[:] = sel_values[:]

                
            ncfile.flush()
            ncfile.close()
            
            # Finally, if we did perform a re-name, and
            # created a new variable. We need to delete the old one
            # This is not necessary for a dimension only change
            if(not isinstance(newparamname,type(None))):
               removevar(base_nc,paramname_o)
               paramname = paramname_o
               
            print("parameter: {}, modified".format(paramname_o))

               
    # Sort the new file
    new_nc = os.popen('mktemp').read().rstrip('\n')
    os.system("../tools/ncvarsort.py --silent --fin "+base_nc+" --fout "+new_nc+" --overwrite")

    # Dump the new file to the cdl
    os.system("ncdump "+new_nc+" > "+new_cdl)

    print("\nAPI update complete, see file: {}\n".format(new_cdl))
    
        
# This is the actual call to main

if __name__ == "__main__":
    main()
