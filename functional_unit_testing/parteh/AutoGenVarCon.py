# =============================================================================
# Walk through lines of a file, if a line contains
# the string of interest (EDPftvarcon_inst), then
# parse the string to find the variable name, and save that
# to the list
# =============================================================================

import imp
import code  # For development: code.interact(local=dict(globals(), **locals()))

F90ParamParse = imp.load_source('F90ParamParse','py_modules/F90ParamParse.py')
CDLParse = imp.load_source('CDLParse','py_modules/CDLParse.py')


from F90ParamParse import f90_param_type, GetParamsInFile, GetPFTParmFileSymbols, MakeListUnique
from CDLParse import CDLParseDims, CDLParseParam, cdl_param_type


# -------------------------------------------------------------------------------------
# Check through the fortran Code we are coupling with, determine the list of parameters
# that we need.
# The procedure GetSymbolUsage() returns a list of strings (non-unique)
# -------------------------------------------------------------------------------------

var_list = GetParamsInFile('../../parteh/PRTParametersMod.F90')


# Now look through EDPftvarcon.F90 to determine the variable name in file
# that is associated with the variable pointer

var_list = GetPFTParmFileSymbols(var_list,'../../parteh/PRTParamsFATESMod.F90')

# -------------------------------------------------------------
# We can now cross reference our list of parameters against
# the parameter file. This will create a new list of parameters
# however in the form of a dictionary. This dictionary of
# entries is accessible by its symbol name, and will also
# read in and store the actual parameter values from the file.
# We will use the default file to get the dimensionality.
#
# NOTE: THE CDLPARSE PROCEDURE WILL LOAD IN THE DATA,
# BUT WE DONT NEED IT. THE CDLPARSE PARAM ROUTINE
# IS JUST USED TO GET THE CORRECT DIMENSIONS.  THUS WE
# CAN JUST POINT TO THE DEFAULT CDL FILE IN VERSION CONTROL
#
# -------------------------------------------------------------


default_file_relpath = '../../parameter_files/fates_params_default.cdl'

dims = CDLParseDims(default_file_relpath)

parms = {}
for elem in var_list:
    parms[elem.var_sym] = CDLParseParam(default_file_relpath,cdl_param_type(elem.var_name,True),dims)
    print('Finished loading PFT parameters')



f = open("f90src/UnitWrapMod.F90_in", "r")
contents = f.readlines()
f.close()

# ADD ARGUMENTS TO EDPFTVARCONALLOC
# ---------------------------------

for i,str in enumerate(contents):
    if 'ARGUMENT_IN1' in str:
        index0=i

str=''
icount=0
for key, value in dims.iteritems():
    print('{}'.format(key))
    if(icount==0):
        str+=key
    else:
        str+=(', & \n   '+key)
    icount+=1

strsplit = contents[index0].split('ARGUMENT_IN1')
strreplace = strsplit[0]+str+strsplit[1]

contents[index0] = strreplace


for i,str in enumerate(contents):
    if 'ARGUMENT_DEF1' in str:
        index0=i

str=''
for key, value in dims.iteritems():
    str+=('   integer,intent(in) :: '+key+'\n')


contents[index0] = str





# Identify where we define the variables, and insert the variable definitions

for i,str in enumerate(contents):
    if 'VARIABLE-DEFINITIONS-HERE' in str:
        index0=i

# Identify where we do the pointer assignments, and insert the pointer assignments


for i,str in enumerate(contents):
    if 'POINTER-SPECIFICATION-HERE' in str:
        index0=i

index=index0+2
for symbol, var in parms.iteritems():

    # Generate the dimension names

    dim_alloc_str=''
    icount=0
    for dimname in reversed(var.dim_namelist):
        if(icount==0):
            dim_alloc_str+=dimname
        else:
            dim_alloc_str+=(','+dimname)
        icount+=1


    if(var.ndims==1):
        ins_l1='\t allocate(prt_params%{}({}))\n'.format(symbol,dim_alloc_str)
        ins_l2='\t prt_params%{}(:) = fates_unset_r8\n'.format(symbol)
        ins_l3='\t iv1 = iv1 + 1\n'
        ins_l4='\t prt_params_ptr%var1d(iv1)%var_name = "{}"\n'.format(var.symbol)
        ins_l5='\t prt_params_ptr%var1d(iv1)%var_rp   => prt_params%{}\n'.format(symbol)
        ins_l6='\t prt_params_ptr%var1d(iv1)%vtype    = 1\n'
        ins_l7='\n'
    elif(var.ndims==2):
        ins_l1='\t allocate(prt_params%{}({}))\n'.format(symbol,dim_alloc_str)
        ins_l2='\t prt_params%{}(:,:) = fates_unset_r8\n'.format(symbol)
        ins_l3='\t iv2 = iv2 + 1\n'
        ins_l4='\t prt_params_ptr%var2d(iv2)%var_name = "{}"\n'.format(var.symbol)
        ins_l5='\t prt_params_ptr%var2d(iv2)%var_rp   => prt_params%{}\n'.format(symbol)
        ins_l6='\t prt_params_ptr%var2d(iv2)%vtype    = 1\n'
        ins_l7='\n'
    else:
        print('Auto-generating FORTRAN parameter code does not handle >2D')
        print(symbol)
        print(var.ndims)
        exit(2)

    contents.insert(index,ins_l1)
    contents.insert(index+1,ins_l2)
    contents.insert(index+2,ins_l3)
    contents.insert(index+3,ins_l4)
    contents.insert(index+4,ins_l5)
    contents.insert(index+5,ins_l6)
    contents.insert(index+6,ins_l7)
    index=index+7


f = open("f90src/UnitWrapMod.F90", "w+")
contents = "".join(contents)
f.write(contents)
f.close()
