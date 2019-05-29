

# Walk through lines of a file, if a line contains
# the string of interest (EDPftvarcon_inst), then
# parse the string to find the variable name, and save that
# to the list


class ParamType:

    def __init__(self,var_sym,n_dims):

        self.var_sym  = var_sym
        self.n_dims   = n_dims
        self.var_name = ''




def CheckFile(filename,check_str):
    file_ptr = file(filename)
    var_list = []
    found = False
    for line in file_ptr:
        if check_str in line:
            line_split = line.split()
            #            substr = [i for i in line_split if check_str in i][0]
            substr = line
            p1 = substr.find('%')+1
            if(p1>0):
                substr=substr[p1:]
                p2 = substr.find('(')
                p3 = substr.find(')')
                # Count the number of commas between p2 and p3
                n_dims = substr[p2:p3].count(',')+1
                if(p2>0):
                    var_list.append(ParamType(substr[:p2],n_dims))

    unique_list = []
    for var in var_list:
        found = False
        for uvar in unique_list:
            if (var.var_sym == uvar.var_sym):
                found = True
        if(not found):
            unique_list.append(var)

    return(unique_list)



check_str = 'EDPftvarcon_inst%'
filename  = '../../biogeochem/FatesAllometryMod.F90'

var_list = CheckFile(filename,check_str)


# Add symbols here

var_list.append(ParamType('hgt_min',1))


# Now look through EDPftvarcon.F90 to determine the variable name in file
# that is associated with the variable pointer

filename = '../../main/EDPftvarcon.F90'

f = open(filename,"r")
contents = f.readlines()


var_name_list = []
for var in var_list:
    for i,line in enumerate(contents):
        if (var.var_sym in line) and ('data' in line) and ('=' in line):
            var.var_name = contents[i-2].split()[-1].strip('\'')
            print("{} {} {}".format(var.var_sym,var.var_name,var.n_dims))


f = open("f90src/AllomUnitWrap.F90_in", "r")
contents = f.readlines()
f.close()

# Identify where we define the variables, and insert the variable definitions

for i,str in enumerate(contents):
    if 'VARIABLE-DEFINITIONS-HERE' in str:
        index0=i

index=index0+2
for var in var_list:
    if(var.n_dims==1):
        contents.insert(index,'    real(r8),pointer :: {}(:)\n'.format(var.var_sym))
    elif(var.n_dims==2):
        contents.insert(index,'    real(r8),pointer :: {}(:,:)\n'.format(var.var_sym))
    else:
        print('Incorrect number of dims...')
        exit(-2)
    index=index+1

# Identify where we do the pointer assignments, and insert the pointer assignments


for i,str in enumerate(contents):
    if 'POINTER-SPECIFICATION-HERE' in str:
        index0=i

index=index0+2
for ivar,var in enumerate(var_list):
    if(var.n_dims==1):
        ins_l1='\t allocate(EDPftvarcon_inst%{}(1:numpft))\n'.format(var.var_sym)
        ins_l2='\t EDPftvarcon_inst%{}(:) = nan\n'.format(var.var_sym)
        ins_l3='\t iv1 = iv1 + 1\n'
        ins_l4='\t EDPftvarcon_ptr%var1d(iv1)%var_name = "{}"\n'.format(var.var_name)
        ins_l5='\t EDPftvarcon_ptr%var1d(iv1)%var_rp   => EDPftvarcon_inst%{}\n'.format(var.var_sym)
        ins_l6='\t EDPftvarcon_ptr%var1d(iv1)%vtype    = 1\n'
        ins_l7='\n'
    if(var.n_dims==2):
        ins_l1='\t allocate(EDPftvarcon_inst%{}(1:numpft,1))\n'.format(var.var_sym)
        ins_l2='\t EDPftvarcon_inst%{}(:,:) = nan\n'.format(var.var_sym)
        ins_l3='\t iv2 = iv2 + 1\n'
        ins_l4='\t EDPftvarcon_ptr%var2d(iv2)%var_name = "{}"\n'.format(var.var_name)
        ins_l5='\t EDPftvarcon_ptr%var2d(iv2)%var_rp   => EDPftvarcon_inst%{}\n'.format(var.var_sym)
        ins_l6='\t EDPftvarcon_ptr%var2d(iv2)%vtype    = 1\n'
        ins_l7='\n'

    contents.insert(index,ins_l1)
    contents.insert(index+1,ins_l2)
    contents.insert(index+2,ins_l3)
    contents.insert(index+3,ins_l4)
    contents.insert(index+4,ins_l5)
    contents.insert(index+5,ins_l6)
    contents.insert(index+6,ins_l7)
    index=index+7


f = open("f90src/AllomUnitWrap.F90", "w+")
contents = "".join(contents)
f.write(contents)
f.close()
