
# =======================================================================================
# This will look through a CDL file for the provided parameter and determine
# the parameter's type, as well as fill an array with the data
# =======================================================================================

import re              # This is a heftier string parser
import code  # For development: code.interact(local=dict(globals(), **locals()))
import numpy as np

# Global identifiers for the type of data
# ---------------------------------------------------------------------------------------

char_type   = 0
int_type    = 1
float_type  = 2
double_type = 3

# If we encounter a "_", ie no data?
no_data_fill='1.e-32'


# This is base object for a parameter
# ===================================
class cdl_param_type:

    def __init__(self,symbol):

        self.datatype = -9
        self.dim_namelist = []
        self.dim_sizelist = []
        self.ndims    = -9
        self.symbol   = symbol
        self.units    = 'NA'

    def Add1DToXD(self,val,indx):

        if(self.ndims==0):
            self.data[indx] = val

        elif(self.ndims==1):
            n1 = self.dim_sizelist[0]
            if((indx<0) or (indx>=n1)):
                print('Problem in CDLParse filling data array')
                print('index must be between {} {}, value = {}'.format(0,n1,indx))
                print('param: {}'.format(self.symbol))
                exit(2)
            else:
                self.data[indx] = val

        elif(self.ndims==2):
            n1 = self.dim_sizelist[0]
            n2 = self.dim_sizelist[1]
            i2 = np.mod(indx,n2)
            i1 = int(indx/n2)
            self.data[i1,i2] = val

        else:
            print('No more than 2 dimensions can be processed by Add1dToXd()')
            exit(2)



# This routine adds a new parameter to the list of cdl_param_types
# ================================================================
def CDLParseParam(file_name,param,dim_dic):

    fp = open(file_name,"r")
    contents = fp.readlines()
    fp.close()

    # Look in the file for the definition for the parameter
    # of interest, note its specified dimensions and cross
    # ref against the dictionary of known dimensions
    # ---------------------------------------------------------
    isfound = False
    for i,line in enumerate(contents):
        if((param.symbol in line) and \
           (not isfound) and \
           (('double' in line) or \
            ('char' in line) or \
            ('float' in line) or \
            ('int' in line))):

            isfound = True

            print('Filling {}'.format(param.symbol))

            datatype = line.split()[0]
            if(datatype.strip()=="float"):
                param.datatype = float_type
            elif(datatype.strip()=="double"):
                param.datatype = double_type
            elif(datatype.strip()=="char"):
                param.datatype = char_type
            elif(datatype.strip()=="int"):
                param.datatype = int_type
            else:
                print('An unknown datatype: {}'.format(datatype.strip()))
                print(' was encountered for parameter: {}'.format(param.symbol))
                exit(2)


            p1=line.find('(')+1
            if(p1>0):
                p2=line.find(')')
                dims_str = line[p1:p2]
                dims_splt = dims_str.split(',')

                for dimname in dims_splt:
                    dimsize = dim_dic.get(dimname.strip())
                    if dimsize:
                        param.dim_namelist.append(dimname.strip())
                        param.dim_sizelist.append(dimsize)
                    else:
                        print('An unknown dimension was requested:')
                        print(' parameter: {}'.format(param.symbol))
                        print(' dimension name: {}'.format(dimname.strip()))
                        exit(2)

                param.ndims = len(param.dim_namelist)

            else:
                param.ndims = 0



            # Allocate and initialize the data space
            if(param.ndims>0):

                param.data = -999*np.ones((param.dim_sizelist))
            else:
                param.data = -999*np.ones((1))




    if(not isfound):
        print('An unknown parameter was requested:')
        print(' parameter: {}'.format(param.symbol))
        exit(2)

    # -----------------------------------------------------------
    # Now that the metadata has been read in, and we
    # know the type of data and its dimensions, lets go retrieve
    # and fill the values in
    # -----------------------------------------------------------


    # First step is to identify the start of the data section:
    # Also, identify the whatever line is next with a ':'
    # ---------------------------------------------------

    iline0=-1
    for i,line in enumerate(contents):
        if('data:' in line):
            iline0 = i
            break

    if(iline0==-1):
        print('Could not find the data section of the CDL file?')
        exit(2)

    # Look for the symbol, again, but now in the "data" section:
    # -----------------------------------------------------------

    isfound = False
    contents=contents[iline0:]
    for i,line in enumerate(contents):

        if(param.symbol in line):

            search_field=True
            lcount=0
            multi_line=''
            while(search_field and (lcount<100)):
                multi_line+=contents[i+lcount]
                if(multi_line.count(';')>0):
                    search_field=False
                else:
                    search_field=True
                lcount=lcount+1

            # Parse the line
            line_split = re.split(',|=',multi_line)
            # Remove the variable name entry
            del line_split[0]

            # This is for real numbers
            if((param.datatype == float_type) or \
               (param.datatype == double_type)):
                ival=0
                indx=0
                for str0 in line_split:
                    str=""
                    isnum=False
                    for s in str0:
                        if (s.isdigit() or s=='.' or s=='-'):
                            str+=s
                            isnum=True
                        elif(s == '_'):
                            str+=no_data_fill
                            isnum=True
                    if(isnum):
                        param.Add1DToXD(float(str),indx)
                        indx=indx+1
                    else:
                        print('No-data values encountered during parameter read in')
                        print('for parameter {}'.format(param.symbol))
                        print('bad value: {}'.format(str0))
                        print('data: {}'.format(line_split))
                        exit(2)

            # This is a string
            #            elif(param.datatype == 1):
            #                for str0 in line_split:
            #                    # Loop several times to trim stuff off
            #                    for i in range(5):
            #                        str0=str0.strip().strip('\"').strip(';').strip()
            #                        param.vals.append(str0)


#            if(param.symbol == 'fates_hydr_thetas_node'):




    return(param)





# This routine returns a dictionary with dimension names and sizes
# ====================================================================

def CDLParseDims(file_name):

    fp = open(file_name,"r")
    contents = fp.readlines()
    fp.close()

    # Identify the line with the "dimensions:" tag
    # Also, identify the whatever line is next with a ':'
    # ---------------------------------------------------

    iline0=-1
    for i,line in enumerate(contents):
        if('dimensions:' in line):
            iline0 = i
            break

    if(iline0==-1):
        print("The CDL Parser could not find the dimensions section")
        print(" in your output file")
        print(" exiting...")
        exit(2)

    iline1=-1
    for i,line in enumerate(contents):
        if((':' in line) and \
           (i > iline0) and \
           ('"' not in line)):
            iline1 = i
            break

    if(iline1==-1):
        print("The CDL Parser could not find a section")
        print(" following the dimensions section.")
        print(" exiting...")
        exit(2)

    # Loop between the two and save the dimensions
    # --------------------------------------------

    dim_dic = {}
    for i in range(iline0+1,iline1):

        # If there is an equals sign, then there is data
        if('=' in contents[i]):

            # Split string into chunks
            sline = contents[i].split()
            dim_dic[sline[0]] = int(sline[2])


    if(len(dim_dic)==0):
        print("No valid dimensions found in your CDL file")
        print(" exiting...")
        exit(2)

    return(dim_dic)
