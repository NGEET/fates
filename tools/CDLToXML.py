
# =======================================================================================
# This will look through a CDL file for the provided parameter and determine
# the parameter's type, as well as fill an array with the data
# =======================================================================================

import re              # This is a heftier string parser
import code  # For development: code.interact(local=dict(globals(), **locals()))
import numpy as np
import sys
import argparse
from datetime import date
import xml.etree.ElementTree as ET

debug = True

# This is saved as a dictionary with
# symbol names as keys
class param_type:

    # vartypes = ["float","double","char","integer"]
    
    def __init__(self,symbol,dim_names,vartype,dims):
        self.symbol = symbol  # somewhat redundant, also the key
        self.dim_names = dim_names
        self.vartype = vartype
        self.meta = {}

        # allocate data based on dimensions
        
        if(dim_names[0]=='scalar'):
            self.data = np.zeros(1)

        else:
            dshape = []
            for dname in dim_names:
                dshape.append(dims[dname])
            if(vartype=='float'):
                self.data = np.zeros(dshape,dtype=np.float32)
            if(vartype=='double'):
                self.data = np.zeros(dshape,dtype=np.float64)
            if(vartype=='integer'):
                self.data = np.zeros(dshape,dtype=np.int32)
            if(vartype=='char'):
                # create empty array to initialize
                # The last dimension should be "fates_string_length"
                if(dim_names[-1]!='fates_string_length'):
                    print('Problem parsing char parameter: {}'.format(symbol))
                    print('This method assumes the last dimension is fates_string_length')
                    exit(2)
                nstr = dims[dim_names[0]]
                if len(dim_names)>2:
                    nstr = nstr * dims[dim_names[1]]
                elif len(dim_names)>3:
                    print('Three dimensional string arrays are not handled')
                    exit(2)
                stra = []
                for i in range(nstr):
                    stra.append(' ')
                self.data = np.array(stra,dtype=object)
        
    def add_meta(self,key,val):
        self.meta[key] = val

    def add_data(self,strdata):
        
        # Check if the arrays are the same total size
        if( len(strdata) != np.prod(self.data.shape)):
            print("Inconsistency between parameter allocation size")
            print("and the amount of data found in the file")
            print("parameter: {}, allocated: {}, found: {}".format(self.symbol,np.sum(self.data.shape),len(strdata)))
            print("dims: {}".format(self.dim_names))
            print("strdata: {}".format(strdata))
            exit(2)

        if(self.vartype=='float' or self.vartype=='double'):
            adata = []
            for i,str in enumerate(strdata):
                if '_' in str:
                    adata.append(np.nan)
                else:
                    adata.append(float(str))

        elif(self.vartype=='integer'):
            adata = []
            for i,str in enumerate(strdata):
                if '_' in str:
                    adata.append(-99999)
                else:
                    adata.append(int(str))
        else:
            adata = strdata
            
        ii=0
        if(len(self.data.shape)==1):
            for i in range(self.data.shape[0]):
                self.data[i] = adata[ii]
                ii=ii+1
                
        elif(len(self.data.shape)==2):
            for i in range(self.data.shape[0]):
                for j in range(self.data.shape[1]):
                    self.data[i,j] = adata[ii]
                    ii=ii+1
                        
        elif(len(self.data.shape)==3):
            for i in range(self.data.shape[0]):
                for j in range(self.data.shape[1]):
                    for k in range(self.data.shape[2]):
                        self.data[i,j,k] = adata[ii]
                        ii=ii+1

                        
    def print_var(self):
        print('\nVariable: {}\n  type: {}'.format(self.symbol,self.vartype))
        for key, val in self.meta.items():
            print("  {} = {}".format(key,val))
        for val in self.dim_names:
            print("  dim: {}".format(val))
        
def main(argv):

    parser = argparse.ArgumentParser(description='Parse command line arguments to this script.')
    parser.add_argument('--cdlfile', dest='cdlfile', type=str, help="cdl file path, required.", required=True)
    parser.add_argument('--outfile', dest='outfile', type=str, help="xml+json file path, no suffix! required.", required=True)
    parser.set_defaults(feature=False)
    parser.add_argument('--verbose', dest='verbose', action='store_true')
   
    
    args = parser.parse_args()

    params,dims = CDLParse(args.cdlfile,args.verbose)

    # Load database of data types (float,int,string)
    # for each variable
    dtypes = {}

    dtype_str = ['float','integer','string']
    
    with open('datatypes.txt', 'r') as file:
        for line in file:
            linevec = line.split(' ')
            dtypes[linevec[0].strip()] = int(linevec[1].strip())
            #print(f'{dtypes[linevec[0].strip()]}')

    xmlfile = args.outfile+'.xml'
    if(True):
        with open(xmlfile,"w") as file:
            file.write('<xml version="1.0" encoding="us-ascii" >\n')
            file.write('  <history>{} : {}</history>\n'.format(\
                date.today().strftime("%d/%m/%y"), \
                'First instantation, copied from: {}.'.format(args.cdlfile)))
            file.write('  <dimensions>\n')
            for key, val in dims.items():
                file.write('    <dim name="{}"> {} </dim>\n'.format(key,val))
            file.write('  </dimensions>\n')
            file.write('  <parameters>\n')
            for key, val in params.items():
                file.write('    <par name="{}">\n'.format(key))
                file.write('      <dtype> {} </dtype>\n'.format(dtype_str[dtypes[key]-1]))
                file.write('      <dims> {} </dims>\n'.format(', '.join(val.dim_names)))
                file.write('      <long> {} </long>\n'.format(val.meta['long_name']))
                file.write('      <units> {} </units>\n'.format(val.meta['units']))
                if(len(val.data.shape)>1):
                    data_strs = ''
                    for k in range(val.data.shape[0]):
                        data_strs = data_strs + ', '.join([str(c) for c in val.data[k,:]])
                        if(k<val.data.shape[0]-1):
                            data_strs = data_strs + ', '

                    file.write('      <data> {} </data>\n'.format(data_strs))
                else:
                    data_strs = [str(c) for c in val.data]
                    file.write('      <data> {} </data>\n'.format(', '.join(data_strs)))
                #if(len(val.dim_names)>1):
                #code.interact(local=dict(globals(), **locals()))
                
                file.write('    </par>\n')
            file.write('  </parameters>\n')
            file.write('</xml>\n')

    jsonfile = args.outfile+'.json'
    if(True):

        ob = "{{{}}"
        cb = "{{}}}"
        
        with open(jsonfile,"w") as file:
            file.write('{\n')
            file.write('  "attributes": {')
            file.write('"history": "{}, {} "'.format(\
                date.today().strftime("%d/%m/%y"), \
                'First instantation, copied from: {}.'.format(args.cdlfile)))
            file.write('},\n')
            file.write('  "dimensions": ')
            file.write('{\n')
            it=0
            for key, val in dims.items():
                it=it+1
                if(key.strip() != 'fates_string_length'):
                    file.write('    "{}": {}'.format(key,val))
                    if(it<len(dims.items())-1):
                        file.write(',\n')
                    else:
                        file.write('\n')
                        
            file.write('  },\n')
            file.write('  "variables": {\n')
            item_list = list(params.items())
            num_items = len(item_list)
            
            for i,(key, val) in enumerate(item_list):
                #print('---{}---'.format(key))
                file.write('    "{}":'.format(key))
                file.write(' {\n')
                file.write('      "dtype": "{}",\n'.format(dtype_str[dtypes[key]-1]))
                joinlist = [f'"{name.strip()}"' for name in val.dim_names if name != "fates_string_length"]
                output_string = ', '.join(joinlist)
                file.write('      "dims": [{}],\n'.format(output_string))
                file.write('      "long_name": "{}",\n'.format(val.meta['long_name']))
                file.write('      "units": "{}",\n'.format(val.meta['units']))

                if(len(val.data.shape)>1):
                    data_strs = '['
                    for k in range(val.data.shape[0]):
                        row_strs = [str(c).strip() for c in val.data[k,:]]
                        new_strs = []
                        for c in row_strs:
                            if(c=='nan'):
                                new_strs.append('null')
                            elif(dtypes[key]==2):
                                new_strs.append('{}'.format(int(float(c))))
                            elif(dtypes[key]==1):
                                new_strs.append(c.strip())
                            else:
                                new_strs.append('"{}"'.format(c.strip('"').strip()))

                        #row_strs = ['null' if c=='nan' else c for c in row_strs]
                        data_strs = data_strs + '[' + ', '.join(new_strs) + ']'
                        if(k<val.data.shape[0]-1):
                            data_strs = data_strs + ','
                    data_strs = data_strs + ']'
                    file.write('      "data": {}\n'.format(data_strs))
                else:
                    data_strs = [str(c).strip() for c in val.data]
                    new_strs = []
                    for c in data_strs:
                        if(c=='nan'):
                            new_strs.append('null')
                        elif(dtypes[key]==2):
                            new_strs.append('{}'.format(int(float(c))))
                        elif(dtypes[key]==1):
                            # the numbers come in as unquoted floats, don't change
                            new_strs.append(c.strip())
                        else:
                            new_strs.append('"{}"'.format(c.strip('"').strip()))
                    #data_strs = ['null' if c=='nan' else c for c in data_strs]
                    file.write('      "data": [{}]\n'.format(', '.join(new_strs)))
                #file.write('    },\n')
                if i == num_items - 1:
                    file.write('    }\n')
                else:
                    file.write('    },\n')
            file.write('  }\n')
            file.write('}\n')            
            
    # Test the file write
    if(debug):
        xmlroot = ET.parse(xmlfile).getroot()
        #code.interact(local=dict(globals(), **locals()))

    
    # Write out the xml file
    # Example template
    # <?xml version="1.0"?>
    # <all>
    #  <notes> This is the default FATES parameter file </notes>
    #  <history> timestamp: note </history>
    #  <dimensions>
    #   <dim name='numpft'> 14 </dim>
    #  </dimensions>
    #  <parameters>
    #   <par name='fates_stoich_nitr'
    #    <d>fates_pft,other</d>
    #    <l>long description</l>
    #    <u>units</u>
    #    <v>1,2,3,4</v>
    #   </par>
    #  </parameters>
    # </all>
    # 
    

    
def CDLParse(file_name,verbose):

    fp = open(file_name,"r")
    contents = fp.readlines()
    fp.close()

    il_dim = -1
    il_var = -1
    
    # contents is a list of strings, each index is a different line
    for il,line in enumerate(contents):
        if "dimensions:" in line:
            il_dim = il
            break

    dims = {}
    params = {}
        
    find_dims = True
    il = il_dim+1
    while(find_dims):
        line = contents[il]
        line = re.sub('\s+','',line)
        line = re.sub(';','',line)
        if "variables:" in line:
            il_var = il
            break
        else:
            lsplit = line.split('=')
            if(len(lsplit) != 2):
                print('failure reading dimensions')
                exit(2)

            dims[lsplit[0]] = int(lsplit[1])
            il=il+1

    if(verbose):
        print("Found the following dimensions:")
        for key, val in dims.items():
            print("  {} = {}".format(key,val))

    
    # lets look through the params

    if not ("variables" in contents[il_var]):
        print("did not find the variables section?")
        exit(2)

    vartypes = ["float","double","char","integer"]
    find_vars = True
    il = il_var+1
    while(find_vars):
        line = contents[il]
        line = re.sub('\s+',' ',line)
        line = re.sub(';','',line)
        line = line.strip()
        sline = line.split(' ')

        # determine if the any of the "type" characters match the first index
        # if there are, then we found a variable
        vartype = sline[0]
        if not (vartype in vartypes):
            find_vars = False
        else:
            if(len(sline)<2):
                print('strange variable descriptor')
                exit(0)
            varcombo = ''.join(sline[1:])
            vartype = sline[0]
            varkey = varcombo.split('(')[0]
            
            if(len(varcombo.split('('))>1):
                vardims = varcombo.split('(')[1].split(')')[0].split(',')
            else:
                # Scalar
                vardims = ['scalar']

            # Initialize the variable
            params[varkey] = param_type(varkey,vardims,vartype,dims)

            # Process the metadata
            find_var_meta = True
            while(find_var_meta):
                il=il+1
                line = contents[il]
                line = re.sub('\s+',' ',line)
                line = re.sub(r"&",'and',line)
                line = re.sub(';','',line)
                line = line.strip()
                sline = line.split(':')
                if(len(sline)>1 and sline[0] == varkey):
                    #code.interact(local=dict(globals(), **locals()))
                    part_sline = sline[1:]
                    #metacombo = sline[1]
                    metacombo = " ".join(part_sline)
                    metakey = metacombo.split('=')[0].strip()
                    metavar = ''.join(metacombo.split('=')[1:])
                    metavar = re.sub('\"',' ',metavar).strip()
                    params[varkey].add_meta(metakey,metavar)
                else:
                    # this is either the next variable new section
                    find_var_meta = False

            if(debug):
                params[varkey].print_var()


    # Add the data to the variable
    # -----------------------------------------------------------------------------------

    # Earmark the start of the data section
    find_data_start = True
    while(find_data_start):
        if not ("data:" in contents[il]):
            il = il+1
            if(il==10000):
                print('could not find start of data section')
                exit(2)
        else:
            il_dat = il
            find_data_start = False
        
    for symbol, var_obj in params.items():

        il = il_dat+1

        # search the file in the data section
        find_var = True
        while(find_var):
            if(il>len(contents)):
                print('Could not find the data for symbol {}'.format(symbol))
                exit(2)
            
            if symbol == contents[il].split('=')[0].strip():

                # its possible that the variable line of interest
                # is for a variable with a different variable's name embedded in it
                # for instance fates_dev_arbitrary can be found WITHIN fates_dev_arbitrary_pft
                # so lets try to find the exact symbol here

                if(debug):
                    print("Filling data for {}".format(symbol))
                    print("found: {}".format(contents[il]))
                
                find_data = True
                ilp = 0
                jline = ''
                while(find_data):
                    ilp=ilp+1
                    if(ilp>1000):
                        print('Having trouble finding data for symbol: {}'.format(symbol))
                        exit(2)
                        
                    # Take everything to the right of the equals sign
                    # keep going line by line unitl you hit a semi-colon
                
                    line = contents[il]
                    line = re.sub('\s+',' ',line)
                    line = line.strip()

                    # remove an equal-sign keep everything to the right
                    if('=' in line):
                        line = line.split('=')[1]
                    
                    if ';' in line:
                        find_data = False
                        find_var = False
                        line = re.sub(';','',line)
                    else:
                        il = il+1

                    # concatenate the string
                    jline=jline+line

                dlist = jline.split(',')

                # push the string of values to the parameter structure
                if(verbose):
                    print('Processing {}'.format(symbol))
                var_obj.add_data(dlist)
                
            il = il + 1
                    

    return params,dims

if __name__ == "__main__":
    main(sys.argv)
    
