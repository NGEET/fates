#!/usr/bin/env python

#### this script modifies the default FATES parameter file to generate
#    a file used in testing E3SM
#    Parser code was based off of modify_fates_paramfile.py

import os
import argparse
import code  # For development: code.interact(local=dict(globals(), **locals()))
from scipy.io import netcdf


# ---------------------------------------------------------------------------------------

class param_type:
    def __init__(self,name,values_text):
        self.name = name
        self.values = values_text.replace(" ","") #[float(x) for x in values_text.split(',')]
        
# ---------------------------------------------------------------------------------------


def load_xml(xmlfile): 

    import xml.etree.ElementTree as et

    xmlroot = et.parse(xmlfile).getroot()
    print("\nOpenend: "+xmlfile)

    base_cdl = xmlroot.find('base_file').text
    new_cdl = xmlroot.find('new_file').text

    pftparams = xmlroot.find('pft_list').text.replace(" ","")
    
    paramroot = xmlroot.find('parameters') 
    paramlist = []
    for param in paramroot:
        print("parsing "+param.tag)
        paramlist.append(param_type(param.tag,param.text))
        
        

    return(base_cdl,new_cdl,pftparams,paramlist)

    

# Little function for assembling the call to the system to make the modification
# ----------------------------------------------------------------------------------------

def parse_syscall_str(fnamein,fnameout,param_name,param_val):

    sys_call_str = "../tools/modify_fates_paramfile.py"+" --fin " + fnamein + \
                   " --fout " + fnameout + " --var " + param_name + " --silent " +\
                   " --val " + "\" "+param_val+"\"" + " --overwrite --all"


    print(sys_call_str)
    
    return(sys_call_str)



def main():

    # Parse arguments
    parser = argparse.ArgumentParser(description='Parse command line arguments to this script.')
    parser.add_argument('--f', dest='xmlfile', type=str, help="XML control file  Required.", required=True)
    args = parser.parse_args()


    # Load the xml file, which contains the base cdl, the output cdl,
    #  and the parameters to be modified
    [base_cdl,new_cdl,pftlist,paramlist] = load_xml(args.xmlfile)


    # Convert the base cdl file into a temp nc binary
    base_nc = os.popen('mktemp').read().rstrip('\n')
    gencmd = "ncgen -o "+base_nc+" "+base_cdl
    os.system(gencmd)
    
    # Generate a temp output file name
    new_nc = os.popen('mktemp').read().rstrip('\n')


    # Use FatesPFTIndexSwapper.py to prune out unwanted PFTs
    swapcmd="../tools/FatesPFTIndexSwapper.py --pft-indices="+pftlist+" --fin="+base_nc+" --fout="+new_nc   #+" 1>/dev/null"
    os.system(swapcmd)

    # We open the new parameter file. We only use this
    # to do some dimension checking. 
    fp_nc  = netcdf.netcdf_file(base_nc, 'r')
    
    # On subsequent parameters, overwrite the file
    for param in paramlist:

        change_str = parse_syscall_str(new_nc,new_nc,param.name,param.values)
        os.system(change_str)

    # Sort the new file
    newer_nc = os.popen('mktemp').read().rstrip('\n')
    os.system("../tools/ncvarsort.py --fin "+new_nc+" --fout "+newer_nc+" --overwrite")

    # Dump the new file to the cdl
    os.system("ncdump "+newer_nc+" > "+new_cdl)

    fp_nc.close()

    print("\nBatch parameter transfer complete\n")
    
        
# This is the actual call to main

if __name__ == "__main__":
    main()
