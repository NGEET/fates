#!/usr/bin/env python

#### this script modifies the default FATES parameter file to generate
#    a file used in testing E3SM
#    Parser code was based off of modify_fates_paramfile.py

import os
import argparse
import code  # For development: code.interact(local=dict(globals(), **locals()))
import xml.etree.ElementTree as et
# Newer versions of scipy have dropped the netcdf module and
# netcdf functions are part of the io parent module
try:
    from scipy import io as nc

except ImportError:
    from scipy.io import netcdf as nc
	 
debug = True

# ---------------------------------------------------------------------------------------

class param_type:
    def __init__(self,name,values_text):
        self.name = name
        self.values = values_text.replace(" ","") #[float(x) for x in values_text.split(',')]
        
# ---------------------------------------------------------------------------------------


def load_xml(xmlfile): 

    

    xmlroot = et.parse(xmlfile).getroot()
    print("\nOpenend: "+xmlfile)

    base_cdl = xmlroot.find('base_file').text
    new_cdl = xmlroot.find('new_file').text

    pftparams = xmlroot.find('pft_trim_list').text.replace(" ","")
    
    paramroot = xmlroot.find('parameters') 

    grouplist = []
    for group in paramroot:

        if(not('ids' in group.attrib.keys())):
            print("pft_mod_group must have an ids attribute with comma delimited pft indices")
            print("exiting")
            exit(2)

        pft_str = group.attrib['ids'].strip()
        pft_strvec = pft_str.split(',')
        pft_ids = [int(s) for s in pft_strvec]
        
        print("Processing PFT group: "+group.tag+" ids: "+pft_str)
            
        #paramlist = []
        #for param in paramroot:
        #    print("parsing "+param.tag)
        #    paramlist.append(param_type(param.tag,param.text))
        
        
    code.interact(local=dict(globals(), **locals()))
    return(base_cdl,new_cdl,pftparams,paramlist)

    

# Little function for assembling the call to the system to make the modification
# ----------------------------------------------------------------------------------------

def parse_syscall_str(fnamein,fnameout,pft_index,param_name,param_val):

    if(pft_index==0):
        sys_call_str = "../tools/modify_fates_paramfile.py"+" --fin " + fnamein + \
            " --fout " + fnameout + " --var " + param_name + " --silent " +\
            " --val " + "\" "+param_val+"\"" + " --nohist --overwrite --all"
    else:
        pft_str_index="{}".format(pft_index)
        sys_call_str = "../tools/modify_fates_paramfile.py"+" --fin " + fnamein + \
            " --fout " + fnameout + " --var " + param_name + " --silent " +\
            " --val " + "\" "+param_val+"\"" + " --nohist --overwrite --pft "+pft_str_index

    if(debug):
        print(sys_call_str)
    
    return(sys_call_str)



def main():

    # Parse arguments
    parser = argparse.ArgumentParser(description='Parse command line arguments to this script.')
    parser.add_argument('--f', dest='xmlfile', type=str, help="XML control file  Required.", required=True)
    args = parser.parse_args()


    # Load the xml file, which contains the base cdl, the output cdl,
    #  and the parameters to be modified
    #[base_cdl,new_cdl,pftlist,paramlist] = load_xml(args.xmlfile)

    xmlroot = et.parse(args.xmlfile).getroot()
    print("\nOpened: "+args.xmlfile)

    base_cdl = xmlroot.find('base_file').text
    new_cdl = xmlroot.find('new_file').text

    # Append extension nc to temporary files 
    # (in some netcdf versions, the lack of extension causes failures)
    ext_nc  = ".nc"

    # Convert the base cdl file into a temp nc binary
    base_nc = os.popen('mktemp').read().rstrip('\n')+ext_nc
    gencmd = "ncgen -o "+base_nc+" "+base_cdl
    os.system(gencmd)
	 
    # Generate a temp output file name
    new_nc = os.popen('mktemp').read().rstrip('\n')+ext_nc

    os.system("ls "+base_nc)
    os.system("ls "+new_nc)
	 
    # Use FatesPFTIndexSwapper.py to prune out unwanted PFTs
    pft_trim_list = xmlroot.find('pft_trim_list').text.replace(" ","")
    swapcmd="../tools/FatesPFTIndexSwapper.py --pft-indices="+pft_trim_list+" --fin="+base_nc+" --fout="+new_nc+" --nohist"   #+" 1>/dev/null"
    os.system(swapcmd)
	 
    # On subsequent parameters, overwrite the file

    paramroot = xmlroot.find('parameters') 

    grouplist = []
    for group in paramroot:

        if(group.tag.strip() == 'non_pft_group'):

            print("Processing non_pft_group")
            
            for param in group:
                print("parsing "+param.tag)
            
                change_str = parse_syscall_str(new_nc,new_nc,0,param.tag,param.text.replace(" ",""))
                os.system(change_str)

        elif(group.tag.strip() == 'pft_group'):
        
            if(not('ids' in group.attrib.keys())):
                print("pft_mod_group must have an ids attribute with comma delimited pft indices")
                print("exiting")
                exit(2)

            pft_str = group.attrib['ids'].strip()
            pft_strvec = pft_str.split(',')
            pft_ids = [int(s) for s in pft_strvec]
        
            print("Processing PFT group: "+group.tag+" ids: "+pft_str)
            
            for param in group:
                print("parsing "+param.tag)
                
                param_vec = [str_val for str_val in param.text.replace(" ", "").split(',')]

                # The number of parameters does not need to equal to the number of PFTs
                # but it does need to be equally divisible by it
                v_per_p = float(len(param_vec))/float(len(pft_ids))

                if( int(100.0*v_per_p) != 100*int(v_per_p) ):
                    print("inconsistent number of parameter values and number of pfts specified")
                    print("exiting")
                    exit(2)
                                
                for j,pft_id in enumerate(pft_ids):
                    j0 = j*int(v_per_p)
                    j1 = (j+1)*int(v_per_p)
                    #code.interact(local=dict(globals(), **locals()))
                    param_str = ",".join(param_vec[j0:j1])
                    change_str = parse_syscall_str(new_nc,new_nc,pft_id,param.tag,param_str)
                    os.system(change_str)

        else:
            print("Unidentified group, should be: non_pft_group or pft_group")
            print("exiting")
            exit(2)


    # Append history
    fp_nc  = nc.netcdf_file(new_nc, 'a')
    fp_nc.history = "This file was generated by BatchPatchParams.py:\n"\
                    "CDL Base File = {}\n"\
                    "XML patch file = {}"\
                     .format(base_cdl,args.xmlfile)
    fp_nc.close()

    # Sort the new file
    newer_nc = os.popen('mktemp').read().rstrip('\n')+ext_nc
    os.system("../tools/ncvarsort.py --fin "+new_nc+" --fout "+newer_nc+" --overwrite")


    # Dump the new file to the cdl
    os.system("ncdump "+newer_nc+" > "+new_cdl)

    print("\nBatch parameter transfer complete\n")
    print("\nGenerated: {}\n".format(new_cdl))
        
# This is the actual call to main

if __name__ == "__main__":
    main()
