import numpy as np
import os
import sys
import getopt
import code  # For development: code.interact(local=locals())
import time
import imp
PartehTypes = imp.load_source('PartehTypes', 'py_modules/PartehTypes.py')


# ========================================================================================
# Interpret the XML file

def load_xml(xmlfile):  #

    import xml.etree.ElementTree as et

    xmlroot = et.parse(xmlfile).getroot()
    print("\nOpenend: "+xmlfile)

    # Time control
    # -----------------------------------------------------------------------------------

    time_control      = PartehTypes.timetype()  # Initialize the time structure

    elem              = xmlroot.find('time_control')
    date_start_str    = elem.find('date_start').text
    date_stop_str     = elem.find('date_stop').text
    timestep_str      = elem.find('timestep_sec').text
    max_trunc_err_str = elem.find('max_trunc_error').text
    time_control.InitializeTime(date_start_str,date_stop_str,timestep_str,max_trunc_err_str)


    # FATES-PARTEH model parameters
    # -----------------------------------------------------------------------------------

    fates_cdl_file = xmlroot.find('fates_cdl_file').text

    # List of PFTs for the plants

    use_pfts_text = xmlroot.find('use_pfts').text
    use_pfts = []
    for use_pft in use_pfts_text.strip().split(','):
        use_pfts.append(int(use_pft))

    # Specify the boundary condition
    # -----------------------------------------------------------------------------------

    boundary_c_check = {}
    boundary_c_check['AllometricCarbon']=['DailyCFromCArea']
    boundary_c_check['AllometricCNP']=['DailyCNPFromCArea','DailyCNPFromStorageSinWaveNoMaint']

    boundary_method = xmlroot.find('boundary_formulation').text.strip()

#    if ( not any(x in boundary_method for x in boundary_c_check[parameters.hypothesis]) ):
#        print("A boundary condition formulation was not associated\n")
#        print(" with your hypothesis in the XML. Exiting.")
#        print("hypothesis: {}".format(parameters.hypothesis))
#        print("boundary formulation: {}".format(parameters.boundary_method))
#        exit(2)


    # Load up all the pft parameters that are specific to the Boundary Condition method
    # Must add a check to see if all correct parameters are loaded
    # -----------------------------------------------------------------------------------

    driver_params = {}

    driver_params_root = xmlroot.find('driver_parameters')

    for par_idx, par_elem in enumerate(driver_params_root.iter('pft_par')):

        pft_param_name = par_elem.attrib['name'].strip()
        pft_param_val  = par_elem.text.strip()
        pft_vector = [float(i) for i in pft_param_val.split(',')]
        if (len(pft_vector) != len(use_pfts)):
            print('PFT parameters in xml file must have as as many entries as use_pfts')
            print('{} is: {}'.format(pft_param_name,pft_param_val))
            print('exiting')
            exit(1)

        driver_params[pft_param_name] = PartehTypes.driver_param_type()
        #        driver_params.append(driver_param)
        for idx,value in enumerate(pft_vector):
            driver_params[pft_param_name].param_vals.append(value)




    print("\n\n Completed Interpreting: "+xmlfile)


    return(time_control, fates_cdl_file, driver_params, boundary_method, use_pfts)
