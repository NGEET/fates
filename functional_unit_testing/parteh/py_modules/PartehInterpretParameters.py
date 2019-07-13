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

def load_xml(xmlfile, time_control, parameters ):

    import xml.etree.ElementTree as et


    xmlroot = et.parse(xmlfile).getroot()
    print("\nOpenend: "+xmlfile)



    # Time control
    # -----------------------------------------------------------------------------------

    elem              = xmlroot.find('time_control')
    date_start_str    = elem.find('date_start').text
    date_stop_str     = elem.find('date_stop').text
    timestep_str      = elem.find('timestep_sec').text
    max_trunc_err_str = elem.find('max_trunc_error').text
    time_control.InitializeTime(date_start_str,date_stop_str,timestep_str,max_trunc_err_str)

    # PARTEH model parameters

    # Read in the hypothesis we are testing

    hypotheses = ('AllometricCarbon','AllometricCNP')

    hypothesis_root       = xmlroot.find('hypothesis')
    parameters.hypothesis = hypothesis_root.text.strip()

    try:
        parameters.prt_model  = hypotheses.index(parameters.hypothesis) + 1
    except ValueError:
        print('Attempted to identify PARTEH model type: {}'.format(parameters.hypothesis))
        print('Not in the list: {}'.format(hypotheses))
        exit(1)

    boundary_c_check = {}
    boundary_c_check['AllometricCarbon']=['DailyCFromCArea']
    boundary_c_check['AllometricCNP']=['DailyCNPFromCArea','DailyCNPFromStorageSinWaveNoMaint']

    boundary_root = xmlroot.find('boundary_formulation')
    parameters.boundary_method = boundary_root.text.strip()

    if ( not any(x in parameters.boundary_method for x in boundary_c_check[parameters.hypothesis]) ):
        print("A boundary condition formulation was not associated\n")
        print(" with your hypothesis in the XML. Exiting.")
        print("hypothesis: {}".format(parameters.hypothesis))
        print("boundary formulation: {}".format(parameters.boundary_method))
        exit(2)

    # Load up all the pft parameters that are specific to the Boundary Condition method
    # Must add a check to see if all correct parameters are loaded
    # -----------------------------------------------------------------------------------

    for ptype_idx, ptype_elem in enumerate(parameters_root.iter('boundary_parameters')):

        for par_idx, par_elem in enumerate(ptype_elem.iter('pft_par')):

            pft_param_name = par_elem.attrib['name'].strip()
            pft_param_val  = par_elem.text.strip()
            pft_vector = [float(i) for i in pft_param_val.split(',')]
            if (len(pft_vector) != parameters.num_pfts):
                print('parameter was given no value?')
                print('{} is: {}'.format(pft_param_name,pft_param_val))
                print('exiting')
                exit(1)

            for idx,value in enumerate(pft_vector):
                parameters.boundary_pfts[idx].param_dic[pft_param_name] = value



    print("\n\n Completed Interpreting: "+xmlfile)
    print("\n Found {} PFT(s)".format(parameters.num_pfts))
