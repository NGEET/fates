# =======================================================================================
#
# This python module contains routines and utilities which will interpret the
# FATES fortran code base to return information on the use of parameters.
# This does not parse the CDL or NC files, this only parses the fortran code.
#
# This module will help:
#   1) List the parameters found in a given file
#   2) Determine the parameter names found therein
#   3) Determine the parameter's name in the parameter file
#
# Note: This module can be used to determine usage of any sybmol associated with
#       the instantiation of a structure.  Ie, you can search for all parameters
#       in the 'EDPftvarcon_inst%' structure.  In FATES, the EDParamsMod and SFParamsMod
#       don't use a structure to hold their parameters though.
#
# =======================================================================================

import code  # For development: code.interact(local=dict(globals(), **locals()))


class f90_param_type:

    # -----------------------------------------------
    # PFTParamType stucture. A list of these will be
    # generated that denotes the PFT parameters used.
    # -----------------------------------------------

    def __init__(self,var_sym):

        self.var_sym  = var_sym   # Name of parameter in FORTRAN code
        self.var_name = ''        # Parameter's name in the parameter file


def GetSymbolUsage(filename,checkstr_in):

    # ---------------------------------------------------------------------
    # This procedure will check a fortran file and return a list (non-unique)
    # of all the PFT parameters found in the code.
    # Note: This will only determine the symbol name in code, this will
    #       not determine the symbol name in the parameter file.
    # ---------------------------------------------------------------------

    checkstr = checkstr_in.lower()

    f = open(filename,"r")
    contents = f.readlines()
    f.close()

    strclose = ',)( '

    var_list = []
    found = False


    for line in contents:
        if checkstr in line.lower():

            if(checkstr[-1] != '%'):
                print('The GetSymbolUsage() procedure requires')
                print(' that a structure ending with % is passed in')
                print(' check_str: --{}--'.format(check_str))
                exit(2)

            # We compare all in lower-case
            # There may be more than one parameter in a line,
            # so evaluate, pop-off, and try again

            substr   = line.lower()

            search_substr=True

            while(search_substr):

                p1 = substr.find(checkstr)+len(checkstr)

                pcomment = substr.find('!')
                if(pcomment<0):
                    pcomment=1000

                # This makes sure that if the line
                # has a comment, that it does not come before
                # the parameter symbol

                if( (p1>len(checkstr)) and (p1 < pcomment)):
                    found = True

                    # Identify the symbol by starting at the first
                    # character after the %, and ending at a list
                    # of possible symols including space
                    substr2=substr[p1:]
                    pend0=-1
                    for ch in strclose:
                        pend = substr2.find(ch)
                        if(pend>0):
                            substr2=substr2[:pend]
                            pend0=pend

                    var_list.append(f90_param_type(substr2))
                    if(pend0!=-1):
                        substr=substr[pend0:]
                    else:
                        print('Could not correctly identify the parameter string')
                        exit(2)

                else:
                    search_substr=False



    if(not found):
        print('No parameters with prefix: {}'.format(checkstr))
        print('were found in file: {}'.format(filename))
        print('If this is expected, remove that file from search list.')
        exit(2)


    return(var_list)




def GetPFTParmFileSymbols(var_list,pft_filename):

    #---------------------------------------------------------------
    # This procedure will determine the parameter file symbol/name
    # for a given PFT parameter name.  This relies on specific
    # file syntax in the PFT definitions file, so this is specific
    # only to PFT parameters.
    # --------------------------------------------------------------

    f = open(pft_filename,"r")
    contents = f.readlines()
    f.close()

    var_name_list = []
    for var in var_list:
        for i,line in enumerate(contents):
            if (var.var_sym in line) and ('data' in line) and ('=' in line):
                var.var_name = contents[i-2].split()[-1].strip('\'')

    return(var_list)


def MakeListUnique(list_in):

    # This procedure simply filters
    # an input list and returns the unique entries

    unique_list = []
    for var in list_in:
        found = False
        for uvar in unique_list:
            if (var.var_sym == uvar.var_sym):
                found = True
        if(not found):
            unique_list.append(var)

    return(unique_list)
