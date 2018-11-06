import numpy as np
import os
import sys
import getopt
import code  # For development: code.interact(local=locals())
import time


# =======================================================================================
#
# Global Parameters
#
# =======================================================================================


os.environ['TZ'] = 'UTC'
time.tzset()

time_precision = 1.0e-10  # Acceptable time error for the 
                          # adaptive time-stepper

class param_type:

    def __init__(self):

        # Initialize the list of parameters

        self.hypothesis = ""

        self.boundary_method = ""

        # These are passed to the PARTEH Fortran code
        # This is a list
        self.parteh_pfts = []

        # This is a list of the organ names
        # These names must be consistent
        # with the indices provided in the parameter file
        # and that those indices should match the global
        # indices in PRTGenericMod.F90
        self.parteh_organs = []

        # These are used in the boundary conditions
        self.boundary_pfts = []

        # Save the number of pfts (as a convencience)
        self.numpfts = -9
        
        # Add other parameter groups as we go

class pft_type:

    def __init__(self,pft_name):

        # Initialize a dictionary of parameters for any pft
        self.name      = pft_name
        self.param_dic = {}


class diagnostics_type:

    def __init__(self):

        self.dates    = []
        self.dbh      = []
        self.dailyc   = []
        self.leaf_c   = []
        self.fnrt_c   = []
        self.sapw_c   = []
        self.store_c  = []
        self.struct_c = [] 
        self.repro_c  = []
        self.leaf_cturn   = []
        self.fnrt_cturn   = []
        self.sapw_cturn   = []
        self.store_cturn  = []
        self.struct_cturn = []
        
        self.leaf_n   = []
        self.fnrt_n   = []
        self.sapw_n   = []
        self.store_n  = []
        self.struct_n = [] 
        self.repro_n  = []
        self.leaf_nturn   = []
        self.fnrt_nturn   = []
        self.sapw_nturn   = []
        self.store_nturn  = []
        self.struct_nturn = []

        self.leaf_p   = []
        self.fnrt_p   = []
        self.sapw_p   = []
        self.store_p  = []
        self.struct_p = [] 
        self.repro_p  = []
        self.leaf_pturn   = []
        self.fnrt_pturn   = []
        self.sapw_pturn   = []
        self.store_pturn  = []
        self.struct_pturn = []

        self.crown_area   = []
        self.root_c_exudate = []
        self.root_n_exudate = []
        self.root_p_exudate = []
        self.growth_resp = []


## Define the state variables and state terms types

class timetype:

    def __init__(self):

        self.datetime_start = np.datetime64("1600-01-01")
        self.datetime_stop  = np.datetime64("1400-01-01")
        self.datetime       = np.datetime64("1300-01-01")
        self.dt_fullstep    = np.timedelta64(int(86400),'s')
        self.sim_complete   = False
        self.max_err        = -9.9
        self.id_substep     = -9
        self.dt_substep     = np.timedelta64(int(3600),'s')
        self.dt_optsubstep  = np.timedelta64(int(3600),'s')


    def InitializeTime(self,date_start_str,date_stop_str,timestep_str,max_trunc_err_str):

        # Perform checks here as well
        date_start_str    = date_start_str.strip()
        date_stop_str     = date_stop_str.strip()
        timestep_str      = timestep_str.strip()
        max_trunc_err_str = max_trunc_err_str.strip()

        # Timing for the main time loop
        # -------------------------------------------------------------------------------
        self.datetime_start = np.datetime64(date_start_str)
        self.datetime_stop  = np.datetime64(date_stop_str)
        self.datetime       = self.datetime_start
        self.dt_fullstep    = float(timestep_str)
        self.sim_complete   = False
   
        # Maximum allowable truncation error on iterator
        self.max_err        = float(max_trunc_err_str)
        
        
        # Timing for the integrator
        # -------------------------------------------------------------------------------
        self.id_substep  = 0
        self.dt_substep    = self.dt_fullstep
        self.dt_optsubstep = self.dt_fullstep

    def ResetTime(self):

        self.datetime      = self.datetime_start
        self.id_substep    = 0
        self.dt_substep    = self.dt_fullstep
        self.dt_optsubstep = self.dt_fullstep
        self.sim_complete  = False

    def UpdateTime(self):

        self.datetime += np.timedelta64(int(self.dt_fullstep),'s')
        if(self.datetime >= self.datetime_stop):
            self.sim_complete = True
            
    def CheckFullStepTime(self,targettime):
        if(np.abs(self.datetime-targettime)>time_precision):
            print('The adaptive time-stepper finished')
            print(' on a time-stamp that does not match')
            print(' the projected timestep')
            print(' projected: {}'.format(targettime))
            print(' actual:    {}'.format(self.datetime))
            print(' exiting')
            exit(2)
        else:
            self.datetime = targettime
            


    def UpdatePartialTime(self,dt_seconds):
        self.datetime += np.timedelta64(int(dt_seconds),'s')

        
