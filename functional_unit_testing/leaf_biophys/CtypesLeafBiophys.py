import re
import subprocess
import sys
import ctypes
from ctypes import *
from operator import add
sys.path.append('../shared/py_src')
from PyF90Utils import c8, ci, cchar, c8_arr, ci_arr, ccharnb

# Subroutines
# =======================================================================================

def GetModSymbol(mod_path,symbol):

    # This routine uses some minor string magic and a system call using
    # "nm" to list the symbols inside a fortran object.  Symbols are simply
    # the text strings the define the routines and data structures inside
    # the compiled objects. When symbols are created inside a compiled object,
    # different compilers do different things, like force to lower or upper
    # case, or attach trailing or preceding underscores.  This routine
    # helps identify the exact symbol name, from the name of the fortran
    # routine (as it appears in code) that we want to match with.
    
    subprocess.run(["nm", mod_path],capture_output=True)
    
    list = str(subprocess.run(["nm",mod_path],capture_output=True)).split()
    p = re.compile(symbol, re.IGNORECASE)
    modlist = [ s for s in list if p.search(s)]

    # If nothing came up, thats bad
    if (len(modlist)==0):
        print('Failed to find the right module symbol for:{} in module {}'.format(symbol,mod_path))
        print(list)
        print('Exiting')
        exit(2)

    # Its possible a routine name could also be part of another routine name,
    # such as alloc_such_and_such  vs dealloc_such_and_such, we want the shorter
    mod_symbols = []
    slen = 10000
    for item in modlist:
        ms_str = item.split('\\')[0]
        if len(ms_str)<slen:
            mod_symbol = ms_str
            slen = len(ms_str)

        
    return mod_symbol


# =======================================================================================

# Instantiate the F90 modules

f90_const_obj = ctypes.CDLL('bld/FatesConstantsMod.o',mode=ctypes.RTLD_GLOBAL)
f90_shr_obj = ctypes.CDLL('bld/WrapShrMod.o',mode=ctypes.RTLD_GLOBAL)
f90_fatesutils_obj = ctypes.CDLL('bld/FatesUtilsMod.o',mode=ctypes.RTLD_GLOBAL)
f90_leaf_biophys_obj = ctypes.CDLL('bld/LeafBiophysicsMod.o',mode=ctypes.RTLD_GLOBAL)
f90_leaf_biophys_supp_obj = ctypes.CDLL('bld/LeafBiophysSuppMod.o',mode=ctypes.RTLD_GLOBAL)

# Identify subroutine objects, so we can call them
f90_set_leaf_param_sub = getattr(f90_leaf_biophys_supp_obj, GetModSymbol('bld/LeafBiophysSuppMod.o','setleafparam'))
f90_alloc_leaf_param_sub = getattr(f90_leaf_biophys_supp_obj, GetModSymbol('bld/LeafBiophysSuppMod.o','allocleafparam'))
f90_dealloc_leaf_param_sub = getattr(f90_leaf_biophys_supp_obj, GetModSymbol('bld/LeafBiophysSuppMod.o','deallocleafparam'))
f90_dump_param_sub =  getattr(f90_leaf_biophys_supp_obj, GetModSymbol('bld/LeafBiophysSuppMod.o','DumpParams'))
f90_set_leaf_param_sub.argtypes = [POINTER(c_double),POINTER(c_int),c_char_p,c_long]
f90_biophysrate_sub = getattr(f90_leaf_biophys_obj,GetModSymbol('bld/LeafBiophysicsMod.o','LeafLayerBiophysicalRate'))
f90_leaflayerphoto_sub = getattr(f90_leaf_biophys_obj,GetModSymbol('bld/LeafBiophysicsMod.o','LeafLayerPhotosynthesis'))
f90_qsat_sub = getattr(f90_leaf_biophys_obj,GetModSymbol('bld/LeafBiophysicsMod.o','QSat'))
f90_cangas_sub = getattr(f90_leaf_biophys_obj,GetModSymbol('bld/LeafBiophysicsMod.o','GetCanopyGasParameters'))
f90_lmr_ryan_sub = getattr(f90_leaf_biophys_obj,GetModSymbol('bld/LeafBiophysicsMod.o','LeafLayerMaintenanceRespiration_Ryan_1991'))
f90_lmr_atkin_sub = getattr(f90_leaf_biophys_obj,GetModSymbol('bld/LeafBiophysicsMod.o','LeafLayerMaintenanceRespiration_Atkin_etal_2017'))
f90_agross_rubiscoc3  = getattr(f90_leaf_biophys_obj,GetModSymbol('bld/LeafBiophysicsMod.o','AgrossRubiscoC3'))
f90_agross_rubpc3  = getattr(f90_leaf_biophys_obj,GetModSymbol('bld/LeafBiophysicsMod.o','AgrossRuBPC3'))
f90_agross_rubpc4  = getattr(f90_leaf_biophys_obj,GetModSymbol('bld/LeafBiophysicsMod.o','AgrossRuBPC4'))
f90_agross_pepc4  = getattr(f90_leaf_biophys_obj,GetModSymbol('bld/LeafBiophysicsMod.o','AgrossPEPC4'))

f90_gs_medlyn = getattr(f90_leaf_biophys_obj,GetModSymbol('bld/LeafBiophysicsMod.o','StomatalCondMedlyn'))
f90_gs_ballberry = getattr(f90_leaf_biophys_obj,GetModSymbol('bld/LeafBiophysicsMod.o','StomatalCondBallBerry'))
f90_velotomolarcf_sub = getattr(f90_leaf_biophys_obj,GetModSymbol('bld/LeafBiophysicsMod.o','VeloToMolarCF'))
f90_cifunc_sub = getattr(f90_leaf_biophys_obj,GetModSymbol('bld/LeafBiophysicsMod.o','CiFunc'))
f90_cibisection_sub = getattr(f90_leaf_biophys_obj,GetModSymbol('bld/LeafBiophysicsMod.o','CiBisection'))

# For functions, define the return value
f90_agross_rubiscoc3.restype = c_double
f90_agross_rubpc3.restype = c_double
f90_agross_rubpc4.restype = c_double
f90_agross_pepc4.restype = c_double
f90_velotomolarcf_sub.restype = c_double
