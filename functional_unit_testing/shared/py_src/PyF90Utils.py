import ctypes
from ctypes import *


# These are shortcuts that simply convert real, ints and stings
# into refernced cytype compliant arguments

def c8(r8):
    return(byref(c_double(r8)))


def ci(i8):
    return(byref(c_int(i8)))


def cchar(fchar):
    return(byref(c_char(fchar)))


def c8_arr(r8_list):
    return(byref((len(r8_list) * c_double)(r8_list)))


def ci_arr(int_list):
    return(byref((len(int_list) * c_int)(int_list)))
