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
