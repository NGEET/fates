import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib as mp
import ctypes
from ctypes import * #byref, cdll, c_int, c_double, c_char_p, c_long
import xml.etree.ElementTree as ET
import argparse
import re              # This is a heftier string parser
import code  # For development: code.interact(local=dict(globals(), **locals()))


# =======================================================================================
# Set some constants. If they are used as constant arguments to the F90 routines,
#  define them with their ctype identifiers
# =======================================================================================

ndbh = 200
maxdbh = 50
ccanopy_trim = c_double(1.0)                # Crown Trim (0=0% of target, 1=100% of targ)
csite_spread = c_double(0.0)                # Canopy spread (0=closed, 1=open)
cnplant      = c_double(1.0)                # Number of plants (don't change)
cilayer      = c_int(1)                     # Index of the plant's canopy layer
ccanopy_lai  = (2 * c_double)(1.0,1.0)      # The LAI of the different canopy layers
                                            # THIS VECTOR MUST MATCH ncanlayer
cdo_reverse  = c_bool(0)                    # DO NOT GET REVERSE CROWN AREA

# =======================================================================================
# Setup references to fortran shared libraries
# =======================================================================================

allom_const_object = "./include/FatesConstantsMod.o"
allom_wrap_object = "./include/AllomUnitWrap.o"
allom_lib_object = "./include/FatesAllometryMod.o"

# ==============================================================================
# Instantiate fortran allometry and other libraries
# ==============================================================================

f90constlib= ctypes.CDLL(allom_const_object,mode=ctypes.RTLD_GLOBAL)
f90wraplib = ctypes.CDLL(allom_wrap_object,mode=ctypes.RTLD_GLOBAL)
f90funclib = ctypes.CDLL(allom_lib_object,mode=ctypes.RTLD_GLOBAL)

# =======================================================================================
# Create aliases to all of the different routines, set return types for functions
# =======================================================================================

f90_pftalloc  = f90wraplib.__edpftvarcon_MOD_edpftvarconalloc    #(numpft)
f90_pftset    = f90wraplib.__edpftvarcon_MOD_edpftvarconpyset
f90_pftset.argtypes = [POINTER(c_int),POINTER(c_double),POINTER(c_int),c_char_p,c_long]
f90_h2d       = f90funclib.__fatesallometrymod_MOD_h2d_allom     #(h,ipft,d,dddh)
f90_h         = f90funclib.__fatesallometrymod_MOD_h_allom       #(d,ipft,h,dhdd)
f90_bagw      = f90funclib.__fatesallometrymod_MOD_bagw_allom    #(d,ipft,bagw,dbagwdd)
f90_bleaf     = f90funclib.__fatesallometrymod_MOD_bleaf         #(d,ipft,canopy_trim,bl,dbldd)
f90_bsap      = f90funclib.__fatesallometrymod_MOD_bsap_allom    #(d,ipft,canopy_trim,asapw,bsap,dbsapdd)
f90_bstore    = f90funclib.__fatesallometrymod_MOD_bstore_allom  #(d,ipft,canopy_trim,bstore,dbstoredd)
f90_bbgw      = f90funclib.__fatesallometrymod_MOD_bbgw_allom    #(d,ipft,canopy_trim,bbgw,dbbgwdd)
f90_bfineroot = f90funclib.__fatesallometrymod_MOD_bfineroot     #(d,ipft,canopy_trim,bfr,dbfrdd)
f90_bdead     = f90funclib.__fatesallometrymod_MOD_bdead_allom   #(bagw,bbgw,bsap,ipft,bdead,dbagwdd,dbbgwdd,dbsapdd,dbdeaddd)
f90_carea     = f90funclib.__fatesallometrymod_MOD_carea_allom   #(d,nplant,site_spread,ipft,c_area)(d,nplant,site_spread,ipft,c_area)
f90_treelai   = f90funclib.__fatesallometrymod_MOD_tree_lai      #(leaf_c, pft, c_area, nplant, cl, canopy_lai, vcmax25top)
f90_treelai.restype = c_double


# This is the object type that holds our parameters
# =======================================================================================
class parameter:

    def __init__(self,symbol):

        self.dtype  = -9
        self.symbol = symbol
        self.vals   = []

    def setval(self,val,ipft):

        self.vals[ipft] = val

# This is just a helper script that generates random colors
# =======================================================================================
def DiscreteCubeHelix(N):

    base = plt.cm.get_cmap('cubehelix')
    np.random.seed(1)
    color_list = base(np.random.randint(0,high=255,size=N))
    cmap_name = base.name + str(N)
    return base.from_list(cmap_name, color_list, N)


# This will look through a CDL file for the provided parameter and determine
# the parameter's type, as well as fill an array with the data
# =======================================================================================
def CDLParse(file_name,parm):

    fp = open(file_name,"r")
    contents = fp.readlines()
    fp.close()

    # Look in the file for the parameters
    # symbol/name, record the line number
    iline=-1
    isfirst = True
    for i,line in enumerate(contents):
        if(parm.symbol in line):
            iline=i
            if(isfirst):
                dtype = line.split()[0]
                if(dtype.strip()=="float" or (dtype.strip()=="double")):
                    parm.dtype = 0
                elif(dtype.strip()=="char"):
                    parm.dtype = 1
                isFirst=False

    if(iline==-1):
        print('Could not find symbol: {} in file: {}'.format(parm.symbol,file_name))
        exit(2)
    else:
        search_field=True
        line=""
        lcount=0
        while(search_field and (lcount<100)):
            line+=contents[iline]
            if(line.count(';')>0):
                search_field=False
            else:
                search_field=True
            lcount=lcount+1
            iline=iline+1

        # Parse the line
        line_split = re.split(',|=',line)
        # Remove the variable name entry
        del line_split[0]

        # This is for read numbers
        if(parm.dtype == 0):
            ival=0
            for str0 in line_split:
                str=""
                isnum=False
                for s in str0:
                    if (s.isdigit() or s=='.'):
                        str+=s
                        isnum=True
                if(isnum):
                    parm.vals.append(float(str))

        # This is a sting
        elif(parm.dtype == 1):
            for str0 in line_split:
                # Loop several times to trim stuff off
                for i in range(5):
                    str0=str0.strip().strip('\"').strip(';').strip()
                parm.vals.append(str0)

    return(parm)



# Read in the arguments
# =======================================================================================

parser = argparse.ArgumentParser(description='Parse command line arguments to this script.')
parser.add_argument('--fin', '--input', dest='fnamein', type=str, help="Input CDL filename.  Required.", required=True)
args = parser.parse_args()


# Read in the parameters of interest that are used in the fortran objects. These
# parameters will be passed to the fortran allocation.
# =======================================================================================

parms = {}
parms['dbh_maxheight'] = CDLParse(args.fnamein,parameter('fates_allom_dbh_maxheight'))
parms['hmode'] = CDLParse(args.fnamein,parameter('fates_allom_hmode'))
parms['amode'] = CDLParse(args.fnamein,parameter('fates_allom_amode'))
parms['lmode'] = CDLParse(args.fnamein,parameter('fates_allom_lmode'))
parms['smode'] = CDLParse(args.fnamein,parameter('fates_allom_smode'))
parms['cmode'] = CDLParse(args.fnamein,parameter('fates_allom_cmode'))
parms['fmode'] = CDLParse(args.fnamein,parameter('fates_allom_fmode'))
parms['stmode'] = CDLParse(args.fnamein,parameter('fates_allom_stmode'))
parms['cushion'] = CDLParse(args.fnamein,parameter('fates_alloc_storage_cushion'))
parms['d2h1'] = CDLParse(args.fnamein,parameter('fates_allom_d2h1'))
parms['d2h2'] = CDLParse(args.fnamein,parameter('fates_allom_d2h2'))
parms['d2h3'] = CDLParse(args.fnamein,parameter('fates_allom_d2h3'))
parms['agb1'] = CDLParse(args.fnamein,parameter('fates_allom_agb1'))
parms['agb2'] = CDLParse(args.fnamein,parameter('fates_allom_agb2'))
parms['agb3'] = CDLParse(args.fnamein,parameter('fates_allom_agb3'))
parms['agb4'] = CDLParse(args.fnamein,parameter('fates_allom_agb4'))
parms['d2bl1'] = CDLParse(args.fnamein,parameter('fates_allom_d2bl1'))
parms['d2bl2'] = CDLParse(args.fnamein,parameter('fates_allom_d2bl2'))
parms['d2bl3'] = CDLParse(args.fnamein,parameter('fates_allom_d2bl3'))
parms['wood_density'] = CDLParse(args.fnamein,parameter('fates_wood_density'))
parms['c2b'] = CDLParse(args.fnamein,parameter('fates_c2b'))
parms['la_per_sa_int'] = CDLParse(args.fnamein,parameter('fates_allom_la_per_sa_int'))
parms['la_per_sa_slp'] = CDLParse(args.fnamein,parameter('fates_allom_la_per_sa_slp'))
parms['slatop'] = CDLParse(args.fnamein,parameter('fates_leaf_slatop'))
parms['slamax'] = CDLParse(args.fnamein,parameter('fates_leaf_slamax'))
parms['l2fr'] = CDLParse(args.fnamein,parameter('fates_allom_l2fr'))
parms['agb_frac'] = CDLParse(args.fnamein,parameter('fates_allom_agb_frac'))
parms['blca_expnt_diff'] = CDLParse(args.fnamein,parameter('fates_allom_blca_expnt_diff'))
parms['d2ca_coeff_min'] = CDLParse(args.fnamein,parameter('fates_allom_d2ca_coefficient_min'))
parms['d2ca_coeff_max'] = CDLParse(args.fnamein,parameter('fates_allom_d2ca_coefficient_max'))
parms['sai_scaler'] = CDLParse(args.fnamein,parameter('fates_allom_sai_scaler'))

# Read in the parameters that are not necessary for the F90 allometry algorithms,
# but are useful for these scripts (e.g. the name of the parameter, and minimum height)
# =======================================================================================

eparms = {}
eparms['recruit_hgt_min'] = CDLParse(args.fnamein,parameter('fates_recruit_hgt_min'))
eparms['name'] = CDLParse(args.fnamein,parameter('fates_pftname'))
eparms['vcmax25top'] = CDLParse(args.fnamein,parameter('fates_leaf_vcmax25top'))


# Determine how many PFTs are here, also check to make sure that all parameters
# have the same number
# =======================================================================================
numpft=-1
for key, parm in parms.items():
    if( (len(parm.vals) == numpft) or (numpft==-1) ):
        numpft=len(parm.vals)
    else:
        print('Bad length in PFT parameter')
        print('parameter: {}, vals:'.format(parm.symbol),parm.vals)


# ==============================================================================
# Allocate fortran PFT arrays
# ==============================================================================

iret=f90_pftalloc(byref(c_int(numpft)))

# ==============================================================================
# Populate the Fortran PFT structure
# ==============================================================================

# First set the arg types
f90_pftset.argtypes = \
    [POINTER(c_int),POINTER(c_double),POINTER(c_int),c_char_p,c_long]

for ipft in range(numpft):
    for key, parm in parms.items():
        #print 'py: sending to F90: {0} = {1}'.format(parm.symbol,parm.vals[ipft])
        iret=f90_pftset(c_int(ipft+1), \
                       c_double(parm.vals[ipft]), \
                       c_int(0), \
                       c_char_p(parm.symbol), \
                       c_long(len(parm.symbol)))


# =========================================================================
# Initialize Output Arrays
# =========================================================================

blmaxi  = np.zeros((numpft,ndbh))
blmaxd  = np.zeros((numpft,ndbh))
bfrmax = np.zeros((numpft,ndbh))
hi     = np.zeros((numpft,ndbh))
hd     = np.zeros((numpft,ndbh))
bagwi   = np.zeros((numpft,ndbh))
bagwd   = np.zeros((numpft,ndbh))
dbh    = np.zeros((numpft,ndbh))
bbgw    = np.zeros((numpft,ndbh))
bsapi  = np.zeros((numpft,ndbh))
bsapd  = np.zeros((numpft,ndbh))
asapd  = np.zeros((numpft,ndbh))
bstore = np.zeros((numpft,ndbh))
bdead  = np.zeros((numpft,ndbh))
dbhe   = np.zeros((numpft,ndbh))
camin  = np.zeros((numpft,ndbh))
ldense = np.zeros((numpft,ndbh))
treelai = np.zeros((numpft,ndbh))
blmax_o_dbagwdh = np.zeros((numpft,ndbh))
blmax_o_dbagwdd = np.zeros((numpft,ndbh))


for ipft in range(numpft):

    print 'py: Solving for pft: {}'.format(ipft+1)

    # Initialize Height   #(d,ipft,h,dhdd)
    ch_min = c_double(eparms['recruit_hgt_min'].vals[ipft])

    cd     = c_double(-9.0)
    cdddh  = c_double(-9.0)
    cipft  = c_int(ipft+1)
    cinit  = c_int(0)

    # Calculate the minimum dbh
    iret=f90_h2d(byref(ch_min),byref(cipft),byref(cd),byref(cdddh))

    # Generate a vector of diameters (use dbh)
    dbh[ipft,:] = np.linspace(cd.value,maxdbh,num=ndbh)

    # Initialize various output vectors
    cd = c_double(dbh[ipft,0])
    ch = c_double(-9.0)
    cdhdd = c_double(-9.0)
    cbagw = c_double(-9.0)
    cdbagwdd = c_double(-9.0)
    cblmax = c_double(-9.0)
    cdblmaxdd = c_double(-9.0)
    cbfrmax = c_double(-9.0)
    cdbfrmaxdd = c_double(-9.0)
    cbbgw = c_double(-9.0)
    cdbbgwdd = c_double(-9.0)
    cbsap = c_double(-9.0)
    cdbsapdd = c_double(-9.0)
    cbdead = c_double(-9.0)
    cdbdeaddd = c_double(-9.0)
    ccamin = c_double(-9.0)
    casapw = c_double(-9.0)   # Sapwood area
    cbstore = c_double(-9.0)
    cdbstoredd = c_double(-9.0)

    iret=f90_h(byref(cd),byref(cipft),byref(ch),byref(cdhdd))
    hi[ipft,0] = ch.value
    hd[ipft,0] = ch.value
    print 'py: initialize h[{},0]={}'.format(ipft+1,ch.value)

    # Initialize AGB      #(d,ipft,bagw,dbagwdd)
    iret=f90_bagw(byref(cd),byref(cipft),byref(cbagw),byref(cdbagwdd))
    bagwi[ipft,0] = cbagw.value
    print 'py: initialize bagwi[{},0]={}'.format(ipft+1,cbagw.value)

    # Initialize bleaf    #(d,ipft,canopy_trim,bl,dbldd)
    iret=f90_bleaf(byref(cd),byref(cipft),byref(ccanopy_trim),byref(cblmax),byref(cdblmaxdd))
    blmaxi[ipft,0] = cblmax.value
    blmaxd[ipft,0] = cblmax.value
    print 'py: initialize blmaxi[{},0]={}'.format(ipft+1,cblmax.value)

    # Initialize bstore #(d,ipft,canopy_trim,bstore,dbstoredd)
    iret=f90_bstore(byref(cd),byref(cipft),byref(ccanopy_trim),byref(cbstore),byref(cdbstoredd))
    bstore[ipft,0] = cbstore.value

    # calculate crown area (d,nplant,site_spread,ipft,c_area)  Using nplant = 1, generates units of m2
    #  spread is likely 0.0, which is the value it tends towards when canopies close
    #              (dbh,      nplant,        site_spread,         ipft,       c_area,inverse)
    iret= f90_carea(byref(cd),byref(cnplant),byref(csite_spread),byref(cipft),byref(ccamin),byref(cdo_reverse))
    camin[ipft,0]  = ccamin.value
    ldense[ipft,0] = blmaxi[ipft,0]/camin[ipft,0]
    print 'py: initialize careai[{},0]={}'.format(ipft+1,ccamin.value)

    #f90_treelai(leaf_c, pft, c_area, nplant, cl, canopy_lai, vcmax25top)
    cvcmax=c_double(eparms['vcmax25top'].vals[ipft])
    treelai[ipft,0]=f90_treelai(byref(cblmax),byref(cipft),byref(ccamin), \
                                byref(cnplant),byref(cilayer),byref(ccanopy_lai),byref(cvcmax))

    # Initialize fine roots  #(d,ipft,canopy_trim,bfr,dbfrdd)
    iret=f90_bfineroot(byref(cd),byref(cipft),byref(ccanopy_trim), \
                       byref(cbfrmax),byref(cdbfrmaxdd))
    bfrmax[ipft,0] = cbfrmax.value
    print 'py: initialize bfrmax[{},0]={}'.format(ipft+1,cbfrmax.value)

    # Initialize coarse roots #(d,ipft,bbgw,dbbgwdd)
    iret=f90_bbgw(byref(cd),byref(cipft),byref(c_double(1.0)), \
                 byref(cbbgw),byref(cdbbgwdd))
    bbgw[ipft,0] = cbbgw.value
    print 'py: initialize bbgw[{},0]={}'.format(ipft+1,cbbgw.value)


    # Initialize bsap  (d,ipft,canopy_trim,asapw,bsap,dbsapdd)
    iret=f90_bsap(byref(cd),byref(cipft),byref(ccanopy_trim),byref(casapw),byref(cbsap),byref(cdbsapdd))
    bsapi[ipft,0] = cbsap.value
    bsapd[ipft,0] = cbsap.value
    asapd[ipft,0] = casapw.value
    print 'py: initialize bsapi[{},0]={}'.format(ipft+1,cbsap.value)

    # bdead #(bagw,bbgw,bsap,ipft,bdead,dbagwdd,dbbgwdd,dbsapdd,dbdeaddd)
    iret=f90_bdead(byref(cbagw),byref(cbbgw),byref(cbsap),byref(cipft), \
                   byref(cbdead),byref(cdbagwdd),byref(cdbbgwdd), \
                   byref(cdbsapdd),byref(cdbdeaddd))

    bdead[ipft,0] = cbdead.value
    print 'py: initialize bdead[{},0]={}'.format(ipft+1,cbdead.value)

    # the metric that shan't be spoken
    blmax_o_dbagwdh[ipft,0]  = blmaxi[ipft,0]/(cdbagwdd.value/cdhdd.value)

    # the metric that shan't be spoken
    blmax_o_dbagwdd[ipft,0]  = blmaxi[ipft,0]/(cdbagwdd.value)

    for idi in range(1,ndbh):

        dp = dbh[ipft,idi-1]  # previous position
        dc = dbh[ipft,idi]    # current position
        dd = dc-dp

        cdp = c_double(dp)
        cdc = c_double(dc)
        cdbhe = c_double(-9.0)
        cddedh = c_double(-9.0)

        if(ipft==2):
            print("===")

        # integrate height  #(d,ipft,h,dhdd)
        iret=f90_h(byref(cdc),byref(cipft),byref(ch),byref(cdhdd))
        hi[ipft,idi] = hi[ipft,idi-1] + cdhdd.value*dd

        # diagnosed height
        hd[ipft,idi] = ch.value

        # diagnose AGB  #(d,h,ipft,bagw,dbagwdd)
        iret=f90_bagw(byref(cdc),byref(cipft),byref(cbagw),byref(cdbagwdd))
        bagwd[ipft,idi] = cbagw.value

        # integrate AGB #(d,h,ipft,bagw,dbagwdd)
        iret=f90_bagw(byref(cdp),byref(cipft),byref(cbagw),byref(cdbagwdd))
        bagwi[ipft,idi] = bagwi[ipft,idi-1] + cdbagwdd.value*dd

        # diagnose bleaf #(d,ipft,blmax,dblmaxdd)
        iret=f90_bleaf(byref(cdc),byref(cipft),byref(c_double(1.0)),byref(cblmax),byref(cdblmaxdd))
        blmaxd[ipft,idi] = cblmax.value

        # bstore #(d,ipft,canopy_trim,bstore,dbstoredd)
        iret=f90_bstore(byref(cdc),byref(cipft),byref(ccanopy_trim),byref(cbstore),byref(cdbstoredd))
        bstore[ipft,idi] = cbstore.value

        # calculate crown area (d,nplant,site_spread,ipft,c_area)  Using nplant = 1, generates units of m2
        iret= f90_carea(byref(cdc),byref(cnplant),byref(csite_spread),byref(cipft),byref(ccamin),byref(cdo_reverse))
        camin[ipft,idi]  = ccamin.value

        #f90_treelai(leaf_c, pft, c_area, nplant, cl, canopy_lai, vcmax25top)
        cvcmax=c_double(eparms['vcmax25top'].vals[ipft])
        treelai[ipft,idi]=f90_treelai(byref(cblmax),byref(cipft),byref(ccamin), \
                                byref(cnplant),byref(cilayer),byref(ccanopy_lai),byref(cvcmax))

        # integrate bleaf #(d,ipft,blmax,dblmaxdd)
        iret=f90_bleaf(byref(cdp),byref(cipft),byref(c_double(1.0)),byref(cblmax),byref(cdblmaxdd))
        blmaxi[ipft,idi] = blmaxi[ipft,idi-1] + cdblmaxdd.value*dd

        # leaf mass per square meter of crown
        ldense[ipft,idi] = blmaxd[ipft,idi]/camin[ipft,idi]

        # integrate bfineroot #(d,ipft,canopy_trim,bfr,dbfrdd)
        iret=f90_bfineroot(byref(cdp),byref(cipft),byref(c_double(1.0)),byref(cbfrmax),byref(cdbfrmaxdd))
        bfrmax[ipft,idi] = bfrmax[ipft,idi-1] + cdbfrmaxdd.value*dd

        # integrate bbgw #(d,h,ipft,bbgw,dbbgwdd)
        iret=f90_bbgw(byref(cdp),byref(cipft),byref(cbbgw),byref(cdbbgwdd))
        bbgw[ipft,idi] = bbgw[ipft,idi-1] + cdbbgwdd.value*dd

        # diagnose bsap  # (d,ipft,canopy_trim,asapw,bsap,dbsapdd)
        iret=f90_bsap(byref(cdc),byref(cipft),byref(ccanopy_trim),byref(casapw),byref(cbsap),byref(cdbsapdd))
        bsapd[ipft,idi] = cbsap.value   # Biomass
        asapd[ipft,idi] = casapw.value  # Area

        # integrate bsap
        iret=f90_bsap(byref(cdp),byref(cipft),byref(ccanopy_trim),byref(casapw),byref(cbsap),byref(cdbsapdd))
        bsapi[ipft,idi] = bsapi[ipft,idi-1] + cdbsapdd.value*dd

        # the metric that shan't be spoken
        # previous t-step derivatives are used for simplicity
        if cdhdd.value<0.000001:
            blmax_o_dbagwdh[ipft,idi] = None
        else:
            blmax_o_dbagwdh[ipft,idi] = blmaxi[ipft,idi-1]/(cdbagwdd.value/cdhdd.value)

        # the metric that shan't be spoken
        # previous t-step derivatives are used for simplicity
        blmax_o_dbagwdd[ipft,idi]  = blmaxi[ipft,idi-1]/(cdbagwdd.value)

        # Diagnose bdead (bagw,bbgw,bsap,ipft,bdead,dbagwdd,dbbgwdd,dbsapdd,dbdeaddd)

        iret=f90_bdead(byref(c_double(bagwi[ipft,idi])), \
                       byref(c_double(bbgw[ipft,idi])), \
                       byref(c_double(bsapi[ipft,idi])), \
                       byref(cipft), byref(cbdead), \
                       byref(cdbagwdd),byref(cdbbgwdd), \
                       byref(cdbsapdd),byref(cdbdeaddd))
        bdead[ipft,idi] = cbdead.value


# Create the appropriate number of line-styles, colors and widths
linestyles_base = ['-', '--', '-.', ':']
linestyles=[]
for i in range(int(math.floor(float(numpft)/float(len(linestyles_base))))):
    linestyles.extend(linestyles_base)
for i in range(numpft-len(linestyles)):
    linestyles.append(linestyles_base[i])

my_colors = DiscreteCubeHelix(numpft)


mp.rcParams.update({'font.size': 14})
mp.rcParams["savefig.directory"] = ""    #os.chdir(os.path.dirname(__file__))

legfs = 12
lwidth = 2.0

#code.interact(local=dict(globals(), **locals()))

if(True):
    fig0   = plt.figure()
    figleg = plt.figure()
    ax = fig0.add_subplot(111)
    ax.axis("off")
    ax.set_axis_off()
    proxies = ()
    for ipft in range(numpft):
        proxies = proxies + (mp.lines.Line2D([],[], \
                                             linestyle=linestyles[ipft], \
                                             color=my_colors(ipft), \
                                             label=eparms['name'].vals[ipft], \
                                             linewidth=lwidth),)
    figleg.legend(handles=proxies,fontsize=12,frameon=False,labelspacing=0.25,loc='center')
    plt.show(block=False)
    plt.close(fig0)


if(True):
    fig1 = plt.figure()
    figleg = plt.figure()
    for ipft in range(numpft):
        plt.plot(dbh[ipft,:],hi[ipft,:],linestyle=linestyles[ipft],color=my_colors(ipft),linewidth=lwidth)
    plt.xlabel('diameter [cm]')
    plt.ylabel('height [m]')
    plt.title('Integrated Heights')
    plt.grid(True)
    plt.tight_layout()

if(False):
    fig1_0 = plt.figure()
    for ipft in range(numpft):
        plt.plot(dbh[ipft,0:15],hi[ipft,0:15],linestyle=linestyles[ipft],color=my_colors(ipft),linewidth=lwidth)
    plt.xlabel('diameter [cm]')
    plt.ylabel('height [m]')
    plt.title('Integrated Heights')
    plt.grid(True)
    plt.tight_layout()

if(False):
    fig1_1 = plt.figure()
    for ipft in range(numpft):
        plt.plot(hd[ipft,:],hi[ipft,:],linestyle=linestyles[ipft],color=my_colors(ipft),linewidth=lwidth)
    plt.xlabel('height (diagnosed) [m]')
    plt.ylabel('height (integrated) [m]')
    plt.title('Height')
    plt.grid(True)
    plt.savefig("plots/hdhi.png")

if(False):
    fig2=plt.figure()
    for ipft in range(numpft):
        plt.plot(blmaxd[ipft,:],blmaxi[ipft,:],linestyle=linestyles[ipft],color=my_colors(ipft),linewidth=lwidth)
    plt.xlabel('diagnosed [kgC]')
    plt.ylabel('integrated [kgC]')
    plt.title('Maximum Leaf Biomass')
    plt.grid(True)
    plt.tight_layout()

if(True):
    fig3=plt.figure()
    for ipft in range(numpft):
        plt.plot(dbh[ipft,:],blmaxi[ipft,:],linestyle=linestyles[ipft],color=my_colors(ipft),linewidth=lwidth)
    plt.xlabel('diameter [cm]')
    plt.ylabel('mass [kgC]')
    plt.title('Maximum Leaf Biomass')
    plt.grid(True)
    plt.tight_layout()

if(True):
    fig3_1=plt.figure()
    for ipft in range(numpft):
        plt.plot(dbh[ipft,1:15],blmaxi[ipft,1:15],linestyle=linestyles[ipft],color=my_colors(ipft),linewidth=lwidth)
    plt.xlabel('diameter [cm]')
    plt.ylabel('mass [kgC]')
    plt.title('Maximum Leaf Biomass (saplings)')
    plt.grid(True)
    plt.tight_layout()


if(True):
    fig4=plt.figure()
    for ipft in range(numpft):
        plt.plot(dbh[ipft,:],camin[ipft,:],linestyle=linestyles[ipft],color=my_colors(ipft),linewidth=lwidth)
    plt.xlabel('diameter [cm]')
    plt.ylabel('[m2] (closed canopy)')
    plt.title('Crown Area')
    plt.grid(True)
    plt.tight_layout()

if(True):
    fig4_1=plt.figure()
    for ipft in range(numpft):
        plt.plot(dbh[ipft,:],ldense[ipft,:],linestyle=linestyles[ipft],color=my_colors(ipft),linewidth=lwidth)
    plt.xlabel('diameter [cm]')
    plt.ylabel('[kgC/m2] (closed canopy)')
    plt.title('Leaf Mass Per Crown Area')
    plt.grid(True)
    plt.tight_layout()


if(True):
    fig6=plt.figure()
    for ipft in range(numpft):
        plt.plot(dbh[ipft,:],bagwi[ipft,:]/1000,linestyle=linestyles[ipft],color=my_colors(ipft),linewidth=lwidth)
    plt.xlabel('diameter [cm]')
    plt.ylabel('AGB [MgC]')
    plt.title('Above Ground Biomass')
    plt.grid(True)
    plt.tight_layout()

if(False):
    fig6_1=plt.figure()
    for ipft in range(numpft):
        plt.plot(bagwd[ipft,:]/1000,bagwi[ipft,:]/1000,linestyle=linestyles[ipft],color=my_colors(ipft),linewidth=lwidth)
    plt.xlabel('AGBW deterministic [MgC]')
    plt.ylabel('AGBW integrated [MgC]')
    plt.title('Above Ground Biomass')
    plt.grid(True)
    plt.tight_layout()

if(False):
    fig5=plt.figure()
    for ipft in range(numpft):
        gpmask  = np.isfinite(blmax_o_dbagwdh[ipft,:])
        plt.plot(dbh[ipft,gpmask],blmax_o_dbagwdh[ipft,gpmask],linestyle=linestyles[ipft],color=my_colors(ipft),linewidth=lwidth)
    plt.xlabel('diameter [cm]')
    plt.ylabel('growth potential: bl/(dAGB/dh) [m]')
    plt.title('Height Growth Potential')
    plt.grid(True)
    plt.tight_layout()

if(False):
    fig6=plt.figure()
    for ipft in range(numpft):
        plt.plot(dbh[ipft,:],blmax_o_dbagwdd[ipft,:],linestyle=linestyles[ipft],color=my_colors(ipft),linewidth=lwidth)
    plt.xlabel('diameter [cm]')
    plt.ylabel('growth potential: bl/(dAGB/dd) [cm]')
    plt.title('Diameter Growth Potential')
    plt.grid(True)
    plt.tight_layout()

if(False):
    fig7=plt.figure()
    for ipft in range(numpft):
        plt.plot(bsapd[ipft,:],bsapi[ipft,:],linestyle=linestyles[ipft],color=my_colors(ipft),linewidth=lwidth)
    plt.xlabel('deterministic [kgC]')
    plt.ylabel('integrated [kgC]')
    plt.title('Sapwood Biomass')
    plt.grid(True)
    plt.tight_layout()

if(False):
    fig7_0=plt.figure()
    for ipft in range(numpft):
        plt.plot(dbh[ipft,:],bsapd[ipft,:],linestyle=linestyles[ipft],color=my_colors(ipft),linewidth=lwidth)
    plt.xlabel('Diameter [cm]')
    plt.ylabel('[kgC]')
    plt.title('Sapwood Biomass')
    plt.grid(True)
    plt.tight_layout()

if(True):
    fig7_2=plt.figure(figsize=(8,6))
    # Sapwood
    ax = fig7_2.add_subplot(221)
    for ipft in range(numpft):
        ax.plot(dbh[ipft,:],bsapd[ipft,:]/(bsapd[ipft,:]+blmaxi[ipft,:]+bfrmax[ipft,:]+bstore[ipft,:]), \
                 linestyle=linestyles[ipft],color=my_colors(ipft),linewidth=lwidth)
    ax.set_xlabel('diameter [cm]')
    ax.set_ylabel('[kgC/kgC]')
    ax.set_title('Sapwood (fraction of total live)')
    ax.grid(True)
    # Leaf
    ax = fig7_2.add_subplot(222)
    for ipft in range(numpft):
        ax.plot(dbh[ipft,:],blmaxi[ipft,:]/(bsapd[ipft,:]+blmaxi[ipft,:]+bfrmax[ipft,:]+bstore[ipft,:]), \
                 linestyle=linestyles[ipft],color=my_colors(ipft),linewidth=lwidth)
    ax.set_xlabel('diameter [cm]')
    ax.set_ylabel('[kgC/kgC]')
    ax.set_title('Leaf (fraction of total live)')
    ax.grid(True)
    # Fine Root
    ax = fig7_2.add_subplot(223)
    for ipft in range(numpft):
        ax.plot(dbh[ipft,:],bfrmax[ipft,:]/(bsapd[ipft,:]+blmaxi[ipft,:]+bfrmax[ipft,:]+bstore[ipft,:]), \
                 linestyle=linestyles[ipft],color=my_colors(ipft),linewidth=lwidth)
    ax.set_xlabel('diameter [cm]')
    ax.set_ylabel('[kgC/kgC]')
    ax.set_title('Fine-Root (fraction of total live)')
    ax.grid(True)
    # Storage
    ax = fig7_2.add_subplot(224)
    for ipft in range(numpft):
        ax.plot(dbh[ipft,:],bstore[ipft,:]/(bsapd[ipft,:]+blmaxi[ipft,:]+bfrmax[ipft,:]+bstore[ipft,:]), \
                 linestyle=linestyles[ipft],color=my_colors(ipft),linewidth=lwidth)
    ax.set_xlabel('diameter [cm]')
    ax.set_ylabel('[kgC/kgC]')
    ax.set_title('Storage (fraction of total live)')
    ax.grid(True)

    plt.tight_layout()



if(True):
    fig8=plt.figure()
    for ipft in range(numpft):
        plt.semilogy(dbh[ipft,:],treelai[ipft,:],linestyle=linestyles[ipft],color=my_colors(ipft),linewidth=lwidth)
    plt.xlabel('diameter [cm]')
    plt.ylabel('[m2/m2]')
    plt.title('In-Crown LAI')
    plt.grid(True)
    plt.tight_layout()


#    print(blmaxi[2,:])
#    print(bfrmax[2,:])
#    print(bstore[2,:])
#    print(bsapd[2,:])

plt.show()
