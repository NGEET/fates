import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.dates as mdates
import sys
import code  # For development: code.interact(local=locals())  code.interact(local=dict(globals(), **locals()))
import argparse
import math
from scipy.io import netcdf as nc
import xml.etree.ElementTree as et

# The user specifies a couplet of domain/surface files
# from which they want to base their new unstructured grid. Then they provide
# a list of geographic coordinates in latitude and longitude. These coordinates
# are sampled from the base dataset, and written to an output dataset.
#
# This method is certainly useful for generating small sets of unstructured
# grid-cells. It may not be the best method for generating large sets. One
# may want to use "ncks" (ie Charlie Zender's NCO tools) for subsetting large
# grids. This method will arrange the new grid cells in a 1D vector, and assumes
# the input grids are based on a 2d array of cells.
#
# This method will assume that the grid-cell extents of the new unstructured grids
# match the extents of the base files. If you want finer or coarser resolution,
# just dig up a different base file.
#
# This method uses nearest neighbor.
#
# This may have trouble on newer surface files, particularly if they have topo
# unit information.  It won't be difficult to add that type of functionality
# if it doesn't work, I (Ryan) just haven't tried it.
#
# All controls over this process can be found in the xml control file. See
# andes7x7.xml for an example.

# Usage MakeUnstructGrid.py --fin=xmlfile.xml

def TransferData(da_key,ds_base,ds_unst,minis,minjs,dset_type):
        
    print('  Transferring: {}'.format(da_key))

    if(dset_type=='domain'):
        xname = 'nj'
        yname = 'ni'
        ny = len(minis)  #nj = len(minis)
        nx = 1           #ni
        nv = 4
    elif(dset_type=='surface'):
        xname = 'lsmlat'
        yname = 'lsmlon'
        ny = len(minis)
        nx = 1
        
        #numurbl = 3 ;
        #nlevurb = 5 ;
        #numrad = 2 ;
        #nchar = 256 ;
        #nlevsoi = 10 ;
        #time = UNLIMITED ; // (12 currently)
        #lsmpft = 17 ;
        #natpft = 17 ;


    
    # Determin the data type
    if(ds_base[da_key].dtype == 'float64'):
        dtype_out = np.float64
    elif(ds_base[da_key].dtype == 'int32'):
        dtype_out = np.int32
    else:
        print('unknown data type: {}.\n Exiting'.format(ds_base[da_key].dtype))
        exit(2)

    

    # The lat-lon is always the last two dimensions
    # Time is always the first dimension
    
    # Check to see if this has spatial dimensions  
    dimlist = list(ds_base[da_key].dims)

    # 2D (nj,ni)
    #    (lsmlat,lsmlon)
    # 3D (nj, ni, nv)
    #    (x,  lsmlat, lsmlon)
    # 4D (x,y,lsmlat, lsmlon)

    if(any([dim==xname for dim in dimlist]) and (len(dimlist)==2)):

        # This is 2D and they are or use geographic coordinates

        ds_unst[da_key] = \
            xr.DataArray(np.empty((nx,ny), dtype=dtype_out),dims=dimlist)

        for k in range(len(minis)):
            i = minis[k]
            j = minjs[k]
            ds_unst[da_key].loc[0,k] = ds_base[da_key].data[j,i]

    elif(any([dim==xname for dim in dimlist]) and any([dim=='nv' for dim in dimlist]) ):

        # This is the 3D coordinate in the domain file for vertices

        ds_unst[da_key] = \
            xr.DataArray(np.empty((nx,ny,nv), dtype=dtype_out),dims=dimlist)
        for k in range(len(minis)):
            i = minis[k]
            j = minjs[k]
            ds_unst[da_key].loc[0,k,:] = ds_base[da_key].data[j,i,:]
                
    elif(any([dim==xname for dim in dimlist]) and len(dimlist)==3):

        # This has dim==3, surface file and contains coordinates (x,lsmlat, lsmlon)

        dim0 = ds_base.dims[dimlist[0]]
        ds_unst[da_key] = \
            xr.DataArray(np.empty((dim0,nx,ny), dtype=dtype_out),dims=dimlist)
        for k in range(len(minis)):
            i = minis[k]
            j = minjs[k]
            ds_unst[da_key].loc[:,0,k] = ds_base[da_key].data[:,j,i]
            
    elif(any([dim==xname for dim in dimlist]) and len(dimlist)==4):

        # This has dim==4, surface file and contains coordinates (x,y,lsmlat, lsmlon)

        dim0 = ds_base.dims[dimlist[0]]
        dim1 = ds_base.dims[dimlist[1]]
        ds_unst[da_key] = \
            xr.DataArray(np.empty((dim0,dim1,nx,ny), dtype=dtype_out),dims=dimlist)
        for k in range(len(minis)):
            i = minis[k]
            j = minjs[k]
            ds_unst[da_key].loc[:,:,0,k] = ds_base[da_key].data[:,:,j,i]

    elif(len(dimlist)==1):   
        # If there is no spatial coordinate, then just copy over what is there
        #dimsizes = tuple([ds_base.dims[txt] for txt in dimlist])
        #ds_unst[da_key] = \
        #    xr.DataArray(np.empty(dimsizes,dtype=dtype_out),dims=dimlist)
        #ds_unst[da_key].loc[:] = ds_base[da_key].data[:]
        ds_unst[da_key] = ds_base[da_key]
    else:
        # If there is no spatial coordinate, then just copy over what is there

        ds_unst[da_key] = ds_base[da_key]
        

    # Once the new dataarray is created, transfer over metadata from original
    ds_unst[da_key].attrs = ds_base[da_key].attrs

    return(ds_unst)

    
def main(argv):

    parser = argparse.ArgumentParser(description='Parse command line arguments to this script.')
    parser.add_argument('--fin', dest='xmlfile', type=str, help="path to the xml control file",required=True)
    args = parser.parse_args()

    xmlroot = et.parse(args.xmlfile).getroot()

    print('  --------------------------------------------------------------  ')
    print('\n          Creating a new domain/surface couplet \n')
    print('  --------------------------------------------------------------\n  ')
    
    # Get the domain base name
    try:
        domain_base = xmlroot.find('domain_base').text.strip()
        domain_base_file = domain_base.split('/')[-1]
        
    except:
        print('Could not find xml entry: {}'.format('domain_base'))
        exit(2)

    # Get the new unstructured domain name (ie output)
    try:
        domain_unst = xmlroot.find('domain_unst').text.strip()
    except:
        print('Could not find xml entry: {}'.format('domain_unst'))
        exit(2)


    # Get the surface base name
    try:
        surface_base = xmlroot.find('surface_base').text.strip()
        surface_base_file = surface_base.split('/')[-1]
        
    except:
        print('Could not find xml entry: {}'.format('surface_base'))
        exit(2)

    # Get the new unstructured surface name (ie output)
    try:
        surface_unst = xmlroot.find('surface_unst').text.strip()
    except:
        print('Could not find xml entry: {}'.format('surface_unst'))
        exit(2)


        
    # Get a list of lon coordinates (force them into 0-360 convention)
    try:
        lon_subset_text = xmlroot.find('lon_list').text.strip().split(',')
        lon_subset = []
        for txt in lon_subset_text:
            lon = float(txt)
            if(lon<0.0):
                lon = 360.0+lon
            lon_subset.append(lon)

    except:
        print('Could not find xml entry: {}'.format('lon_list'))
        exit(2)
        
    # Get a list of lat coordinates
    try:
        lat_subset_text = xmlroot.find('lat_list').text.strip().split(',')
        lat_subset = [float(txt) for txt in lat_subset_text]

    except:
        print('Could not find xml entry: {}'.format('lat_list'))
        exit(2)

    # Check to make sure that the lat and lons are same length

    if( len(lat_subset) != len(lon_subset) ):
        print('number of latitude subset points must match number of longitude subsets')
        exit(2)
    else:
        nj = len(lat_subset)
        print('  Found N={} lat/lon coordinates'.format(nj))
        
    
    #code.interact(local=dict(globals(), **locals()))
    # ------------------------------------------------------------------------------
    # >>> ds_domain_base.data_vars
    # Data variables:
    #xv       (nj, ni, nv) float64 358.8 1.25 1.25 358.8 ... 358.7 358.7 356.2
    #yv       (nj, ni, nv) float64 -90.0 -90.0 -89.05 -89.05 ... 89.05 90.0 90.0
    #mask     (nj, ni) int32 1 1 1 1 1 1 1 1 1 1 1 1 ... 0 0 0 0 0 0 0 0 0 0 0 0
    #area     (nj, ni) float64 5.964e-06 5.964e-06 ... 5.964e-06 5.964e-06
    #frac     (nj, ni) float64 1.0 1.0 1.0 1.0 1.0 1.0 ... 0.0 0.0 0.0 0.0 0.0
    # >>> ds_domain_base.coords
    #Coordinates:
    #    xc       (nj, ni) float64 0.0 2.5 5.0 7.5 10.0 ... 350.0 352.5 355.0 357.5
    #    yc       (nj, ni) float64 -90.0 -90.0 -90.0 -90.0 ... 90.0 90.0 90.0 90.0
    # ------------------------------------------------------------------------------

    lon_subset_faux_2d = np.reshape(lon_subset, (-1, 1))
    lat_subset_faux_2d = np.reshape(lat_subset, (-1, 1))
    
    ds_domain_base = xr.open_dataset(domain_base)
    
    # Lets find the indices for xc and yc that most closely match our coordinates
    minis = []
    minjs = []
    for j in range(nj):
        lat = lat_subset[j]
        lon = lon_subset[j]
        delt = (ds_domain_base['xc'].data-lon)**2.0 +  (ds_domain_base['yc'].data-lat)**2.0
        minj,mini = np.unravel_index(delt.argmin(), delt.shape)
        minis.append(mini)
        minjs.append(minj)

    # Domain Processing
    # ===========================================================================================
    # Initialize the new dataset
    ds_domain_unst = xr.Dataset(
        attrs=ds_domain_base.attrs,        
    )
    ds_domain_unst.attrs["modification"]="Modified with SurfToVec.py, based on {}.".format(domain_base_file)

    #ode.interact(local=dict(globals(), **locals()))
    
    # Loop through existing datasets, allocate new arrays
    # and transfer over point data
    for da_key in ds_domain_base.data_vars:
        ds_domain_unst = TransferData(da_key,ds_domain_base,ds_domain_unst,minis,minjs,'domain')

    for da_key in ds_domain_base.coords:
        ds_domain_unst = TransferData(da_key,ds_domain_base,ds_domain_unst,minis,minjs,'domain')
        
    print('\n  Writing: {}'.format(domain_unst))
    ds_domain_unst.to_netcdf(domain_unst)  #,mode='a')

    # Surface Processing
    # ===========================================================================================
    
    ds_surface_base = xr.open_dataset(surface_base)

    # Initialize the new dataset
    ds_surface_unst = xr.Dataset(
        attrs=ds_surface_base.attrs,        
    )
    ds_surface_unst.attrs["modification"]="Modified with SurfToVec.py, based on {}.".format(surface_base_file)

    
    # Loop through existing datasets, allocate new arrays
    # and transfer over point data
    for da_key in ds_surface_base.data_vars:
        ds_surface_unst = TransferData(da_key,ds_surface_base,ds_surface_unst,minis,minjs,'surface')

    for da_key in ds_surface_base.coords:
        ds_surface_unst = TransferData(da_key,ds_surface_base,ds_surface_unst,minis,minjs,'surface')
        
    print('\n  Writing: {}'.format(surface_unst))
    ds_surface_unst.to_netcdf(surface_unst)  #,mode='a')
    
    print('\n')

# This is the actual call to main
if __name__ == "__main__":
    main(sys.argv)
