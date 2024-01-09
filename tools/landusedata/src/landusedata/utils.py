import xarray as xr

def RenameLatLon(regrid_target):
    """
    Rename the CLM surface data file regrid target lat/lon dimensions

    The CLM surface datasets use a dimensions name that is not recognizable
    by xesmf.  As such, this function renames the dimensions.
    """
    ds = regrid_target.rename_dims(dims_dict={'lsmlat':'lat','lsmlon':'lon'})
    return ds
