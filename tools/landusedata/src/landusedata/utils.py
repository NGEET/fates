import xarray

def RegridTargetPrep(regrid_target):
    """
    Rename the CLM surface data file regrid target lat/lon dimensions

    The CLM surface datasets use a dimensions name that is not recognizable
    by xesmf.  As such, this function renames the dimensions.  It also adds
    lat/lon coordinates based on the LONGXY and LATIXY variables.
    """
    ds = regrid_target.rename_dims(dims_dict={'lsmlat':'lat','lsmlon':'lon'})
    ds['lon'] = ds.LONGXY.isel(lat=0)
    ds['lat'] = ds.LATIXY.isel(lon=0)

    return ds
