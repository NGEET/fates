import xarray as xr

def ImportRegridTarget(filename):
    dataset = xr.open_dataset(filename)

    # Check the file type
    dim_list = list(dataset.dims)
    if ('lsmlat' in list(dataset.dims)) != True:
        raise TypeError("incorrect file, must be surface dataset")

    # Prepare the the regrid dataset
    dataset = _RegridTargetPrep(dataset)

    return dataset

def _RegridTargetPrep(regrid_target):
    """
    Rename the CLM surface data file regrid target lat/lon dimensions

    The CLM surface datasets use a dimensions name that is not recognizable
    by xesmf.  As such, this function renames the dimensions.  It also adds
    lat/lon coordinates based on the LONGXY and LATIXY variables.
    """
    regrid_target = regrid_target.rename_dims(dims_dict={'lsmlat':'lat','lsmlon':'lon'})
    regrid_target['lon'] = regrid_target.LONGXY.isel(lat=0)
    regrid_target['lat'] = regrid_target.LATIXY.isel(lon=0)

    return regrid_target

# TODO: add the follow common functions
# - write to files
# - read dataset into list of datasets
