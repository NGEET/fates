"""functions for using fates and xarray"""
import xarray as xr
import numpy as np


def _get_check_dim(dim_short, dataset):
    """Get dim name from short code and ensure it's on Dataset

    Probably only useful internally to this module; see deduplex().

    Args:
        dim_short (string): The short name of the dimension. E.g., "age"
        dataset (xarray Dataset): The Dataset we expect to include the dimension

    Raises:
        NameError: Dimension not found on Dataset

    Returns:
        string: The long name of the dimension. E.g., "fates_levage"
    """

    dim = "fates_lev" + dim_short
    if dim not in dataset.dims:
        raise NameError(f"Dimension {dim} not present in Dataset with dims {dataset.dims}")
    return dim


def _get_dim_combined(dim1_short, dim2_short):
    """Get duplexed dimension name, given two short names

    Args:
        dim1_short (string): Short name of first duplexed dimension. E.g., when de-duplexing
                             fates_levscpf, dim1_short=scls.
        dim2_short (string): Short name of second duplexed dimension. E.g., when de-duplexing
                             fates_levscpf, dim2_short=pft.

    Returns:
        string: Duplexed dimension name
    """
    dim_combined = "fates_lev" + dim1_short + dim2_short

    # Handle further-shortened dim names
    if dim_combined == "fates_levcanleaf":
        dim_combined = "fates_levcnlf"
    elif dim_combined == "fates_levcanpft":
        dim_combined = "fates_levcapf"
    elif dim_combined == "fates_levcdamscls":
        dim_combined = "fates_levcdsc"
    elif dim_combined == "fates_levsclsage":
        dim_combined = "fates_levscag"
    elif dim_combined == "fates_levsclspft":
        dim_combined = "fates_levscpf"
    elif dim_combined == "fates_levlandusepft":
        dim_combined = "fates_levlupft"
    elif dim_combined == "fates_levcanleaf":
        dim_combined = "fates_levcnlf"

    return dim_combined


def deduplex(dataset, this_var, dim1_short, dim2_short, preserve_order=True):
    """Reshape a duplexed FATES dimension into its constituent dimensions

    For example, given a variable with dimensions
        (time, fates_levagepft, lat, lon),
    this will return a DataArray with dimensions
        (time, fates_levage, fates_levpft, lat, lon)
    Or with reorder=False:
        (time, fates_levpft, lat, lon, fates_levage).

    Args:
        dataset (xarray Dataset): Dataset containing the variable with dimension to de-duplex
        this_var (string or xarray DataArray): (Name of) variable with dimension to de-duplex
        dim1_short (string): Short name of first duplexed dimension. E.g., when de-duplexing
                             fates_levagepft, dim1_short=age.
        dim2_short (string): Short name of second duplexed dimension. E.g., when de-duplexing
                             fates_levagepft, dim2_short=pft.
        preserve_order (bool, optional): Preserve order of dimensions of input DataArray? Defaults
                                         to True. Might be faster if False. See examples above.

    Raises:
        RuntimeError: dim1_short == dim2_short (not yet handled)
        TypeError: Incorrect type of this_var
        NameError: Dimension not found on Dataset

    Returns:
        xarray DataArray: De-duplexed variable
    """

    if dim1_short == dim2_short:
        raise RuntimeError("deduplex() can't currently handle dim1_short==dim2_short")

    # Get DataArray
    if isinstance(this_var, xr.DataArray):
        da_in = this_var
    elif isinstance(this_var, str):
        da_in = dataset[this_var]
    else:
        raise TypeError("this_var must be either string or DataArray, not " + type(this_var))

    # Get combined dim name
    dim_combined = _get_dim_combined(dim1_short, dim2_short)
    if dim_combined not in da_in.dims:
        raise NameError(f"Dimension {dim_combined} not present in DataArray with dims {da_in.dims}")

    # Get individual dim names
    dim1 = _get_check_dim(dim1_short, dataset)
    dim2 = _get_check_dim(dim2_short, dataset)

    # Split multiplexed dimension into its components
    n_dim1 = len(dataset[dim1])
    da_out = (
        da_in.rolling({dim_combined: n_dim1}, center=False)
        .construct(dim1)
        .isel({dim_combined: slice(n_dim1 - 1, None, n_dim1)})
        .rename({dim_combined: dim2})
        .assign_coords({dim1: dataset[dim1]})
        .assign_coords({dim2: dataset[dim2]})
    )

    # Reorder so that the split dimensions are together and in the expected order
    if preserve_order:
        new_dim_order = []
        for dim in da_out.dims:
            if dim == dim2:
                new_dim_order.append(dim1)
            if dim != dim1:
                new_dim_order.append(dim)
        da_out = da_out.transpose(*new_dim_order)

    return da_out


def agefuel_to_age_by_fuel(agefuel_var, dataset):
    """function to reshape a fates multiplexed age and fuel size indexed variable to one indexed by age and fuel size
    first argument should be an xarray DataArray that has the FATES AGEFUEL dimension
    second argument should be an xarray Dataset that has the FATES FUEL dimension
    (possibly the dataset encompassing the dataarray being transformed)
    returns an Xarray DataArray with the size and pft dimensions disentangled"""
    return deduplex(dataset, agefuel_var, "age", "fuel", preserve_order=False)


def scpf_to_scls_by_pft(scpf_var, dataset):
    """function to reshape a fates multiplexed size and pft-indexed variable to one indexed by size class and pft
    first argument should be an xarray DataArray that has the FATES SCPF dimension
    second argument should be an xarray Dataset that has the FATES SCLS dimension
    (possibly the dataset encompassing the dataarray being transformed)
    returns an Xarray DataArray with the size and pft dimensions disentangled"""
    return deduplex(dataset, scpf_var, "scls", "pft", preserve_order=False)


def scag_to_scls_by_age(scag_var, dataset):
    """function to reshape a fates multiplexed size and pft-indexed variable to one indexed by size class and pft
     first argument should be an xarray DataArray that has the FATES SCAG dimension
     second argument should be an xarray Dataset that has the FATES age dimension
    (possibly the dataset encompassing the dataarray being transformed)                                                                                                                                                     returns an Xarray DataArray with the size and age dimensions disentangled"""
    return deduplex(dataset, scag_var, "scls", "age", preserve_order=False)

def lupft_to_landuse_by_pft(lupft_var, dataset):
    """function to reshape a fates multiplexed land use and pft-indexed variable to one indexed by land use class and pft
     first argument should be an xarray DataArray that has the FATES LUPFT dimension
     second argument should be an xarray Dataset that has the FATES pft dimension
    (possibly the dataset encompassing the dataarray being transformed)                                                                                                                                                     returns an Xarray DataArray with the land use and pft dimensions disentangled"""
    return deduplex(dataset, lupft_var, "landuse", "pft", preserve_order=False)

def cnlf_to_canopy_by_leaf(cnlf_var, dataset):
    """function to reshape a fates multiplexed canopy and leaf-indexed variable to one indexed by canopy and leaf layer
     first argument should be an xarray DataArray that has the FATES LUPFT dimension
     second argument should be an xarray Dataset that has the FATES pft dimension
    (possibly the dataset encompassing the dataarray being transformed)                                                                                                                                                     returns an Xarray DataArray with the land use and pft dimensions disentangled"""
    return deduplex(dataset, cnlf_var, "can", "leaf", preserve_order=False)


def monthly_to_annual(array):
    """calculate annual mena from monthly data, using unequal month lengths fros noleap calendar.
    originally written by Keith Lindsay."""
    mon_day = xr.DataArray(
        np.array([31.0, 28.0, 31.0, 30.0, 31.0, 30.0, 31.0, 31.0, 30.0, 31.0, 30.0, 31.0]),
        dims=["month"],
    )
    mon_wgt = mon_day / mon_day.sum()
    return (
        array.rolling(time=12, center=False)  # rolling
        .construct("month")  # construct the array
        .isel(
            time=slice(11, None, 12)
        )  # slice so that the first element is [1..12], second is [13..24]
        .dot(mon_wgt, dims=["month"])
    )


def monthly_to_month_by_year(array):
    """go from monthly data to month x year data (for calculating climatologies, etc"""
    return (
        array.rolling(time=12, center=False)  # rolling
        .construct("month")  # construct the array
        .isel(time=slice(11, None, 12))
        .rename({"time": "year"})
    )
