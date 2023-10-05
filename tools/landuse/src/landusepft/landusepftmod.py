import xarray as xr

# Open the LUH2 static data file
def ImportStaticLUH2File(filename):
    dataset = xr.open_dataset(filename)

    # Check to see if the imported dataset has correct variables
    listcheck = ['ptbio', 'fstnf', 'carea', 'icwtr', 'ccode', 'lat_bounds', 'lon_bounds']
    if list(dataset.var()) != listcheck:
        raise TypeError("incorrect file, must be LUH2 static file")

    # Convert all data from single to double precision
    dataset = dataset.astype('float64')
    return dataset

# Open the CLM5 landuse x pft data file
def ImportLandusePFTFile(filename):
    dataset = xr.open_dataset(filename)

    # Check to see if the imported dataset has the the percent
    # natural pft variable present.
    if 'PCT_NAT_PFT' not in list(dataset.var()):
        raise TypeError("incorrect file, must be CLM5 landuse file")

    # change the percent natural pft from single to double precision
    # for downstream calculations
    dataset['PCT_NAT_PFT'] = dataset.PCT_NAT_PFT.astype('float64')
    return dataset

# Add lat/lon coordinates to the CLM5 landuse dataset
# While the lat and lon are available as variables, they are
# not defined as 'coords' in the imported dataset
def AddLatLonCoordinates(dataset):
    dataset['lon'] = dataset.lon * 0.25 - 180. + 1./8.
    dataset['lat'] = dataset.lat * 0.25 - 90. + 1./8.
    return dataset

# Define the land/ocean mask based on the ice/water data
# from the LUH2 static data set
def DefineMask(dataset):
    try:
        mask = (1.-dataset.icwtr) / (1.-dataset.icwtr)
    except AttributeError:
        raise AttributeError("incorrect dataset, must be static luh2 dataset")
    else:
        return mask

# Renormalize the pft percentages without the bareground
def RenormalizePFTs(dataset, mask):
    # Remove the bareground pft index from the dataset
    percent = dataset.PCT_NAT_PFT.isel(natpft=slice(1,None))
    # Normalize and apply mask
    percent = percent / percent.sum(dim='natpft') * mask
    return percent

# Steps
# - concatenate all this information together (including mask)
# - add/adjust lat/lon names and time input for regridding
# - regrid using luh2mod.py tooling

