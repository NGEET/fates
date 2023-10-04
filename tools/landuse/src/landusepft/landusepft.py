import xarray as xr

# Open the LUH2 static data file
def ImportStaticLUH2File(filename):
    dataset = xr.open_dataset(filename)
    listcheck = ['ptbio', 'fstnf', 'carea', 'icwtr', 'ccode', 'lat_bounds', 'lon_bounds']
    if list(dataset.var()) != listcheck:
        raise TypeError("incorrect file, must be LUH2 static file")
    dataset = dataset.astype('float64')
    return dataset

# Open the CLM5 landuse x pft data file
def ImportLandusePFTFile(filename):
    dataset = xr.open_dataset(filename)
    if 'PCT_NAT_PFT' not in list(dataset.var()):
        raise TypeError("incorrect file, must be CLM5 landuse file")
    dataset = dataset.astype('float64')
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

def RenormalizePFTs(dataset):
    percent = dataset.PCT_NAT_PFT.isel(natpft=slice(1,None))
    # percent = percent / percent * 100.0 * mask
    percent = percent / percent.sum(dim='natpft') * 100.0
    return percent

# Steps
# - import clm landuse-pft data (1/4 degree)
# - import luh2 static data file (1/4 degree)
# - set the mask based on static data file `icwtr` data
# - calculate the percentage of forest, pasture, and other (i.e. range) pft percentages
# - calculate the primary and secondard forest percent using the static data `fstnf`
# primary_secondary_percent = luh2_staticdata.fstnf * forest_pft_percent + (1.- luh2_staticdata.fstnf) * other_pft_percent
# - concatenate all this information together (including mask)
# - add/adjust lat/lon names and time input for regridding
# - regrid using luh2mod.py tooling

