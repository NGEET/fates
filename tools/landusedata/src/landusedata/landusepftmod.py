import xarray as xr

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
    dataset['lon'] = dataset.LON
    dataset['lat'] = dataset.LAT
    return dataset

# Renormalize the pft percentages without the bareground
def RenormalizePFTs(dataset):
    # Remove the bareground pft index from the dataset
    percent = dataset.PCT_NAT_PFT.isel(natpft=slice(1,None))
    # Normalize
    percent = percent / percent.sum(dim='natpft')
    return percent
