import pytest
import numpy as np

from landusepft import landusepft

# Postive test case for importing LUH2 static data file
def test_posImportStaticLUH2File(static_file_location):
    data = landusepft.ImportStaticLUH2File(static_file_location)
    static_variables = list(data.var())
    assert static_variables == ['ptbio', 'fstnf', 'carea', 'icwtr', 'ccode', 'lat_bounds', 'lon_bounds']

# Negative test case for importing incorrect file via static luh2 file open function
def test_negImportStaticLUH2File(landusepft_file_location):
    with pytest.raises(TypeError) as exp:
        landusepft.ImportStaticLUH2File(landusepft_file_location)
    assert str(exp.value) == "incorrect LUH2 file, must be static file"

# Positive test case for importing landuse x pft data file
def test_posImportLandusePFTFile(landusepft_file_location):
    data = landusepft.ImportLandusePFTFile(landusepft_file_location)
    landusepft_variables = list(data.var())
    assert landusepft_variables == ['EDGEN',
                                   'EDGEE',
                                   'EDGES',
                                   'EDGEW',
                                   'LAT',
                                   'LATIXY',
                                   'LON',
                                   'LONGXY',
                                   'LANDMASK',
                                   'LANDFRAC',
                                   'AREA',
                                   'PCT_NAT_PFT']


# Define Mask function unit tests
# - Is there a better way to simply test that the mask matches all gridcells?
# - Should we test what happens when the wrong dataset is used (i.e. landuse x pft)?

# Make sure that the generated mask makes sense in that a known gridcell is nan
def test_posDefineMask(static_dataset):
    maskoutput = landusepft.DefineMask(static_dataset)
    assert np.isnan(maskoutput[0][0])

# Make sure a known gridcell that should be land is not nan
# TODO: change this to use sel indexing
def test_negDefineMask(static_dataset):
    maskoutput = landusepft.DefineMask(static_dataset)
    assert not(np.isnan(maskoutput[360][850]))

# Bare ground fraction removal unit test
# Make sure that the pft fractions sum to 100% on a known land gridcell
# TODO: this needs to take the mask as the input
# This won't work because there is no lat/lon coordinate for the PCT_NAT_PFT dimensions
# even though xarray recognizes that there technically are dimensions specified
def test_posRenormalizeNoBareGround(landusepft_dataset):
    percent = sum(landusepft_dataset.PCT_NAT_PFT.sel(lat=0.125, lon=20.125))
    assert percent == 100.0

# forest_pft_percent = fin_forestdata.PCT_NAT_PFT.isel(natpft=slice(1,None)) / fin_forestdata.PCT_NAT_PFT.isel(natpft=slice(1,None)).sum(dim='natpft') * 100. * landoceanmask

# Negative test case for importing the LUH2 static data file
# def test_negImportStaticLUH2File(mockluh2_file_location):
#     with pytest.raises(ValueError) as exp:
#         data = landusepft.ImportStaticLUH2File(mockluh2_file_location)
#         variables = list(data.var())
#     assert str(exp.value) == "'Error: file imported is not the LUH2 static data file'"

# def test_negImportFile(static_file_location):
#     with pytest.raises(AttributeError) as exp:
#         data = landusepft.ImportStaticLUH2File(static_file_location)
#         static_variables = list(data.var())
#     assert str(exp.value) == "'list' object has no attribute 'var'"

# Possible tests
# - Check that import fails if argument is not given the appropriate landuse pft data
# - Check that import fails if argument is not given the appropriate static file data
# - Check by assertion that mask value is zero over specific points on a 1/4 deg grid?
# - Check by assertion that sum of certain land use and/or pft combinations sum to 100?
