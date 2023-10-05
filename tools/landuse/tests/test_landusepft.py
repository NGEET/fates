import pytest
import numpy as np

from landusepft import landusepft

# Postive test case for importing LUH2 static data file
def test_posImportStaticLUH2File(static_file_location):
    data = landusepft.ImportStaticLUH2File(static_file_location)
    static_variables = list(data.var())
    assert static_variables == ['ptbio',
                                'fstnf',
                                'carea',
                                'icwtr',
                                'ccode',
                                'lat_bounds',
                                'lon_bounds']

# Negative test case for importing incorrect file via static luh2 file open function
def test_negImportStaticLUH2File(landusepft_file_location):
    with pytest.raises(TypeError) as exp:
        landusepft.ImportStaticLUH2File(landusepft_file_location)
    assert str(exp.value) == "incorrect file, must be LUH2 static file"

# Positive test case for importing landuse x pft data file
# TODO: parameterize this call to run all the landusepft files
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

# Negative test case for importing incorrect file via CLM5 landuse file open function
def test_negImportLandusePFTFile(static_file_location):
    with pytest.raises(TypeError) as exp:
        landusepft.ImportLandusePFTFile(static_file_location)
    assert str(exp.value) == "incorrect file, must be CLM5 landuse file"

# Confirm that the data type has been converted upon import from single to double precision
def test_pos_ImportLandusePFTFileDtypeConversion(landusepft_file_location):
    dataset = landusepft.ImportLandusePFTFile(landusepft_file_location)
    assert np.dtype('float64') == dataset['PCT_NAT_PFT'].dtype

# Unit tests to update CLM5 landuse datasets with lat/lon coordinates
def test_latlon_posAddLatLonCoordinates(landusepft_dataset):
    dataset_modified = landusepft.AddLatLonCoordinates(landusepft_dataset)
    assert 'lat' in list(dataset_modified.coords) and 'lon' in list(dataset_modified.coords)

# Test that the math of the function makes sense.  We can sum that lat/lon variables
# to makes sure that the update is correct
def test_latlonsum_posAddLatLonCoordinates(landusepft_dataset):
    dataset_modified = landusepft.AddLatLonCoordinates(landusepft_dataset)
    assert (sum(abs(dataset_modified.lon)).values.item() == 129600.
           and sum(abs(dataset_modified.lat)).values.item() == 32400.)

# Define Mask function unit tests
# - Is there a better way to simply test that the mask matches all gridcells?
# - Should we test what happens when the wrong dataset is used (i.e. landuse x pft)?
# - Should this test that the mask is truly one or zero everywhere?

# Make sure that the generated mask makes sense in that a known gridcell is nan
def test_posDefineMask(static_dataset):
    maskoutput = landusepft.DefineMask(static_dataset)
    assert np.isnan(maskoutput[0][0])

# Make sure that the generated mask is equal to 1 where it isn't nan
def test_pos_ones_DefineMask(static_dataset):
    maskoutput = landusepft.DefineMask(static_dataset)
    assert (maskoutput.to_dataframe().dropna() == 1.).all().values.item()

# Make sure a known gridcell that should be land is not nan
# TODO: change this to use sel indexing
def test_neg_output_DefineMask(static_dataset):
    maskoutput = landusepft.DefineMask(static_dataset)
    assert not(np.isnan(maskoutput[360][850]))

def test_neg_input_DefineMask(landusepft_dataset):
    with pytest.raises(AttributeError) as exp:
        maskoutput = landusepft.DefineMask(landusepft_dataset)
    assert str(exp.value) == "incorrect dataset, must be static luh2 dataset"

# Bare ground fraction removal unit test
# Make sure that the pft fractions sum to 100% on a known land gridcell
def test_posRenormalizePFTs(landusepft_dataset,mock_mask):
    percent = landusepft.RenormalizePFTs(landusepft_dataset, mock_mask)

    # Sum along the natpft dimension only.  Use min_count to avoid
    # NaNs from masking turning into zeros
    percent = percent.sum(dim='natpft',min_count=1)

    # Convert to a dataframe to stack lat and lon in one dimension
    # and drop the NaNs.  Note that since percent is unassociated
    # with a dataset we need to pass the name argument when converting
    # to a dataframe.
    percent = percent.to_dataframe(name="testpercent").dropna(how='all')
    # breakpoint()

    # Check that all the summations are unity within a specific tolerance
    tolerance = 3.0e-16
    assert (abs(percent - 1.0) < tolerance).all().values.item()

# @pytest.mark.skip(reason="this needs more work")
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
