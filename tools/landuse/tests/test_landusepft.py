import pytest

from landusepft import landusepft

# Postive test case for importing LUH2 static data file
def test_posImportStaticLUH2File(static_file_location):
    data = landusepft.ImportStaticLUH2File(static_file_location)
    static_variables = list(data.var())
    assert static_variables == ['ptbio', 'fstnf', 'carea', 'icwtr', 'ccode', 'lat_bounds', 'lon_bounds']

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
