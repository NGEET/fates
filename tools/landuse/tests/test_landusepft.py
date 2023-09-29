# import pytest
from landusepft import landusepft

def test_TempFunc():
    assert landusepft.basicfunc() == "yup"

# Possible tests
# - Check that import fails if argument is not given the appropriate landuse pft data
# - Check that import fails if argument is not given the appropriate static file data
# - Check by assertion that mask value is zero over specific points on a 1/4 deg grid?
# - Check by assertion that sum of certain land use and/or pft combinations sum to 100?
