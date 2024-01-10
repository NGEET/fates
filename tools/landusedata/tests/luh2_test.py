import pytest

from landusedata import luh2mod

# Postive test case for importing LUH2 static data file
def test_posImportLUH2StaticFile(static_file_location):
    data = luh2mod.ImportLUH2StaticFile(static_file_location)
    static_variables = list(data.var())
    assert static_variables == ['ptbio',
                                'fstnf',
                                'carea',
                                'icwtr',
                                'ccode',
                                'lat_bounds',
                                'lon_bounds']

# Negative test case for importing incorrect file via static luh2 file open function
def test_negImportLUH2StaticFile(landusepft_file_location):
    with pytest.raises(TypeError) as exp:
        luh2mod.ImportLUH2StaticFile(landusepft_file_location)
    assert str(exp.value) == "incorrect file, must be LUH2 static file"
