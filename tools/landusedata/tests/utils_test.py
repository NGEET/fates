import pytest

from landusedata import utils

def test_RegridTargetPrep_dims(target_dataset):
    target_dataset = utils._RegridTargetPrep(target_dataset)
    dim_list = list(target_dataset.dims)
    latlon_list = ['lat','lon']
    test_list = [dim for latlon in latlon_list for dim in dim_list
                 if latlon in dim]
    assert latlon_list == test_list

def test_RegridTargetPrep_coords(target_dataset):
    target_dataset = utils._RegridTargetPrep(target_dataset)
    coords_list = list(target_dataset.coords)
    latlon_list = ['lat','lon']
    test_list = [coord for latlon in latlon_list for coord in coords_list
                 if latlon in coord]
    assert latlon_list == test_list

def test_negImportRegridTarget(landusepft_file_location):
    with pytest.raises(TypeError) as exp:
        utils.ImportRegridTarget(landusepft_file_location)
    assert str(exp.value) == "incorrect file, must be surface dataset"
