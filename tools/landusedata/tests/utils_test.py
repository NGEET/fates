import pytest

from landusedata import utils

def test_RenameLatLon(target_dataset):
    target_dataset = utils.RenameLatLon(target_dataset)
    dim_list = list(target_dataset.dims)
    latlon_list = ['lat','lon']
    test_list = [dim for latlon in latlon_list for dim in dim_list
                 if latlon in dim]
    assert latlon_list == test_list
