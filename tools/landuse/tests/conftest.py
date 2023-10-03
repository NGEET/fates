import os
import pytest
import xarray as xr

@pytest.fixture(scope="module")
def static_file_location():
    return 'tests/resources/staticData_quarterdeg.nc'

@pytest.fixture(scope="module")
def landusepft_file_location():
    return 'tests/resources/CLM5_current_luhforest_deg025.nc'

@pytest.fixture(scope="module")
def lupftsurf_file_location():
    return 'tests/resources/CLM5_current_surf_deg025.nc'

@pytest.fixture(scope="module")
def static_dataset(static_file_location):
    dataset = xr.open_dataset(static_file_location)
    return dataset

@pytest.fixture(scope="module")
def landusepft_dataset(landusepft_file_location):
    dataset = xr.open_dataset(landusepft_file_location)
    return dataset

# @pytest.fixture(scope="module")
# def mock_mask(landusepft_file_location):
#     mask = xr.open_dataset(landusepft_file_location)
#     return mask
