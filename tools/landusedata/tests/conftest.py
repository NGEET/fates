import os
import pytest
import xarray as xr

@pytest.fixture(scope="module")
def target_file_location():
    return 'tests/resources/surfdata_4x5.nc'

@pytest.fixture(scope="module")
def static_file_location():
    return 'tests/resources/staticData_quarterdeg.nc'

@pytest.fixture(scope="module")
def landusepft_file_location():
    return 'tests/resources/CLM5_current_luhforest_deg025.nc'

@pytest.fixture(scope="module")
def lupftsurf_file_location():
    return 'tests/resources/CLM5_current_surf_deg025.nc'

@pytest.fixture(scope="function")
def static_dataset(static_file_location):
    dataset = xr.open_dataset(static_file_location)
    dataset = dataset.astype('float64')
    return dataset

@pytest.fixture(scope="function")
def target_dataset(target_file_location):
    dataset = xr.open_dataset(target_file_location)
    dataset = dataset.astype('float64')
    return dataset

@pytest.fixture(scope="function")
def landusepft_dataset(landusepft_file_location):
    dataset = xr.open_dataset(landusepft_file_location)
    dataset = dataset.astype('float64',casting='safe')
    return dataset

@pytest.fixture(scope="function")
def mock_mask(static_file_location):
    dataset = xr.open_dataset(static_file_location)
    # This is duplicative of the DefineMask function.  Should it simply be called?
    mask = (1.-dataset.icwtr) / (1.-dataset.icwtr)
    return mask
