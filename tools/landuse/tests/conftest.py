import os
import pytest

@pytest.fixture(scope="module")
def static_file_location():
    return 'tests/resources/staticData_quarterdeg.nc'
