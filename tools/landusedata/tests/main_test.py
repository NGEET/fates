import pytest
from argparse import ArgumentError
from argparse import ArgumentParser as parser

from landusedata._main import main

# TODO: parameterize this
def test_main(capsys):
    """
    Postive case testing for choice input
    """
    main(['luh2'])
    out, err = capsys.readouterr()
    assert out == 'calling luh2 code\n'
    assert err == ''

# This test only catches the system exit. 
def test_main_neg(capsys):
    """
    Negative case testing for choice input

    This only tests the system exit. Testing the ArgumentError
    is a little tricky to catch as it is upstream of SystemExit.
    """
    with pytest.raises(SystemExit) as exp:
        main(['notallowed'])
    assert str(exp.value) == "2"
