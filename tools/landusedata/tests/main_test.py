from landusedata._main import main

def test_main(capsys):
    main(['--output','testmain'])
    out, err = capsys.readouterr()
    assert out == 'data type is testmain\n'
    assert err == ''
