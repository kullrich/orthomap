import orthomap


def test_title():
    assert orthomap.__title__ == "orthomap"


def test_version():
    assert orthomap.__version__ == "0.0.1"


def test_license():
    assert orthomap.__license__ == "GPL-3"
