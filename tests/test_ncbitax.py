from orthomap import ncbitax
import argparse
import os


def test_define_parse():
    parse = ncbitax.define_parser()
    assert isinstance(parse, argparse.ArgumentParser)

def test_update_ncbi():
    path = os.path.expanduser('~/.etetoolkit/taxa.sqlite')
    path_exist = os.path.exists(path)

    if path_exist == False:
        ncbitax.update_ncbi()
        assert os.path.exists(path)