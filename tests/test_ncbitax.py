from orthomap import ncbitax
import argparse
from os.path import exists


def test_define_parse():
    parse = ncbitax.define_parser()
    assert isinstance(parse, argparse.ArgumentParser)
