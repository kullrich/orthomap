import argparse

from orthomap import eggnog2orthomap


def test_define_parse():
    parse = eggnog2orthomap.define_parser()
    assert isinstance(parse, argparse.ArgumentParser)