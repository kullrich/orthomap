#!/usr/bin/python
# -*- coding: UTF-8 -*-

import argparse
import pandas as pd
from orthomap import cds2aa


def test_define_parser():
    parse = cds2aa.define_parser()
    assert isinstance(parse, argparse.ArgumentParser)