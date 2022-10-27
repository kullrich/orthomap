from orthomap import qlin
import argparse
from ete3 import NCBITaxa


def test_define_parse():
    parse = qlin.define_parser()
    assert isinstance(parse, argparse.ArgumentParser)


def test_get_qlin_q_argument():
    q = "Danio rerio"
    info = qlin.get_qlin(q=q)

    assert isinstance(info, list)
    assert info[0] == q
    assert isinstance(info[3], dict)
    assert isinstance(info[-1], str)

    for info_index in [2, 4, 5]:
        assert isinstance(info[info_index], list)


def test_get_qlin_qt_argument():
    qt = "7955"
    info = qlin.get_qlin(qt=qt)

    assert isinstance(info, list)
    assert info[1] == int(qt)
    assert isinstance(info[3], dict)
    assert isinstance(info[-1], str)

    for info_index in [2, 4, 5]:
        assert isinstance(info[info_index], list)


def test_get_qlin_q_and_qt_argument():
    q = "Danio rerio"
    qt = "7955"
    info = qlin.get_qlin(q=q, qt=qt)

    assert isinstance(info, list)
    assert info[0] == q
    assert info[1] == int(qt)


def test_get_qlin_q_with_wrong_qt_argument():
    """If the name and taxid do not much then `orthomap` returns information
    based on the taxid."""

    q = "Danio rerio"
    qt = "7956"
    info = qlin.get_qlin(q=q, qt=qt)

    assert isinstance(info, list)
    assert info[0] != q
    assert info[1] == int(qt)


def test_get_qtid():
    ncbi = NCBITaxa()
    q = "Carassius"
    qt = "7956"

    info = qlin.get_qtid(ncbi=ncbi, q=q, qt=qt)

    info2 = qlin.get_qlin(q=q, qt=qt)

    assert isinstance(info, list)
    assert info[0] == q
    assert info[1] == int(qt)
    assert info == info2
