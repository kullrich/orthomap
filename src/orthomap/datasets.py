#!/usr/bin/python
# -*- coding: UTF-8 -*-


"""
Author: Kristian K Ullrich
date: November 2022
email: ullrich@evolbio.mpg.de
License: GPL-3
"""


import os
import wget


def ensembl105(datapath='.'):
    """
    OrthoFinder results for all translated coding sequences (CDS) from Ensembl release-105
    (keeping only longest isoforms)
    and Xtropicalisv9.0.Named.primaryTrs.pep.fa from www.xenbase.org.

    All files can be obtained from here:
    https://doi.org/10.5281/zenodo.7242263

    :param datapath:
    :return:
    """
    oc_filename = os.path.join(datapath, 'ensembl_105_orthofinder_Orthogroups.GeneCount.tsv.zip')
    og_filename = os.path.join(datapath, 'ensembl_105_orthofinder_Orthogroups.zip.tsv')
    sl_filename = os.path.join(datapath, 'ensembl_105_orthofinder_species_list.tsv')
    oc_url = 'https://github.com/kullrich/orthomap/raw/main/examples/' \
             'ensembl_105_orthofinder_Orthogroups.GeneCount.tsv.zip'
    og_url = 'https://github.com/kullrich/orthomap/raw/main/examples/' \
             'ensembl_105_orthofinder_Orthogroups.tsv.zip'
    sl_url = 'https://github.com/kullrich/orthomap/raw/main/examples/' \
             'ensembl_105_orthofinder_species_list.tsv'
    wget.download(oc_url)
    wget.download(og_url)
    wget.download(sl_url)
    return [oc_filename, og_filename, sl_filename]
