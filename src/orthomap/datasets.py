#!/usr/bin/python
# -*- coding: UTF-8 -*-


"""
Author: Kristian K Ullrich
date: April 2023
email: ullrich@evolbio.mpg.de
License: GPL-3
"""


import os
import wget
import scanpy as sc


def ensembl105(datapath='.'):
    """
    OrthoFinder results for all translated coding sequences (CDS) from Ensembl release-105
    (keeping only longest isoforms)
    and Xtropicalisv9.0.Named.primaryTrs.pep.fa from www.xenbase.org.

    All files can be obtained from here:
    https://doi.org/10.5281/zenodo.7242263

    :param datapath: Path to safe dataset.
    :return: Path to Orthogroups.GeneCount file, OrthoGroups file and species list file.

    :type datapath: str
    :rtype: list of str

    Example
    -------
    >>> from orthomap import datasets
    >>> datasets.ensembl105(datapath='.')
    """
    if not os.path.exists(datapath):
        print('datapath does not exist, is created now')
        os.makedirs(name=datapath)
    oc_filename = os.path.join(datapath,
                               'ensembl_105_orthofinder_Orthogroups.GeneCount.tsv.zip')
    og_filename = os.path.join(datapath,
                               'ensembl_105_orthofinder_Orthogroups.tsv.zip')
    sl_filename = os.path.join(datapath,
                               'ensembl_105_orthofinder_species_list.tsv')
    oc_url = 'https://zenodo.org/record/7796253/files/ensembl_105_orthofinder_Orthogroups.GeneCount.tsv.zip'
    og_url = 'https://zenodo.org/record/7796253/files/ensembl_105_orthofinder_Orthogroups.tsv.zip'
    sl_url = 'https://zenodo.org/record/7796253/files/ensembl_105_orthofinder_species_list.tsv'

    wget.download(url=oc_url,
                  out=datapath)
    wget.download(url=og_url,
                  out=datapath)
    wget.download(url=sl_url,
                  out=datapath)
    return [oc_filename,
            og_filename,
            sl_filename]


def zebrafish_gtf(datapath='.'):
    """
    Download GTF for species Danio rerio from ensembl release-105
    https://ftp.ensembl.org/pub/release-105/gtf/danio_rerio/

    :param datapath: Path to safe dataset.
    :return: Path to GTF file.

    :type datapath: str
    :rtype: str

    Example
    -------
    >>> from orthomap import datasets
    >>> datasets.zebrafish_gtf(datapath='.')
    """
    if not os.path.exists(datapath):
        print('datapath does not exist, is created now')
        os.makedirs(name=datapath)
    zebrafish_gtf_filename = os.path.join(datapath,
                                          'Danio_rerio.GRCz11.105.gtf.gz')
    zebrafish_gtf_url = 'https://ftp.ensembl.org/pub/release-105/gtf/danio_rerio/Danio_rerio.GRCz11.105.gtf.gz'
    wget.download(url=zebrafish_gtf_url,
                  out=datapath)
    return zebrafish_gtf_filename


def mouse_gtf(datapath='.'):
    """
    Download GTF for species Mus musculus from ensembl release-105
    https://ftp.ensembl.org/pub/release-105/gtf/mus_musculus/

    :param datapath: Path to safe dataset.
    :return: Path to GTF file.

    :type datapath: str
    :rtype: str

    Example
    -------
    >>> from orthomap import datasets
    >>> datasets.mouse_gtf(datapath='.')
    """
    if not os.path.exists(datapath):
        print('datapath does not exist, is created now')
        os.makedirs(name=datapath)
    mouse_gtf_filename = os.path.join(datapath,
                                      'Mus_musculus.GRCm39.105.gtf.gz')
    mouse_gtf_url = 'https://ftp.ensembl.org/pub/release-105/gtf/mus_musculus/Mus_musculus.GRCm39.105.gtf.gz'
    wget.download(url=mouse_gtf_url,
                  out=datapath)
    return mouse_gtf_filename


def sun21_orthomap(datapath='.'):
    """
    Pre-calculated orthomap for Caenorhabditis elegans from:

    Sun, S., Rödelsperger, C. and Sommer, R.J., 2021.
    Single worm transcriptomics identifies a developmental core network of oscillating genes with deep
    conservation across nematodes.
    Genome research, 31(9), pp.1590-1601.

    :param datapath: Path to safe dataset.
    :return: Path to Orthomap file.

    :type datapath: str
    :rtype: str

    Example
    -------
    >>> from orthomap import datasets
    >>> datasets.sun21_orthomap(datapath='.')
    """
    if not os.path.exists(datapath):
        print('datapath does not exist, is created now')
        os.makedirs(name=datapath)
    sun21_orthomap_filename = os.path.join(datapath,
                                           'Sun2021_Orthomap.tsv')
    sun21_orthomap_url = 'https://zenodo.org/record/7783163/files/Sun2021_Orthomap.tsv'
    wget.download(url=sun21_orthomap_url,
                  out=datapath)
    return sun21_orthomap_filename


def zebrafish_orthomap(datapath='.'):
    """
    Pre-calculated and gene ID matched orthomap for Danio rerio extracted from OrthoFinder results:

    OrthoFinder results for all translated coding sequences (CDS) from Ensembl release-105
    (keeping only longest isoforms)
    and Xtropicalisv9.0.Named.primaryTrs.pep.fa from www.xenbase.org.

    All files can be obtained from here:
    https://doi.org/10.5281/zenodo.7242263

    Gene ID matching was done using the following GTF file for species Danio rerio from ensembl release-105:
    https://ftp.ensembl.org/pub/release-105/gtf/danio_rerio/

    :param datapath: Path to safe dataset.
    :return: Path to Orthomap file.

    :type datapath: str
    :rtype: str

    Example
    -------
    >>> from orthomap import datasets
    >>> datasets.zebrafish_orthomap(datapath='.')
    """
    if not os.path.exists(datapath):
        print('datapath does not exist, is created now')
        os.makedirs(name=datapath)
    zebrafish_orthomap_filename = os.path.join(datapath,
                                               'zebrafish_ensembl_105_orthomap.tsv')
    zebrafish_orthomap_url = 'https://zenodo.org/record/8345243/files/zebrafish_ensembl_105_orthomap.tsv'
    wget.download(url=zebrafish_orthomap_url,
                  out=datapath)
    return zebrafish_orthomap_filename


def ma21_fst(datapath='.'):
    """
    Pre-calculated TajimaD, NormalizedPi, FayWu and Fst for Caenorhabditis elegans from:

    Ma, F., Lau, C.Y. and Zheng, C., 2021.
    Large genetic diversity and strong positive selection in F-box and GPCR genes among the wild isolates of
    Caenorhabditis elegans.
    Genome Biology and Evolution, 13(5), p.evab048.

    :param datapath: Path to safe dataset.
    :return: Path to diversity file.

    :type datapath: str
    :rtype: str

    Example
    -------
    >>> from orthomap import datasets
    >>> datasets.ma21_fst(datapath='.')
    """
    if not os.path.exists(datapath):
        print('datapath does not exist, is created now')
        os.makedirs(name=datapath)
    ma21_orthomap_filename = os.path.join(datapath,
                                          'Ma2021_Fst.tsv')
    ma21_orthomap_url = 'https://zenodo.org/record/7802177/files/Ma2021_Fst.tsv'
    wget.download(url=ma21_orthomap_url,
                  out=datapath)
    return ma21_orthomap_filename


def cazet22_orthomap(datapath='.'):
    """
    Pre-calculated orthomap for Hydra vulgaris from:

    Cazet, Jack, Stefan Siebert, Hannah Morris Little, Philip Bertemes, Abby S. Primack, Peter Ladurner,
    Matthias Achrainer et al. (2022)
    New Hydra genomes reveal conserved principles of hydrozoan transcriptional regulation.,
    bioRxiv, 2022.06.21.496857.

    :param datapath: Path to safe dataset.
    :return: Path to Orthomap file.

    :type datapath: str
    :rtype: str

    Example
    -------
    >>> from orthomap import datasets
    >>> datasets.cazet22_orthomap(datapath='.')
    """
    if not os.path.exists(datapath):
        print('datapath does not exist, is created now')
        os.makedirs(name=datapath)
    cazet22_orthomap_filename = os.path.join(datapath,
                                             'Cazet2022_Orthomap.tsv')
    cazet22_orthomap_url = 'https://zenodo.org/record/7404798/files/Cazet2022_Orthomap.tsv'
    wget.download(url=cazet22_orthomap_url,
                  out=datapath)
    return cazet22_orthomap_filename


def packer19(datapath='.'):
    """
    scRNA count data for Caenorhabditis elegans from:

    Packer, J.S., Zhu, Q., Huynh, C., Sivaramakrishnan, P., Preston, E., Dueck, H., Stefanik, D.,
    Tan, K., Trapnell, C., Kim, J. and Waterston, R.H., 2019.
    A lineage-resolved molecular atlas of C. elegans embryogenesis at single-cell resolution.
    Science, 365(6459), p.eaax1971.

    All files can be obtained from here:
    https://doi.org/10.5281/zenodo.7245547

    :param datapath: Path to safe dataset.
    :return: AnnData object.

    :type datapath: str
    :rtype: AnnData

    Example
    -------
    >>> from orthomap import datasets
    >>> datasets.packer19(datapath='.')
    """
    if not os.path.exists(datapath):
        print('datapath does not exist, is created now')
        os.makedirs(name=datapath)
    packer19_filename = os.path.join(datapath,
                                     'GSE126954.h5ad')
    packer19_url = 'https://zenodo.org/record/7496490/files/GSE126954.h5ad'
    wget.download(url=packer19_url,
                  out=datapath)
    adata = sc.read(packer19_filename)
    return adata


def packer19_small(datapath='.'):
    """
    scRNA count data for Caenorhabditis elegans from:

    Packer, J.S., Zhu, Q., Huynh, C., Sivaramakrishnan, P., Preston, E., Dueck, H., Stefanik, D.,
    Tan, K., Trapnell, C., Kim, J. and Waterston, R.H., 2019.
    A lineage-resolved molecular atlas of C. elegans embryogenesis at single-cell resolution.
    Science, 365(6459), p.eaax1971.

    All files can be obtained from here:
    https://doi.org/10.5281/zenodo.7245547

    :param datapath: Path to safe dataset.
    :return: AnnData object.

    :type datapath: str
    :rtype: AnnData

    Example
    -------
    >>> from orthomap import datasets
    >>> datasets.packer19_small(datapath='.')
    """
    if not os.path.exists(datapath):
        print('datapath does not exist, is created now')
        os.makedirs(name=datapath)
    packer19_small_filename = os.path.join(datapath,
                                           'packer19_small.h5ad')
    packer19_small_url = 'https://zenodo.org/record/7496490/files/packer19_small.h5ad'
    wget.download(url=packer19_small_url,
                  out=datapath)
    adata = sc.read(packer19_small_filename)
    return adata


def cazet22(datapath='.'):
    """
    scRNA count data for Hydra vulgaris from:

    Cazet, Jack, Stefan Siebert, Hannah Morris Little, Philip Bertemes, Abby S. Primack, Peter Ladurner,
    Matthias Achrainer et al. (2022)
    New Hydra genomes reveal conserved principles of hydrozoan transcriptional regulation.,
    bioRxiv, 2022.06.21.496857.

    All files can be obtained from here:
    https://doi.org/10.5281/zenodo.7366178

    :param datapath: Path to safe dataset.
    :return: AnnData object.

    :type datapath: str
    :rtype: AnnData

    Example
    -------
    >>> from orthomap import datasets
    >>> datasets.cazet22(datapath='.')
    """
    if not os.path.exists(datapath):
        print('datapath does not exist, is created now')
        os.makedirs(name=datapath)
    cazet22_filename = os.path.join(datapath,
                                    'aepAtlasNonDub.h5ad')
    cazet22_url = 'https://zenodo.org/record/7369647/files/aepAtlasNonDub.h5ad'
    wget.download(url=cazet22_url,
                  out=datapath)
    adata = sc.read(cazet22_filename)
    return adata


def qiu22_zebrafish(datapath='.'):
    """
    combined scRNA count data for Danio rerio from:

    Qiu, C., Cao, J., Martin, B.K., Li, T., Welsh, I.C., Srivatsan, S., Huang, X., Calderon,
    D., Noble, W.S., Disteche, C.M. and Murray, S.A., 2022.
    Systematic reconstruction of cellular trajectories across mouse embryogenesis.
    Nature genetics, 54(3), pp.328-341.

    original scRNA count data from:

    Farrell, J.A., Wang, Y., Riesenfeld, S.J., Shekhar, K., Regev, A. and Schier, A.F., 2018.
    Single-cell reconstruction of developmental trajectories during zebrafish embryogenesis.
    Science, 360(6392), p.eaar3131.

    Wagner, D.E., Weinreb, C., Collins, Z.M., Briggs, J.A., Megason, S.G. and Klein, A.M., 2018.
    Single-cell mapping of gene expression landscapes and lineage in the zebrafish embryo.
    Science, 360(6392), pp.981-987.

    All files can be obtained from here:
    https://doi.org/10.5281/zenodo.7243602

    :param datapath: Path to safe dataset.
    :return: AnnData object.

    :type datapath: str
    :rtype: AnnData

    Example
    -------
    >>> from orthomap import datasets
    >>> datasets.qiu22_zebrafish(datapath='.')
    """
    if not os.path.exists(datapath):
        print('datapath does not exist, is created now')
        os.makedirs(name=datapath)
    qiu22_zebrafish_filename = os.path.join(datapath,
                                            'zebrafish_data.h5ad')
    qiu22_zebrafish_url = 'https://zenodo.org/record/7243603/files/zebrafish_data.h5ad'
    wget.download(url=qiu22_zebrafish_url,
                  out=datapath)
    adata = sc.read(qiu22_zebrafish_filename)
    return adata


def qiu22_frog(datapath='.'):
    """
    combined scRNA count data for Xenopus tropicalis from:

    Qiu, C., Cao, J., Martin, B.K., Li, T., Welsh, I.C., Srivatsan, S., Huang, X., Calderon,
    D., Noble, W.S., Disteche, C.M. and Murray, S.A., 2022.
    Systematic reconstruction of cellular trajectories across mouse embryogenesis.
    Nature genetics, 54(3), pp.328-341.

    original scRNA count data from:

    Briggs, J.A., Weinreb, C., Wagner, D.E., Megason, S., Peshkin, L., Kirschner, M.W. and Klein, A.M., 2018.
    The dynamics of gene expression in vertebrate embryogenesis at single-cell resolution.
    Science, 360(6392), p.eaar5780.

    All files can be obtained from here:
    https://doi.org/10.5281/zenodo.7244440

    :param datapath: Path to safe dataset.
    :return: AnnData object.

    :type datapath: str
    :rtype: AnnData

    Example
    -------
    >>> from orthomap import datasets
    >>> datasets.qiu22_frog(datapath='.')
    """
    if not os.path.exists(datapath):
        print('datapath does not exist, is created now')
        os.makedirs(name=datapath)
    qiu22_frog_filename = os.path.join(datapath,
                                       'frog_data.h5ad')
    qiu22_frog_url = 'https://zenodo.org/record/7244441/files/frog_data.h5ad'
    wget.download(url=qiu22_frog_url,
                  out=datapath)
    adata = sc.read(qiu22_frog_filename)
    return adata


def qiu22_mouse(datapath='.'):
    """
    combined scRNA count data for Mus musculus from:

    Qiu, C., Cao, J., Martin, B.K., Li, T., Welsh, I.C., Srivatsan, S., Huang, X., Calderon,
    D., Noble, W.S., Disteche, C.M. and Murray, S.A., 2022.
    Systematic reconstruction of cellular trajectories across mouse embryogenesis.
    Nature genetics, 54(3), pp.328-341.

    original scRNA count data from:

    Mohammed, H., Hernando-Herraez, I., Savino, A., Scialdone, A., Macaulay, I., Mulas, C., Chandra, T.,
    Voet, T., Dean, W., Nichols, J. and Marioni, J.C., 2017.
    Single-cell landscape of transcriptional heterogeneity and cell fate decisions during mouse early gastrulation.
    Cell reports, 20(5), pp.1215-1228.

    Cheng, S., Pei, Y., He, L., Peng, G., Reinius, B., Tam, P.P., Jing, N. and Deng, Q., 2019.
    Single-cell RNA-seq reveals cellular heterogeneity of pluripotency transition and X chromosome dynamics
    during early mouse development.
    Cell reports, 26(10), pp.2593-2607.

    Pijuan-Sala, B., Griffiths, J.A., Guibentif, C., Hiscock, T.W., Jawaid, W., Calero-Nieto, F.J.,
    Mulas, C., Ibarra-Soria, X., Tyser, R.C., Ho, D.L.L. and Reik, W., 2019.
    A single-cell molecular map of mouse gastrulation and early organogenesis.
    Nature, 566(7745), pp.490-495.

    Cao, J., Spielmann, M., Qiu, X., Huang, X., Ibrahim, D.M., Hill, A.J., Zhang, F., Mundlos,
    S., Christiansen, L., Steemers, F.J. and Trapnell, C., 2019.
    The single-cell transcriptional landscape of mammalian organogenesis.
    Nature, 566(7745), pp.496-502.

    All files can be obtained from here:
    https://doi.org/10.5281/zenodo.7244567

    :param datapath: Path to safe dataset.
    :return: AnnData object.

    :type datapath: str
    :rtype: AnnData

    Example
    -------
    >>> from orthomap import datasets
    >>> datasets.qiu22_mouse(datapath='.')
    """
    if not os.path.exists(datapath):
        print('datapath does not exist, is created now')
        os.makedirs(name=datapath)
    qiu22_mouse_filename = os.path.join(datapath,
                                        'mouse_data.h5ad')
    qiu22_mouse_url = 'https://zenodo.org/record/7244568/files/mouse_data.h5ad'
    wget.download(url=qiu22_mouse_url,
                  out=datapath)
    adata = sc.read(qiu22_mouse_filename)
    return adata


def mytai_example(datapath='.'):
    """
    expression count data for Arabidopsis thaliana from:

    Drost, H.G., Janitza, P., Grosse, I. and Quint, M., 2017.
    Cross-kingdom comparison of the developmental hourglass.
    Current Opinion in Genetics & Development, 45, pp.69-75.

    All files can be obtained from here:
    https://doi.org/10.5281/zenodo.7242263

    :param datapath: Path to safe dataset.
    :return: AnnData object.

    :type datapath: str
    :rtype: AnnData

    Example
    -------
    >>> from orthomap import datasets
    >>> datasets.mytai_example(datapath='.')
    """
    if not os.path.exists(datapath):
        print('datapath does not exist, is created now')
        os.makedirs(name=datapath)
    mytai_example_filename = os.path.join(datapath,
                                          'PhyloExpressionSetExample.h5ad')
    mytai_example_url = 'https://zenodo.org/record/7783163/files/PhyloExpressionSetExample.h5ad'
    wget.download(url=mytai_example_url,
                  out=datapath)
    adata = sc.read(mytai_example_filename)
    return adata
