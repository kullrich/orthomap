# orthomap

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![pypi-badge](https://img.shields.io/pypi/v/orthomap)](https://pypi.org/project/orthomap)
[![docs-badge](https://readthedocs.org/projects/orthomap/badge/?version=latest)](https://orthomap.readthedocs.io/en/latest/?badge=latest)

## orthologous maps - evolutionary age index

[orthomap](https://github.com/kullrich/orthomap) is a python package to extract orthologous maps (in other words the evolutionary age of a given orthologous group) from OrthoFinder results.

## Installing `orthomap`

### Anaconda

The environment is created with `conda create` in which orthomap is installed.

If you do not have a working installation of Python 3.7 (or later), consider
installing [Miniconda] (see [Installing Miniconda]). Then run:

```shell
conda create --file environmen.yml
conda activate orthomap
```

Install `orthomap`:

```shell
pip install orthomap
```

### PyPI

Install `orthomap` into your current python environment:

```shell
pip install orthomap
```

## Documentation

Online documentation can be found [here](https://orthomap.readthedocs.io/en/latest/).

## Quick use

Update/download local ncbi taxonomic database:

``` 
from orthomap import ncbitax
ncbitax.update_ncbi()
```

Get query species lineage information:

example: Danio rerio

``` 
from orthomap import qlin
qlin.get_qlin(q = 'Danio rerio')
qlin.get_qlin(qt = '7955')
```

Extract orthomap from OrthoFinder result:

example: ensembl release-105

OrthoFinder results files can be found [here](https://doi.org/10.5281/zenodo.7242264)

```
from orthomap import orthomap
omap = orthomap.get_orthomap(qname = 'Danio_rerio.GRCz11.cds.longest',
    qt = '7955',
    sl = 'ensembl_105_orthofinder_species_list.tsv',
    oc = 'ensembl_105_orthofinder_Orthogroups.GeneCount.tsv',
    og = 'ensembl_105_orthofinder_Orthogroups.tsv')
```

Calculate transcriptome evolutionary index (TEI) for each cell of a scRNA data set:

example: Danio rerio - [http://tome.gs.washington.edu](http://tome.gs.washington.edu)
([Qui et al. 2022](https://www.nature.com/articles/s41588-022-01018-x))

`AnnData` file can be found [here](https://doi.org/10.5281/zenodo.7243602) 

```
from orthomap import orthomap2tei
```

## Development Version

To work with the latest version [on GitHub]: clone the repository and `cd` into its root directory.

```shell
git clone kullrich/orthomap
cd orthomap
```

Install `orthomap` into your current python environment:

```shell
pip install .
```
## Contributing Code

If you would like to contribute to `orthomap`, please file an issue so that one can establish a statement of need, avoid redundant work, and track progress on your contribution.

Before you do a pull request, you should always file an issue and make sure that someone from the `orthomap` developer team agrees that it's a problem, and is happy with your basic proposal for fixing it.

Once an issue has been filed and we've identified how to best orient your contribution with package development as a whole, [fork](https://docs.github.com/en/github/getting-started-with-github/fork-a-repo) the [main repo](https://github.com/kullrich/orthomap/orthomap.git), branch off a [feature branch](https://docs.github.com/en/github/collaborating-with-issues-and-pull-requests/about-branches) from `master`, [commit](https://docs.github.com/en/desktop/contributing-and-collaborating-using-github-desktop/committing-and-reviewing-changes-to-your-project) and [push](https://docs.github.com/en/github/using-git/pushing-commits-to-a-remote-repository) your changes to your fork and submit a [pull request](https://docs.github.com/en/github/collaborating-with-issues-and-pull-requests/proposing-changes-to-your-work-with-pull-requests) for `orthomap:master`.

By contributing to this project, you agree to abide by the Code of Conduct terms.

## Bug reports

Please report any errors or requests regarding [`orthomap`](https://github.com/kullrich/orthomap) to Kristian Ullrich (ullrich@evolbio.mpg.de)

or use the issue tracker at [https://github.com/kullrich/orthomap/issues](https://github.com/kullrich/orthomap/issues)

## Code of Conduct - Participation guidelines

This repository adhere to [Contributor Covenant](http://contributor-covenant.org) code of conduct for in any interactions you have within this project. (see [Code of Conduct](https://github.com/kullrich/orthomap/-/blob/master/CODE_OF_CONDUCT.md))

See also the policy against sexualized discrimination, harassment and violence for the Max Planck Society [Code-of-Conduct](https://www.mpg.de/11961177/code-of-conduct-en.pdf).

By contributing to this project, you agree to abide by its terms.


## References

Emms, D.M. and Kelly, S. (2019). **OrthoFinder: phylogenetic orthology inference for comparative genomics.** *Genome biology*, **20(1)**. [https://doi.org/10.1186/s13059-019-1832-y](https://doi.org/10.1186/s13059-019-1832-y)

Huerta-Cepas, J., Serra, F. and Bork, P. (2016). **ETE 3: reconstruction, analysis, and visualization of phylogenomic data.** *Molecular biology and evolution*, **33(6)**. [https://doi.org/10.1093/molbev/msw046](https://doi.org/10.1093/molbev/msw046)

Wolf, F.A., Angerer, P. and Theis, F.J. (2018). **SCANPY: large-scale single-cell gene expression data analysis.** *Genome biology*, **19(1)**. [https://doi.org/10.1186/s13059-017-1382-0](https://doi.org/10.1186/s13059-017-1382-0)

Qiu, C., Cao, J., Martin, B.K., Li, T., Welsh, I.C., Srivatsan, S., Huang, X., Calderon, D., Noble, W.S., Disteche, C.M. and Murray, S.A. (2022). **Systematic reconstruction of cellular trajectories across mouse embryogenesis.** *Nature genetics*, **54(3)**. [https://doi.org/10.1038/s41588-022-01018-x](https://doi.org/10.1038/s41588-022-01018-x)

[bioconda]: https://bioconda.github.io/
[from pypi]: https://pypi.org/project/orthomap
[miniconda]: http://conda.pydata.org/miniconda.html
[on github]: https://github.com/kullrich/orthomap