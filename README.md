# orthomap

[![GitHub Workflow Status](https://img.shields.io/github/actions/workflow/status/kullrich/orthomap/build_check.yml?branch=main)](https://github.com/kullrich/orthomap/actions/workflows/build_check.yml)
[![PyPI](https://img.shields.io/pypi/v/orthomap?color=blue)](https://pypi.org/project/orthomap/)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/orthomap)](https://pypi.org/project/orthomap/)
[![PyPI - Wheel](https://img.shields.io/pypi/wheel/orthomap)](https://pypi.org/project/orthomap/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![docs-badge](https://readthedocs.org/projects/orthomap/badge/?version=latest)](https://orthomap.readthedocs.io/en/latest/?badge=latest)

## orthologous maps - evolutionary age index

[`orthomap`](https://github.com/kullrich/orthomap) is a python package to extract orthologous maps
(in other words the evolutionary age of a given orthologous group) from OrthoFinder results.
Orthomap results (gene ages per orthogroup) can be further used to calculate and visualize weighted expression data
from scRNA sequencing objects.

## Documentation

Online documentation can be found [here](https://orthomap.readthedocs.io/en/latest/).

## Installing `orthomap`

More installation options can be found [here](https://orthomap.readthedocs.io/en/latest/installation/index.html).

### orthomap installation using conda and pip

We recommend installing `orthomap` in an independent conda environment to avoid dependent software conflicts.
Please make a new python environment for `orthomap` and install dependent libraries in it.

The environment is created with conda create in which `orthomap` is installed.

If you do not have a working installation of Python 3.8 (or later),
consider installing [Anaconda](https://docs.anaconda.com/anaconda/install/) or
[Miniconda](https://docs.conda.io/en/latest/miniconda.html). Then run:

```shell
$ git clone https://github.com/kullrich/orthomap.git
$ cd orthomap
$ conda env create --file environment.yml
$ conda activate orthomap_env
```

Install `orthomap` via [PyPI](Install orthomap via PyPI:):

```shell
$ pip install orthomap
```

## Quick usage

Detailed tutorials how to use `orthomap` can be found [here](https://orthomap.readthedocs.io/en/latest/tutorials/index.html).

### Update/download local ncbi taxonomic database:

The following command downloads or updates your local copy of the
NCBI's taxonomy database (~300MB). The database is saved at
`~/.etetoolkit/taxa.sqlite`.

```python
>>> from orthomap import ncbitax
>>> ncbitax.update_ncbi()
```

### Step 1 - Get query species taxonomic lineage information:

You can query a species lineage information based on its name or its
taxID. For example `Danio rerio` with taxID `7955`:

```python
>>> from orthomap import qlin
>>> qlin.get_qlin(q = 'Danio rerio')
>>> qlin.get_qlin(qt = '7955')
```

You can get the query species topology as a tree.
For example for `Danio rerio` with taxID `7955`:

```python
>>> from orthomap import qlin
>>> query_topology = qlin.get_lineage_topo(qt = '7955')
>>> query_topology.write()
```

### Step 2 - Get query species orthomap from OrthoFinder results:

The following code extracts the `orthomap` for `Danio rerio` based on pre-calculated 
OrthoFinder results and ensembl release-105:

OrthoFinder results (-S diamond_ultra_sens) using translated, longest-isoform coding sequences
from ensembl release-105 have been archived and can be found
[here](https://zenodo.org/record/7242264#.Y1p19i0Rowc).

```python
>>> from orthomap import datasets, of2orthomap
>>> datasets.ensembl105(datapath='.')
>>> query_orthomap = of2orthomap.get_orthomap(
...     seqname='Danio_rerio.GRCz11.cds.longest',
...     qt='7955',
...     sl='ensembl_105_orthofinder_species_list.tsv',
...     oc='ensembl_105_orthofinder_Orthogroups.GeneCount.tsv',
...     og='ensembl_105_orthofinder_Orthogroups.tsv',
...     out=None, quiet=False, continuity=True, overwrite=True)
>>> query_orthomap
```

### Step 3 - Map OrthoFinder gene names and scRNA gene/transcript names:

The following code extracts the gene to transcript table for `Danio rerio`:

GTF file obtained from [here](https://ftp.ensembl.org/pub/release-105/gtf/danio_rerio/Danio_rerio.GRCz11.105.gtf.gz).

```python
>>> from orthomap import datasets, gtf2t2g
>>> gtf_file = datasets.zebrafish_gtf(datapath='.')
>>> query_species_t2g = gtf2t2g.parse_gtf(
...     gtf=gtf_file,
...     g=True, b=True, p=True, v=True, s=True, q=True)
>>> query_species_t2g
```

Import now, the scRNA dataset of the query species.

example: **Danio rerio** - [http://tome.gs.washington.edu](http://tome.gs.washington.edu)
([Qui et al. 2022](https://www.nature.com/articles/s41588-022-01018-x))

`AnnData` file can be found [here](https://doi.org/10.5281/zenodo.7243602).

```python
>>> import scanpy as sc
>>> from orthomap import datasets, orthomap2tei
>>> # download zebrafish scRNA data here: https://doi.org/10.5281/zenodo.7243602
>>> # or download with datasets.qiu22_zebrafish(datapath='.')
>>> zebrafish_data = datasets.qiu22_zebrafish(datapath='.')
>>> zebrafish_data
>>> # check overlap of transcript table <gene_id> and scRNA data <var_names>
>>> orthomap2tei.geneset_overlap(zebrafish_data.var_names, query_species_t2g['gene_id'])
```

The `replace_by` helper function can be used to add a new column to the `orthomap` dataframe by matching e.g.
gene isoform names and their corresponding gene names.

```python
>>> # convert orthomap transcript IDs into GeneIDs and add them to orthomap
>>> query_orthomap['geneID'] = orthomap2tei.replace_by(
...    x_orig = query_orthomap['seqID'],
...    xmatch = query_species_t2g['transcript_id_version'],
...    xreplace = query_species_t2g['gene_id'])
>>> # check overlap of orthomap <geneID> and scRNA data
>>> orthomap2tei.geneset_overlap(zebrafish_data.var_names, query_orthomap['geneID'])
```

### Step 4 - Get transcriptome evolutionary index (TEI) values and add them to scRNA dataset:

Since now the gene names correspond to each other in the `orthomap` and the scRNA adata object,
one can calculate the transcriptome evolutionary index (TEI) and add them to the scRNA dataset (adata object).

```python
>>> # add TEI values to existing adata object
>>> orthomap2tei.get_tei(adata=zebrafish_data,
...    gene_id=query_orthomap['geneID'],
...    gene_age=query_orthomap['PSnum'],
...    keep='min',
...    layer=None,
...    add=True,
...    obs_name='tei',
...    boot=False,
...    bt=10,
...    normalize_total=False,
...    log1p=False,
...    target_sum=1e6)
```

### Step 5 - Downstream analysis

Once the gene age data has been added to the scRNA dataset,
one can e.g. plot the corresponding transcriptome evolutionary index (TEI) values
by any given observation pre-defined in the scRNA dataset.

#### Boxplot TEI per stage:

```python
>>>sc.pl.violin(adata=zebrafish_data,
...             keys=['tei'],
...             groupby='stage',
...             rotation=90,
...             palette='Paired',
...             stripplot=False,
...             inner='box')
```

## orthomap via Command Line

`orthomap` can also be used via the command line.

Command line documentation can be found [here](https://orthomap.readthedocs.io/en/latest/modules/orthomap.html).

```shell
$ orthomap
```

```
usage: orthomap <sub-command>

orthomap

optional arguments:
  -h, --help            show this help message and exit

sub-commands:
  {cds2aa,gtf2t2g,ncbitax,of2orthomap,qlin}
                        sub-commands help
    cds2aa              translate CDS to AA and optional retain longest
                        isoform <cds2aa -h>
    gtf2t2g             extracts transcript to gene table from GTF <gtf2t2g
                        -h>
    ncbitax             update local ncbi taxonomy database <ncbitax -h>
    of2orthomap         extract orthomap from OrthoFinder output for query
                        species <orthomap -h>
    qlin                get query lineage based on ncbi taxonomy <qlin -h>
```

To retrieve e.g. the lineage information for `Danio rerio` run the following command:

```shell
$ orthomap qlin -q "Danio rerio"
```

## Development Version

To work with the latest version [on GitHub](https://github.com/kullrich/orthomap):
clone the repository and `cd` into its root directory.

```shell
$ git clone kullrich/orthomap
$ cd orthomap
```

Install `orthomap` into your current python environment:

```shell
$ pip install -e .
```

## Contributing Code

If you would like to contribute to `orthomap`, please file an issue so that one can establish a statement of need, avoid redundant work, and track progress on your contribution.

Before you do a pull request, you should always file an issue and make sure that someone from the `orthomap` developer team agrees that it's a problem, and is happy with your basic proposal for fixing it.

Once an issue has been filed and we've identified how to best orient your
contribution with package development as a whole,
[fork](https://docs.github.com/en/github/getting-started-with-github/fork-a-repo)
the [main repo](https://github.com/kullrich/orthomap/orthomap.git), branch off a
[feature
branch](https://docs.github.com/en/github/collaborating-with-issues-and-pull-requests/about-branches)
from `master`,
[commit](https://docs.github.com/en/desktop/contributing-and-collaborating-using-github-desktop/committing-and-reviewing-changes-to-your-project)
and
[push](https://docs.github.com/en/github/using-git/pushing-commits-to-a-remote-repository)
your changes to your fork and submit a [pull
request](https://docs.github.com/en/github/collaborating-with-issues-and-pull-requests/proposing-changes-to-your-work-with-pull-requests)
for `orthomap:master`.

By contributing to this project, you agree to abide by the Code of Conduct terms.

## Bug reports

Please post troubles or questions on the GitHub repository [issue tracker](https://github.com/kullrich/orthomap/issues).
Also, please look at the closed issue pages. This might give an answer to your question.

## Inquiry for collabolation or discussion

Please send e-mail to us if you want a discussion with us.

Principal code developer: Kristian Ullrich

E-mail address can be found [here](https://www.evolbio.mpg.de).

## Code of Conduct - Participation guidelines

This repository adheres to the [Contributor Covenant](http://contributor-covenant.org) code of conduct for in any interactions you have within this project. (see [Code of Conduct](https://github.com/kullrich/orthomap/-/blob/master/CODE_OF_CONDUCT.md))

See also the policy against sexualized discrimination, harassment and violence for the Max Planck Society [Code-of-Conduct](https://www.mpg.de/11961177/code-of-conduct-en.pdf).

By contributing to this project, you agree to abide by its terms.

## References

Emms, D.M. and Kelly, S. (2019). **OrthoFinder: phylogenetic orthology inference for comparative genomics.** *Genome biology*, **20(1)**. [https://doi.org/10.1186/s13059-019-1832-y](https://doi.org/10.1186/s13059-019-1832-y)

Huerta-Cepas, J., Serra, F. and Bork, P. (2016). **ETE 3: reconstruction, analysis, and visualization of phylogenomic data.** *Molecular biology and evolution*, **33(6)**. [https://doi.org/10.1093/molbev/msw046](https://doi.org/10.1093/molbev/msw046)

Wolf, F.A., Angerer, P. and Theis, F.J. (2018). **SCANPY: large-scale single-cell gene expression data analysis.** *Genome biology*, **19(1)**. [https://doi.org/10.1186/s13059-017-1382-0](https://doi.org/10.1186/s13059-017-1382-0)

Qiu, C., Cao, J., Martin, B.K., Li, T., Welsh, I.C., Srivatsan, S., Huang, X., Calderon, D., Noble, W.S., Disteche, C.M. and Murray, S.A. (2022). **Systematic reconstruction of cellular trajectories across mouse embryogenesis.** *Nature genetics*, **54(3)**. [https://doi.org/10.1038/s41588-022-01018-x](https://doi.org/10.1038/s41588-022-01018-x)
