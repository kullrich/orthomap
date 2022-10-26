# Installation

## Anaconda

The environment is created with `conda create` in which orthomap is installed.

If you do not have a working installation of Python 3.7 (or later), consider
installing [Miniconda] (see [Installing Miniconda]). Then run:

```shell
conda env create --file environmen.yml
conda activate orthomap
```

Install `orthomap`:

```shell
pip install orthomap
```

## PyPI

Install `orthomap` into your current python environment:

```shell
pip install orthomap
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

## Installing Miniconda

After downloading [Miniconda], in a unix shell (Linux, Mac), run

```shell
cd DOWNLOAD_DIR
chmod +x Miniconda3-latest-VERSION.sh
./Miniconda3-latest-VERSION.sh
```

[bioconda]: https://bioconda.github.io/
[from pypi]: https://pypi.org/project/orthomap
[miniconda]: http://conda.pydata.org/miniconda.html
[on github]: https://github.com/kullrich/orthomap