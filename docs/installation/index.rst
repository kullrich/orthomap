.. _install:

Installation
============

Please follow the guide below to install orthomap and its dependent software.

.. _require:

Docker image
------------

- Pre-built docker image is available through `Docker Hub <https://hub.docker.com/repository/docker/kenjikamimoto126/celloracle_ubuntu>`_ .

::

    docker pull kenjikamimoto126/celloracle_ubuntu:latest


- This docker image was built based on Ubuntu 20.04.
- Python dependent packages and celloracle are installed in an anaconda environment, celloracle_env. This environment will be activated automatically when you log in.

.. toctree::
   :maxdepth: 1

   docker_additional_information

Install CellOracle
------------------

## Anaconda

The environment is created with `conda create` in which orthomap is installed.

If you do not have a working installation of Python 3.8 (or later), consider
installing [Miniconda] (see [Installing Miniconda]). Then run:

```shell
conda env create --file environment.yml
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