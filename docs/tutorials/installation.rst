Installation
============

**PyPi**

The simplest way to install the package is to obtain it from the PyPi repository::

    $ pip install orthomap

**Source**

You can also build `orthomap` from source. To clone the repository run::

    $ git clone https://github.com/kullrich/orthomap.git

The repository contains a `conda` environment file with all the requirements.
The environment is created with `conda create`::

    $ cd orthomap
    $ conda env create --file environment.yml
    $ conda activate orthomap

To install `orthomap` run::

    $ pip install .

or it install the package in developement mode::

    $ pip install -e .
