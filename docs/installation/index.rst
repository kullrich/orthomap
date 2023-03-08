.. _install:

Installation
============

Please follow the guide below to install orthomap and its dependent software.

.. _require:

Docker image
------------

- Pre-built docker image is available through `Docker Hub <https://hub.docker.com/repository/docker/kkuweb/orthomap_ubuntu>`_ .

::

    docker pull kkuweb/orthomap_ubuntu:latest

- This docker image was built based on Ubuntu 22.04.
- Python dependent packages and orthomap are installed in an anaconda environment, `orthomap_env`. This environment will be activated automatically when you log in.
- See additional information

.. toctree::
   :maxdepth: 1

   docker_additional_information

Singularity image
------------

- Pre-built docker image is available through `Docker Hub <https://hub.docker.com/repository/docker/kkuweb/orthomap_ubuntu>`_ .

::

    singularity pull kkuweb/orthomap_ubuntu:latest

- This docker image was built based on Ubuntu 22.04.
- Python dependent packages and orthomap are installed in an anaconda environment, `orthomap_env`. This environment needs to be activated when you log in.
- See additional information

.. toctree::
   :maxdepth: 1

   singularity_additional_information

Install orthomap
----------------

Python Requirements
^^^^^^^^^^^^^^^^^^^

- orthomap was developed using python 3.8. We do not support python 2.7x or python <=3.7.

orthomap installation using conda and pip
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

  We recommend installing orthomap in an independent conda environment to avoid dependent software conflicts.
  Please make a new python environment for orthomap and install dependent libraries in it.

  The environment is created with `conda create` in which orthomap is installed.

  If you do not have a working installation of Python 3.8 (or later), consider
  installing `Anaconda <https://docs.anaconda.com/anaconda/install/>`_ or `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_. Then run:

  ::

      git clone https://github.com/kullrich/orthomap.git
      cd orthomap
      conda env create --file environment.yml
      conda activate orthomap_env

  Install `orthomap` via `PyPI <https://pypi.org/project/orthomap>`_:

  ::

      pip install orthomap

Development Version
^^^^^^^^^^^^^^^^^^^

To work with the latest version `on GitHub <https://github.com/kullrich/orthomap>`_: clone the repository and `cd` into its root directory.

  ::

      git clone kullrich/orthomap
      cd orthomap

Install `orthomap` into your current python environment:

  ::

      pip install .

Installing Miniconda
^^^^^^^^^^^^^^^^^^^^

After downloading `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_, in a unix shell (Linux, Mac), run

  ::

      cd DOWNLOAD_DIR
      chmod +x Miniconda3-latest-VERSION.sh
      ./Miniconda3-latest-VERSION.sh

