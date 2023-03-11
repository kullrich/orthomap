.. _docker_additional_information:

orthomap docker image quick start and additional notes
======================================================

In this documentation, we provide example commands for our orthomap docker image below. But we are not aiming to provide detailed information about general docker usage.
We highly recommend learning about docker if you are not familiar with it, and make sure you have adequate knowledge about docker prior to start orthomap analysis with docker.

Quick start - Docker - bash
^^^^^^^^^^^^^^^^^^^^^^^^^^^

1. Download orthomap docker image from docker Hub.

::

    docker pull kkuweb/orthomap_ubuntu:latest

2. Make docker container and start running it.

- To access files smoothly, we can use `bind mount <https://docs.docker.com/storage/bind-mounts/>`_. Here, we will mount the `data_folder`, a directory in your local machine.

::

    mkdir data_folder # Create data folder in your local environment.
    docker run -dit \
      --name orthomap_container \
      -v $(pwd)/data_folder:/root/data_folder \
      orthomap_ubuntu:latest

3. Enter the docker container.

::

    docker container exec -it orthomap_container bash

4. In the docker container environment, start using orthomap bash scripts.

::

    cd
    cd data_folder
    orthomap -h
    orthomap qlin -q "Arabidopsis thaliana"

Quick start - Docker - jupyter
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

1. Download orthomap docker image from docker Hub.

::

    docker pull kkuweb/orthomap_ubuntu:latest

2. Make docker container and start running it.

- As we recommend using jupyter notebook (or jupyter lab), we need to set up port. Here, we connect docker container port 8888 with local host port 8888.

- To access files smoothly, we can use `bind mount <https://docs.docker.com/storage/bind-mounts/>`_. Here, we will mount the `data_folder`, a directory in your local machine.

::

    mkdir data_folder # Create data folder in your local environment.
    docker run -dit \
      --name orthomap_container \
      -p 8888:8888 \
      -v $(pwd)/data_folder:/mnt/data_folder \
      kkuweb/orthomap_ubuntu:latest

3. Enter the docker container.

::

    docker container exec -it orthomap_container bash

4. In the docker container environment, start jupyter notebook as follows.

::

    cd /mnt/data_folder
    jupyter notebook --port=8888 --ip=0.0.0.0 --allow-root --no-browser --notebook-dir=/mnt/data_folder

After starting jupyter, please open your browser and enter http://localhost:8888 to access jupyter notebook running in the docker container.
You need to enter a token to access jupyter. The token can be found in your terminal running jupyter notebook.

5. Start using orthomap from jupyter using the `orthomap_env` kernel.

::

    from orthomap import qlin
    qlin.get_qtid(q='Arabidopsis thaliana')

Docker image build information
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We built our docker image using Dockerfile automatic build function.
The Dockerfile is available `here <https://github.com/kullrich/orthomap/blob/main/docs/dockerfile>`_.
You can modify it to create custom docker image by yourself.
If you make custom environment, please do so on your responsibility.
