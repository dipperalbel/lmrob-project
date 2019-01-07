.. _Installation:

System Requirements
===================

Additional to the repository the following packages must be in the system
before installing the package.

* ``lapack``
* ``blas``
* ``R``

Installation
============

Installing the required packages
    .. code-block:: bash

        $ sudo apt install r-base -y

Install a virtualenv
    .. code-block:: bash

        sudo apt install virtualenv


Create virtualev
    .. code-block:: bash

        $ virtualenv -p python3 /path/to/env

Activate the virtualenv
    .. code-block:: bash

        $ source /path/to/env/bin/activate


Installing lmrob amd nlrob

    .. code-block:: bash

        (env) $ python setup.py install_lib

Issues
------

The following compilers must be installed in the system:

* gcc
* gfortran

If the ``Python.h`` is not found install python-devel to the system

