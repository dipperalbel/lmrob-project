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

        $ sudo apt install virtualenv


Create virtualev. Avoid to use the command sudo before virtualenv.
    .. code-block:: bash

        $ virtualenv -p python3 /path/to/env

Activate the virtualenv
    .. code-block:: bash

        $ source /path/to/env/bin/activate

Installing the python modules required for lmrob

    .. code-block:: bash

        (env) $ python setup.py install

Installing lmrob

    .. code-block:: bash

        (env) $ python setup.py install_lib

Issues
------

The following compilers must be installed in the system:

* gcc
* gfortran

If the ``Python.h`` is not found install python-devel to the system.

Installing python-devel

    .. code-block:: bash

	$ sudo apt install python3-dev

