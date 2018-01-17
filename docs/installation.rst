Installing CrysFiPy
====================
There are two ways how to install CrysFiPy package.

Python Package Index - pip
--------------------------
Easier is to use `pip <https://pip.pypa.io/en/latest/installing.html>`_.

If you do not already have ``pip``, to install it first download `get-pip.py <https://bootstrap.pypa.io/get-pip.py>`_ and run it with the following command::

    python get-pip.py

With ``pip`` installed, you can install the latest version of crysfipy with the command::

    pip install crysfipy

New releases will be pushed to the package index automatically. If you wish to install the development version, you will need to follow the instructions for installation from source.

Installation from Source
------------------------
To install from source, either download the `master branch source from Bitbucket <https://bitbucket.org/cermak/crysfipy/get/master.zip>`_ or clone the repository::

    git clone https://bitbucket.org/cermak/crysfipy.git

From inside of the crysfipy directory install the package using::

    python setup.py install