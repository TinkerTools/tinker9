Build the Canonical Tinker
==========================

**Checkout the Required Version of Tinker**

Checkout Tinker from the `tinker` git submodule inside this repository:

.. code-block:: bash

   cd tinker.gpu/tinker
   git submodule update --init

Alternatively, clone
`Tinker from GitHub <https://github.com/tinkertools/tinker>`_,
then checkout the required version 350df099.

**Make libtinker**

- Copy Makefile to the `source` directory.
- Change the value of variable `FFTWDIR` to the top-level FFTW installation.
- Locate the correct flags in the Makefile for your compiler and operating
  system, uncomment them, and comment out the others.
- **Must disable the OpenMP flag.**
- **Change optimization level to -O2.**
- Run `make` command

.. code-block:: bash

   make libtinker.a -j
