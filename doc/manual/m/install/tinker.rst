Build the Canonical Tinker
==========================

**Checkout the Required Version of Tinker**

Build Tinker from the git submodule inside this repository.

Alternatively, clone
`Tinker from GitHub <https://github.com/tinkertools/tinker>`_,
checkout the required version 350df099,
and copy Makefile to the `source` directory.

**Make libtinker**

- Change the value of variable `FFTWDIR` to the top-level FFTW installation.
- Locate the correct flags in the Makefile for your compiler and operating
  system, uncomment them, and comment out the others.
- **Must disable the OpenMP flag.**
- **Change optimization level to -O2.**
- Run `make` command

.. code-block:: bash

   make libtinker.a -j
