Build the Canonical Tinker
==========================

**Checkout the Required Version of Tinker**

If source code was `cloned` by `git`, you can
checkout Tinker from the `tinker` git submodule
and copy Tinker to the build directory:

.. code-block:: bash

   # create build directory
   mkdir -p tinker.gpu/build
   # checkout Tinker
   cd tinker.gpu/tinker
   git submodule update --init
   # copy Tinker
   cd ..
   cp -r tinker build

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

   # current working directory
   # /home/tinker.gpu/build/tinker/source
   make libtinker.a -j
