Build the Canonical Tinker
==========================

**Checkout the Required Version of Tinker**

   **If the incorrect Tinker version was used, the executables would be
   very likely to fail with segfault.**

If this source code was cloned by git, you can
checkout Tinker from the *tinker* git submodule
and copy Tinker to the build directory:

.. code-block:: bash

   # create build directory
   mkdir -p tinker9/build

   # checkout Tinker
   cd tinker9/tinker
   git submodule update --init

   # copy Tinker
   cd ..
   cp -r tinker build

Alternatively, clone
`Tinker from GitHub <https://github.com/tinkertools/tinker>`_,
then checkout the required version **25b7ee7a**.

**Make libtinker**

- Copy the entire Tinker direcotry to the build directory.
- Copy Makefile (under *tinker/make*) to *tinker/source*.
- Change the value of variable *FFTWDIR* to the top-level FFTW installation.
- Locate the correct flags in the Makefile for your compiler and operating
  system, uncomment them, and comment out the others.
- **Must remove the OpenMP flag: -fopenmp, -qopenmp etc.**
- **Change optimization level (-O3, -Ofast etc.) to -O2.**
- Run *make* command

.. code-block:: bash

   # current working directory
   # /home/tinker9/build/tinker/source
   make libtinker.a -j
