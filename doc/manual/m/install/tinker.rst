Build the Canonical Tinker
==========================

**Checkout the Required Version of Tinker**

Using the incorrect Tinker version, the executables would be
very likely to fail with segfault.

If this source code was cloned by Git, you can
checkout Tinker from the *tinker* Git submodule:

.. code-block:: bash

   # create build directory
   mkdir -p tinker9/build

   # checkout Tinker
   cd tinker9/tinker
   git submodule update --init

Alternatively, clone
`Tinker from GitHub <https://github.com/tinkertools/tinker>`_,
then checkout the required version **6a1c6104**.
You should move this *tinker* directory under *tinker9*.

**Make libtinker**

Compiling libtinker is now automated in the next step.
