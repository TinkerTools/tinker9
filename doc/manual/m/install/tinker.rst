Download the Canonical Tinker
=============================

**Recommended**

Using the incorrect Tinker version, the executables would be
very likely to fail with segfault.

If this source code was cloned by Git, you can
checkout Tinker from the *tinker* Git submodule:

.. code-block:: bash

   # checkout Tinker
   cd tinker9
   git submodule update --init

**Deprecated**

Alternatively, remove the directory *tinker9/tinker* and clone
`Tinker from GitHub <https://github.com/tinkertools/tinker>`_
to replace the deleted directory,
then checkout the required version **080e8f1d**.

.. code-block:: bash

   cd tinker9
   rm -rf tinker
   git clone https://github.com/tinkertools/tinker
   cd tinker
   git checkout <TheRequiredVersion>
