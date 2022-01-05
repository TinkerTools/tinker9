Download the Canonical Tinker
=============================

Since **d71f4793** (`Dec. 6, 2021 <https://github.com/TinkerTools/tinker9/commit/d71f4793>`_),
downloading the required Tinker version is automated in the CMake script.
For versions prior to this commit, please refer to the following paragraphs.

**Deprecated 1**

Using the incorrect Tinker version, the executables would be
very likely to fail with segfault.

If this source code was cloned by Git, you can
checkout Tinker from the *tinker* Git submodule:

.. code-block:: bash

   # checkout Tinker
   cd tinker9
   git submodule update --init

**Deprecated 2**

Alternatively, remove the directory *tinker9/tinker* and clone
`Tinker from GitHub <https://github.com/tinkertools/tinker>`_
to replace the deleted directory,
then checkout the required version **5aa9948d**.

.. code-block:: bash

   cd tinker9
   rm -rf tinker
   git clone https://github.com/tinkertools/tinker
   cd tinker
   git checkout <TheRequiredVersion>
