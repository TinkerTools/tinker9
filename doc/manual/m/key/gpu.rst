GPU, OpenACC, CUDA Keywords
===========================

**CUDA-DEIVCE [integer]**
Followed by an integer value starting from 0, this keyword sets the CUDA-enabled
GPU device for the program. Value will be overwritten by environment variable
``cuda_device``.

For instance, a node has four CUDA devices, and the ``CUDA_VISIBLE_DEVICES``
environment variable (part of CUDA library) has been set to
``CUDA_VISIBLE_DEVICES=1,3``. This means only two CUDA devices are avaiable
here, thus the valid values for the **CUDA-DEVICE** are 0 and 1.

**GPU-PACKAGE [CUDA/OpenACC]**
This keyword selects code paths of some GPU algorithms where both CUDA and
OpenACC versions have been implemented.

The default value is ``CUDA``. Value will be overwritten by environment variable
``gpu_package``.
