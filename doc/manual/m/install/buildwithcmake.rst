Build Tinker GPU with CMake
===========================

Configure CMake
---------------
This section can be skipped for those who are familar with CMake.

Suppose the current working directory is "/home/tinker.gpu" and we
want to create a build directory called "build-cmake" in
"/home/tinker.gpu". We can do "mkdir build-cmake" then "cd build-cmake".
Because the top-level CMakeLists.txt file is in the parent directory,
if there were nothing else to configure, command "cmake .." would generate
the Makefile. The alternative way is to specify the build and source
directories to CMake, e.g.,

.. code-block:: bash

   cmake -B /home/tinker.gpu/build-cmake -S /home/tinker.gpu

There is also a command line gui "ccmake" that can replace "cmake" in the
commands above.

Configure Compilers
-------------------
Set "CXX=..." and "FC=..." to use the non-default C++ and Fortran compilers,
respectively.

In general there is only one situation where you have to configure the C++
compiler in command line, regardless of using "cmake" or "ccmake", which
is to compile the GPU OpenACC code with "pgc++". The command
"(c)cmake [...]" will become "CXX=pgc++ (c)cmake [...]"

Configure Tinker GPU
--------------------
The following options are passed to CMake program with their default
values if they have one. **-D** is prefixed to the options. To change
the default value of the options, there are two ways.
You can always change their values interactively in the "ccmake" command
line gui. You can also pass the new value to "cmake" via command line
argument "cmake -DOPTION=NewValue".

In addition to these two canonical methods, default values can also be set
by the corresponding environmental variables, documented as **(env)**.
Note that there is no **-D** prefix for the environmental variables.

Here are two equivalent examples to have Tinker GPU configured as follows

=======================  ===================
Item                     Value
=======================  ===================
CXX                      pgc++
opt                      release
host                     0
prec                     m
cuda_dir                 /usr/local/cuda
compute_capability       75
tinker_dir               /home/tinker/source
fftw_dir                 /usr/local
CMakeLists.txt Location  /home/tinker.gpu
=======================  ===================

.. code-block:: bash

   # use environmental variables
   CXX=pgc++ \
   opt=release host=0 prec=m \
   cuda_dir=/usr/local/cuda compute_capability=75 \
   tinker_dir=/home/tinker/source fftw_dir=/usr/local \
   cmake /home/tinker.gpu

   # use cmake -D option
   CXX=pgc++ cmake -S /home/tinker.gpu \
   -DCMAKE_BUILD_TYPE=Release -DHOST=0 -DPREC=m \
   -DCUDA_DIR=/usr/local/cuda -DCOMPUTE_CAPABILITY=75 \
   -DTINKER_DIR=/home/tinker/source -DFFTW_DIR=/usr/local


**-DCMAKE_BUILD_TYPE (opt) = release**
Build type is case insensitive and can either be "release" or "debug".

**-DTINKER_DIR (tinker_dir) = ${CMAKE_BINARY_DIR}/tinker/source**
Directory in which user compiled "libtinker.a".

**-DFFTW_DIR (fftw_dir) = ${CMAKE_BINARY_DIR}/fftw**
Top-level FFTW3 installation, under which
"include/fftw3.h" and "lib/libfftw3" are expected to be found.

**-DHOST (host) = 1**
Flag to compile to GPU (with value 0 or OFF) or CPU (with value 1 or ON)
version.

**-DPREC (prec) = d**
Precision of the floating-point numbers. With flag "d", all of the
floating-point numbers are treated as real*8/double values,
or real*4/single values if with flag "s". Mixed precision flag "m" will
use real*4 or real*8 numbers in different places. Note that this flag will
not change the hard-coded types.

**-DDETERMINISTIC_FORCE (deterministic_force) = [NO DEFAULT]**
Flag to use deterministic force. There is no default value for this flag.
This feature will be implicitly enabled by mixed and single precisions, but
can be explicitly disabled by setting the flag to 0,
and can be explicitly enabled by value 1.

In general, evaluating energy, forces etc. twice, we don't expect to get
two identical answers, but we may not care as much because the difference
is usually negligible. (See
`Why is cos(x) != cos(y)? <https://isocpp.org/wiki/faq/newbie#floating-point-arith2>`_)
Whereas in MD, two simulations with the same initial configurations can
easily diverge due to the accumulated difference. If, for whatever reason,
you are willing to elongate the process of the inevitable divergence at the
cost of slightly slower simulation speed, a more "deterministic" force
(using fixed-point arithmetic) can help.

**-DCOMPUTE_CAPABILITY (compute_capability) = 60,70**
CUDA compute capability (multiplied by 10) of GPU.
Valid values (noninclusive) are 35, 50, 60, 70, 75 etc., and can be
comma-separated, e.g. 35,60.
Multiple compute capabilites will increase the size of executables.

The full list of compute capabilities can be found on the
`Nvidia website. <https://developer.nvidia.com/cuda-gpus>`_

**-DCUDA_DIR (cuda_dir) = /usr/local/cuda**
Top-level CUDA installation directory, under which "include"
directory can be found.

Sometimes the PGI compiler and the NVCC compiler are not "compatible". For
instance, PGI 19.4 supports CUDA 9.2, 10.0, 10.1, and the default CUDA
version used in PGI 19.4 may be 9.2 and the external NVCC version is 10.1.
One solution is to pass "CUDA_HOME=${cuda_dir}" to the PGI compiler, in
which case, **cuda_dir** should be set to "/usr/local/cuda-10.1".

Make Tinker GPU
---------------
The following targets will be available in the Makefile generated by CMake.
Run "make -j" for the default target(s) and "make TARGET(S) -j" for others.

**tinker.gpu**
Compile and link the tinker.gpu executable.

**all.tests**
Compile and link the all.tests executable.

**default**
Make two targets: tinker.gpu and all.tests executables.

**all**
Same as the default target.

**test**
Run unit tests in a random order. Exit on the first error.

**man**
Generate user's manual.

**doc**
Generate developer's manual.

