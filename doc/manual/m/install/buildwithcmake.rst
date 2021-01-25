Build Tinker9 with CMake
========================

Configure CMake
---------------
You can skip this section if you are familar with CMake.

Suppose the current working directory is */home/tinker9* and we
want to create a build directory called *build-cmake* in
*/home/tinker9*. We can do *mkdir build-cmake* then *cd build-cmake*.
Because the top-level CMakeLists.txt file is in the parent directory,
if there was nothing else to configure, command *cmake ..* would generate
the Makefile. The alternative way is to specify the build and source
directories to CMake, e.g.,

.. code-block:: bash

   cmake -B /home/tinker9/build-cmake -S /home/tinker9

Some CMake installations also provide a command line gui *ccmake* and a
simple gui program *cmake-gui* that can replace *cmake* in the commands
above.

Configure Compilers
-------------------
Set *CXX=...*, *CUDACXX=...*, and *FC=...* to specify the non-default C++,
CUDA, and Fortran compilers, respectively. These environmental variables
are supported by *cmake*.

This cmake script checks a custom environmental variable *ACC=...*
*only* for the GPU code.
If not set, the building script will take a guess at the OpenACC compiler.
It will be set to the default C++ compiler for the CPU code. For instance,
command *(c)cmake [...]* will become *ACC=pgc++ (c)cmake [...]*.

Configure Tinker9
-----------------
The following options are passed to CMake program with their default
values (if there is one). **-D** is prefixed to the options. CMake provides
two standard ways to let users to customize the values:

- Change their values interactively in the *ccmake* command line gui;
- Pass the new value to CMake via command line arguments
  *cmake -DOPTION=NewValue*.

In addition to these two canonical methods, default value can also be set
by its corresponding environmental variable, documented as **(env)** here.
Note that there is no **-D** prefix for the environmental variables.

Here are two equivalent examples to have Tinker9 configured as follows

=======================  ===================
Item                     Value
=======================  ===================
ACC                      pgc++
opt                      release
host                     0
prec                     m
cuda_dir                 /usr/local/cuda
compute_capability       75
tinker_dir               /home/tinker/source
fftw_dir                 /usr/local
CMakeLists.txt Location  /home/tinker9
=======================  ===================

.. code-block:: bash

   # use environmental variables
   ACC=pgc++ \
   opt=release host=0 prec=m \
   cuda_dir=/usr/local/cuda compute_capability=75 \
   tinker_dir=/home/tinker/source fftw_dir=/usr/local \
   cmake /home/tinker9

   # use cmake -DOPTIONS
   ACC=pgc++ cmake -S /home/tinker9 \
   -DCMAKE_BUILD_TYPE=Release -DHOST=0 -DPREC=m \
   -DCUDA_DIR=/usr/local/cuda -DCOMPUTE_CAPABILITY=75 \
   -DTINKER_DIR=/home/tinker/source -DFFTW_DIR=/usr/local


**-DCMAKE_BUILD_TYPE (opt) = release**

Build type is case insensitive and can either be *release* or *debug*.

**-DCMAKE_INSTALL_PREFIX (no env var) = [NO DEFAULT VALUE]**

Install the executables under *${CMAKE_INSTALL_PREFIX}*. If this option is
not set, *make install* is configured not to install anything, which is
different from the default cmake behavior to install the program under */usr/local*.

**-DTINKER_DIR (tinker_dir) = ${CMAKE_BINARY_DIR}/tinker/source**

Directory in which user compiled *libtinker.a*.

**-DFFTW_DIR (fftw_dir) = ${CMAKE_BINARY_DIR}/fftw**

Top-level FFTW3 installation, under which
*include/fftw3.h* and *lib/libfftw3* static libraries are expected to be found.

**-DHOST (host) = 1**

Flag to compile to GPU (with value 0 or OFF) or CPU (with value 1 or ON)
version.

**-DPREC (prec) = d**

Precision of the floating-point numbers. With flag *d*, all of the
floating-point numbers are treated as real*8/double values,
or real*4/single values if with flag *s*. Mixed precision flag *m* will
use real*4 or real*8 numbers in different places. Note that this flag will
not change the precision of the variables hard-coded as *float* or *double*
types.

**-DDETERMINISTIC_FORCE (deterministic_force) = [NO DEFAULT VALUE]**

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
`NVIDIA website. <https://developer.nvidia.com/cuda-gpus>`_

**-DCUDA_DIR (cuda_dir) = /usr/local/cuda**

Top-level CUDA installation directory, under which directories *include*,
*lib* or *lib64* can be found.
This option will supersede the CUDA installation identified by the official
*CUDACXX* environmental variable.

Sometimes the PGI compiler and the NVCC compiler are not "compatible." For
instance, although PGI 19.4 supports CUDA 9.2, 10.0, 10.1, but the default
CUDA version configured in PGI 19.4 may be 9.2 and the external NVCC version
is 10.1. One solution is to pass *CUDA_HOME=${cuda_dir}* to the PGI
compiler, in which case, **cuda_dir** should be set to
*/usr/local/cuda-10.1*.

Make Tinker9
------------
The following Makefile targets will be generated by CMake.
Run *make -j* for the default target(s) and *make TARGET(S) -j* for others.

**tinker9**

Compile and link the *tinker9* executable.

**all.tests**

Compile and link the *all.tests* executable.

**default**

Make two targets: *tinker9* and *all.tests* executables.

**all**

Same as the default target.

**test**

Run unit tests in a random order. Exit on the first error.

**man**

Generate user manual.

**doc**

Generate developer guides.

