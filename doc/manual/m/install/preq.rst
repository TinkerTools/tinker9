Prerequisites
=============

**Hardware**

A relatively recent Nvidia GPU is mandatory for the GPU code.
Nothing special is needed for the CPU code.

**Operating Systems and Compilers**

In order to compile the GPU code, the most recent
`PGI compiler <https://www.pgroup.com/products/community.htm>`_
for OpenACC directives is mandatory. Due to the limitations of the PGI compiler,
the GPU code is unavailable on macOS.

No effort has been spared on building excutables for Windows yet.

For Linux, we need:

- GNU or Intel Fortran compiler.
- Recent C++ compiler that supports C++11 syntax.
- (GPU code only) PGI compiler for OpenACC and nvcc for CUDA.

**FFTW Libraries**

Two prebuilt FFTW libraries: `libfftw3` and `libfftw3_threads` are used by
Tinker. Two more FFTW libraries, `libfftw3f`, `libfftw3f_threads` are
needed by the single precision CPU code.

**More About Using PGI Compiler on the Clusters**

   I recently (in Jan. 2020) worked on a cluster that was still running
   Red Hat with gcc 4.8.5 by default without root privilege. Although several
   more recent gcc and PGI versions were available via the `module` program,
   the most recent PGI compiler (2019) was still configured with gcc 4.8.5
   by default, which had a very bad support for C++11.
   Since I didn't have root privilege on the cluster, I had to use
   a custom `localrc` file by running the following command to
   reconfigure PGI compiler with gcc 7.4.0:

   .. code-block:: bash

      makelocalrc $PGI/linux86-64-llvm/2019 \
      -gcc /usr/local/gcc-7.4.0/bin/gcc \
      -gpp /usr/local/gcc-7.4.0/bin/g++ \
      -g77 /usr/local/gcc-7.4.0/bin/gfortran \
      -o -net > /path/to/new_config

   then added `export PGI_LOCALRC=/path/to/new_config` to my bash resource file.

**Other Nonmandatory Utilities**

- `ClangFormat <https://clang.llvm.org/docs/ClangFormat.html>`_:
  to format the source code.
- `Sphinx <https://www.sphinx-doc.org>`_ and tex: to generate user manual.
- `Doxygen <https://www.doxygen.nl>`_: to generate developers' manual.
