# Prerequisites


## Hardware
A relatively recent Nvidia GPU is mandatory for the GPU code.
Nothing special is needed for the CPU code.


## Compilers and Operating Systems
   - A GNU or Intel Fortran compiler.
   - A recent C++ compiler that supports C++11 syntax.
   - The most recent [PGI compiler](https://www.pgroup.com/products/community.htm)
     for the `OpenACC` directives.
   - Linux and Windows 10 (Windows Subsystem for Linux) are preferred because
     the PGI compiler does not support the GPU code on macOS.


## Libraries
   - Prebuilt FFTW libraries are required for
      - Tinker: `libfftw3` and `libfftw3_threads`.
      - Tinker GPU (single precision CPU version): `libfftw3f`, `libfftw3f_threads`.
   - Tinker
      - Tinker GPU is in sync with [Tinker](https://github.com/TinkerTools/Tinker) of commit `350df099`.
      For details, see `Build the Canonical Tinker`.
   - CUDA
      - Recent `CUDA` libraries are required by the GPU code,
        which may have been included in the PGI compiler.


## More About Using PGI Compiler on the Clusters
I recently (in Jan. 2020) worked on a cluster that was still running
Red Hat with gcc 4.8.5 by default without root privilege. Although several
more recent gcc and PGI versions were available via the `module` program,
the most recent PGI compiler (2019) was still configured with gcc 4.8.5
by default. `makelocalrc` command from PGI compiler required root privilege
to let it cooperate with newer gcc versions.


There is a way to use custom `localrc` file. I run the following command to
reconfigure PGI compiler with gcc 7.4.0.
```
makelocalrc $PGI/linux86-64-llvm/2019 \
-gcc /usr/local/gcc-7.4.0/bin/gcc \
-gpp /usr/local/gcc-7.4.0/bin/g++ \
-g77 /usr/local/gcc-7.4.0/bin/gfortran \
-o -net > /path/to/new_config
```
then added `export PGI_LOCALRC=/path/to/new_config` in my bash resource file.


## Other Nonmandatory Utilities
   - [clang-format](https://clang.llvm.org/docs/ClangFormat.html): to format the source code.
   - [doxygen](http://www.doxygen.nl): to generate the documentation.
