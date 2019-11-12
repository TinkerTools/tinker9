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
      - Tinker GPU is in sync with [Tinker](https://github.com/TinkerTools/Tinker) of commit `11e84c69`.
      For details, see `Build the Canonical Tinker`.
   - CUDA
      - Recent `CUDA` libraries are required by the GPU code,
        which may have been included in the PGI compiler.


## Other Nonmandatory Utilities
   - [clang-format](https://clang.llvm.org/docs/ClangFormat.html): to format the source code.
   - [doxygen](http://www.doxygen.nl): to generate the documentation.
