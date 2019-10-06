# Prerequisites


## Hardware
A relatively recent Nvidia GPU is mandatory for the GPU code.
Nothing special is needed for the CPU platform.


## Compilers and Operating Systems
   - A GNU or Intel Fortran compiler.
   - A recent C++ compiler that supports C++11 syntax. 
   - The most recent [PGI compiler](https://www.pgroup.com/products/community.htm) 
     for the `OpenACC` directives.
   - Linux and Windows 10 (Windows Subsystem for Linux) are preferred because
     the PGI compiler does not support the GPU platform on macOS.


## Libraries
   - Prebuilt FFTW libraries are required for
      - Tinker: `libfftw3` and `libfftw3_threads`.
      - Tinker GPU (single precision CPU version): `libfftw3f`, `libfftw3f_threads`.
   - Tinker
      - Tinker GPU is in sync with [Tinker](https://github.com/TinkerTools/Tinker/tree/291a85c1435feddc835e80bfa340497b67cc1393)
        of commit `291a85c1`. For details, see `Build the Canonical Tinker`.
   - CUDA
      - Recent `CUDA` libraries are required the GPU platform,
        which may have been included in the PGI compiler.
