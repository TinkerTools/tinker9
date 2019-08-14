# Prerequisites

# Table of Contents
* [Hardware](#hardware)
* [Operating Systems and Compilers](#os)
* [Libraries](#lib)

<a name='hardware'></a>
# Hardware
A relatively recent Nvidia GPU is mandatory to compile and execute the code.

The code can also run on CPU, in which case nothing special is needed.

<a name='os'></a>
# Operating Systems and Compilers

## Linux
* a recent GNU C++ compiler that supports `-std=c++11` flag
* the most recent [PGI C++ compiler](https://www.pgroup.com/products/community.htm)
for the `OpenACC` directives
* a recent `cuda` library, may have been included in the PGI C++ compiler

## macOS

## Windows

<a name='lib'></a>
# Libraries

## FFTW
Prebuilt FFTW libraries are required for
* Tinker: `libfftw3` and `libfftw3_threads`
* Tinker GPU (single precision CPU version): `libfftw3f`, `libfftw3f_threads`

## Tinker
Tinker GPU is in sync with
[Tinker](https://github.com/TinkerTools/Tinker/tree/291a85c1435feddc835e80bfa340497b67cc1393)
of commit `291a85c1`.
For details, see section `Build and Install`.
