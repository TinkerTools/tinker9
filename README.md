Tinker GPU: Software Tools for Molecular Design
===============================================
[//]: # (Badges)
[![Build Status](https://travis-ci.com/zhi-wang/tinker.gpu.svg?branch=master)](https://travis-ci.com/zhi-wang/tinker.gpu)


## What I Am Up To
### Build System
<del>Converting to cmake is in progress.</del>
CMake scripts have been added. Tests are ongoing.

**Deprecated:**
In the meanwhile, please always run `make -f make/Makefile create_build` to create a build directory first
in order to generate other files in the build directory required by `Makefile`.

### GNU C++
The version of `g++` cannot be overly recent, in case `nvcc` might fail to support.
I also believe version `4.x` is outdated for some of the C++11 features used here.

## Installation Guide
   - [Prerequisites](doc/manual/m/install/preq.rst)
   - [Build the Canonical Tinker](doc/manual/m/install/tinker.rst)
   - [(Deprecated) Build Tinker GPU with Makefile](doc/manual/m/install/tinkergpu.rst)
   - [(New) Build Tinker GPU with CMake](doc/manual/m/install/buildwithcmake.rst)


## Repository
<a href="https://github.com/zhi-wang/tinker.gpu">
   <img src="https://github.githubassets.com/images/modules/logos_page/GitHub-Mark.png"
   alt="GitHub Repository" width="100"/>
</a>


For more information, visit the user manual
hosted on [github.io](https://zhi-wang.github.io/tinker.gpu).
