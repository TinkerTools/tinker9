Tinker 9: Software Tools for Molecular Design
===============================================
[//]: # (Badges)
[![Build Status](https://travis-ci.com/zhi-wang/tinker.gpu.svg?branch=master)](https://travis-ci.com/zhi-wang/tinker.gpu)

# This repository will be renamed to `tinker9` and move to [TinkerTools](https://github.com/tinkertools) in a few days.

## What I Am Up To
### Build System
CMake scripts have been added.

**(Deprecated) Makefile**
Please always run `make -f make/Makefile create_build` to create a build
directory first in order to generate other files in the build directory
required by `Makefile`.

### GNU C++
The version of `g++` cannot be overly recent, in case `nvcc` might fail to support.
I also believe version `4.x` is outdated for some of the C++11 features used here.


## Getting the Source Code
If you don't plan to send us pull requests on Github, and you don't want to get
the entire history of the source code, you can clone a few recent versions by
```bash
git clone --depth 5 https://github.com/zhi-wang/tinker.gpu
```

To get the full git history,
```bash
git clone https://github.com/zhi-wang/tinker.gpu
```

If you are planning to send us pull requests, please fork this project to your
personal Github account.

Directly downloading zip file from Github webpage would work but is not
recommended because of the following cons:
   * No Git commit information to keep track of issues or bugs.
   * Tinker submodule would not work.
   * Has to download zip file in the furture for new commits.


## Installation Steps
   1. [Prerequisites](doc/manual/m/install/preq.rst)
   2. [Build the Canonical Tinker (important to get THE REQUIRED version of Tinker)](doc/manual/m/install/tinker.rst)
   3. [Build Tinker GPU with CMake](doc/manual/m/install/buildwithcmake.rst)
      or [(Deprecated) Build Tinker GPU with Makefile](doc/manual/m/install/tinkergpu.rst)

Links are only functional on GitHub.


## Repository
<a href="https://github.com/zhi-wang/tinker.gpu">
   <img src="https://github.githubassets.com/images/modules/logos_page/GitHub-Mark.png"
   alt="GitHub Repository" width="100"/>
</a>


For more information, visit the user manual
hosted on [github.io](https://zhi-wang.github.io/tinker.gpu).
