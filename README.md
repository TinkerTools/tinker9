Tinker9: Tools for Molecular Design with GPU Support
====================================================
[//]: # (Badges)
[![Build Status](https://travis-ci.com/tinkertools/tinker9.svg?branch=master)](https://travis-ci.com/tinkertools/tinker9)
[![Docs Status](https://readthedocs.org/projects/tinker9-manual/badge/?version=latest&style=flat)](https://tinker9-manual.readthedocs.io)

## Introduction
Tinker9 is a complete rewrite and new version of the canonical Tinker software, currently Tinker8. Tinker9 consists of C++ code with
OpenACC directives and CUDA kernels allow efficient usage on GPUs. At present, Tinker9 builds against the object library from Tinker8,
and provides GPU versions of the Tinker DYNAMIC, MINIMIZE, ANALYZE and TESTGRAD programs. Over time we plan to port much or all of the
remaining portions of Fortran Tinker8 to the new Tinker9 code base.

## Installation Steps
   1. [Prerequisites](doc/manual/m/install/preq.rst)
   2. [Build the Canonical Tinker (important to get the REQUIRED version)](doc/manual/m/install/tinker.rst)
   3. [Build Tinker9 with CMake](doc/manual/m/install/buildwithcmake.rst)


## User's Manual (In Progress)
The HTML version is hosted on [readthedocs](https://tinker9-manual.readthedocs.io)
and the PDF version is accessible from the same webpage.


## Style Guide
Please read the [style guide](doc/style.md) before contributing to this code base.
