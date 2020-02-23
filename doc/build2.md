# Build Tinker GPU

## Introduction
This is a self-documenting Makefile. After being processed by `doxygen`,
the documentation will appear in the `Related Pages`. Simply run
```
$> make help
```
can also print the documentation in terminal window.

Every option usually has a default value and is documented as
`option=default_value`. In order to override the defaults, set
them in the command line as follows:
```
$> make option1=value1 option2=value2 targets
```

In order to correctly generate this documentations, an empty line is
mandatory at the end of each paragraph.

------

## Maintaining Multiple Build Directories

The design of this Makefile strongly discourages compiling the object files
inside the source code directory, where the situation becomes extreamely
messy if you want to switch back and forth among different configurations.
And as you will see later, quite a few arguments are necessary to be passed
to the `make` file, and the `make` command becomes very long. Therefore, we
recommend you write and keep a shell script in every build directory to
replace the general `make` command, so that a lot of the common options are
saved in the script.

Some examples are given below.

```
#!/bin/bash
# file name: build-v1/z.sh
# GPU/single/release/fftw_dir=$HOME/local/fftw
make ignore=0 host=off prec=s opt=release fftw_dir=$HOME/local/fftw "$@"
```

```
#!/bin/tcsh
# file name: build-v2/z.sh
# CPU/doube/debug/fftw_dir=$HOME/local/fftw
make ignore=0 host=on prec=d opt=debug fftw_dir=$HOME/local/fftw "$argv"
```

------

## Options

### opt=release
Optimization level.
   - Other valid values are `debug` and `profile`.

### prec=real8
Precision of the floating point numbers.
  - These equivalent values all compile the code to double precision:
    `real8`, `8`, `64`, `double`, `d`.
  - These equivalent values all compile the code to single precision:
    `real4`, `4`, `32`, `float`, `single`, `s`.

### host=1
Flag to compile to GPU version or CPU version.
  - These equivalent values all compile the code to GPU version:
    `false`, `off`, `0`.
  - These equivalent values all compile the code to CPU version:
    `true`, `on`, `1`.

### build=build
By default, Makefile can create a `build` directory in the current
working directory. To create a build directory with another name,
change its default value by running

For more information, see `create_build`.

### ignore=1
The default value will minimize the warning messages from the `make`
command, especially when we run `make clean`, `make doc`, etc., where
these warnings and errors are not critical.

`ignore` can also be set to 0, and is recommended to build the excutables.

### tinker_dir (NO DEFAULT)
The directory in which user compiled the `libtinker`.
If this value is not set by a command line argument, the Makefile will
attempt to locate `libtinker` under `$HOME/tinker/source` and emit a fatal
error when it fails to find it. No check will be performed if this option
is set explicitly via command line arguments.

### fftw_dir / fftw_include / fftw_lib (NO DEFAULT)
No default values are set for these three options.

`fftw_dir` is the top-level `FFTW` installation, under which
`include/fftw3.h` and `lib/libfftw3.a` are expected to be found.
If this value is not set by a command line argument, the Makefile will
attempt to locate them under `/usr/local`. However, no fatal error will
be emitted should the Makefile fail to find them.

`fftw_include` and `fftw_lib` will override the values set by `fftw_dir`.
Fatal errors will be emitted if `fftw_x` (`x = include/lib`) is neither set
by `fftw_x` nor `fftw_dir`. No check will be performed on either one of the
explicitly set `fftw_x` in the command line.

### compute_capability=60,70
CUDA Compute Capability
  - 35, 60, 70, 75, etc.
  - Can be comma-separated, e.g. `35,60`.

### cuda_dir=/usr/local/cuda
Top-level `CUDA` installation directory, under which `include` and `lib`
directories are expected to be found.

Sometimes the PGI compiler and the NVCC compiler are not "compatible". For
instance, PGI 19.4 supports CUDA 9.2, 10.0, 10.1, and the default CUDA
version used in PGI 19.4 may be 9.2 and the external NVCC version is 10.1.
One solution is to pass `CUDA_HOME=${cuda_dir}` to the PGI compiler, in
which case, `cuda_dir` should be set to "/usr/local/cuda-10.1".

### fortran_config / cxx_config / acc_config / link_config (UNSPECIFIED)
Addtional files that contain Fortran/C++/OpenAcc/linker related flags.
This is the mechanism we adopt to extend the multi-compiler support.
The default file of each category on different operating systems are
different and there is **NO GUARANTEE** that the defaults will not change.
Tested compilers are tabulated below.

| OS    | Fortran                 | C++                | OpenACC         | Linker          |
|:-----:|:-----------------------:|:------------------:|:---------------:|:---------------:|
| Linux | gfortran 5.4.0          | g++                | pgc++ -ta=telsa | pgc++ -ta=telsa |
| Linux | gfortran 5.4.0          | g++                | g++             | g++             |
| macOS | gfortran 8.3.0 homebrew | clang++ xcode 10.1 | clang++         | clang++         |

------

## Targets

### default
The `tinker.gpu` executable.

### unittest
The `all.tests` executable.

### all
Both `default` and `unittest` targets.

### test
Run the unit tests in a random order.
### dirs
Create sub-directories inside the `$build` directory.
### copy_files (depends on dirs)
Copy files to the `$build` directory.
### create_build (depends on copy_files)
Setup the `$build` directory, including copying the necessary files.
For example, setup a new build directory `build-v3`
inside the `tinker.gpu` directory:

```
$> pwd
#> /home/developer/tinker.gpu
$> make -f make/Makefile create_build build=build-v3
```

### info
Show some of the compiler and linker flags.

### doc
Generate documentation with `doxygen`.

### help
Print the Makefile documentation in the terminal window.

### clean
Clean up the current build directory.

### headers
Test whether the "include" directives in every header file are complete.

