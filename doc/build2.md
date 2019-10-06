# Build Tinker GPU: Makefile, Options, and Targets

## Introduction
This is a self-documenting Makefile. After being processed by `doxygen`,
the documentation will appear in the `Related Pages` section. Simply run
```
$> make help
```
can also print the documentation in terminal window. An empty line is
mandatory at the end of each paragraph of the documentation.

Every option has a default value and is documented as
`option=default_value`. In order to override the default value(s), set
them in the command line as follows:
`$> make option1=value1 option2=value2 targets`

------

## Options

### opt=release
Optimization level. Other valid values are `debug` and `profile`.

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
The Makefile is capable of creating a directory with name `$build` in the
current working directory. This option enables the user to change its
default name.

### tinker_dir (NO DEFAULT)
The directory in which user compiled the `libtinker`.
If this value is not set by a command line argument, the Makefile will
attempt to locate `libtinker` under `$HOME/tinker/source` and emit a fatal
error when it fails to find it. No check will be performed if this option
is set via command line arguments.

### fftw_dir / fftw_include / fftw_lib (NO DEFAULT)
No default values are set for these three options.

`fftw_dir` is the top-level `FFTW` installation, under which
`include/fftw3.h` and `lib/libfftw3` are expected to be found.
If this value is not set by a command line argument, the Makefile will
attempt to locate them under `/usr/local`. However, no fatal error will
be emitted should the Makefile fail to find them.

`fftw_include` and `fftw_lib` will override the values set by `fftw_dir`.
Fatal errors will be emitted if `fftw_x` (`x = include/lib`) is neither set
by `fftw_x` nor `fftw_dir`.

### cuda_dir=/usr/local/cuda
Top-level `CUDA` installation directory, under which `include` and `lib`
directories are expected to be found.

### fortran_config / cxx_config / acc_config / link_config (UNSPECIFIED)
Default configuration files that contain Fortran/C++/OpenAcc/linker
related flags. This is the mechanism we adopt to extend the multi-compiler
support. Tested compilers are tabulated below.

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
For example, setup a new build directory inside the `tinker.gpu` directory
`build-v3`:

```bash
$> pwd
#> /home/developer/tinker.gpu
$> make -f make/Makefile create_build build=build-v3
```

### info
List some of the Makefile information.

### doc
Generate documentation with `doxygen`.

### help
Print the Makefile documentation in the terminal window.

### clean
Clean up the current build directory.

------

## Maintaining Multiple Build Directories

