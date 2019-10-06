# Build the Canonical Tinker

* Clone Tinker from GitHub, checkout the required version (see section `Prerequisites`),
and copy the `Makefile` to the `source` directory, e.g.,
```
$ pwd
> /home/user
$ git clone https://github.com/tinkertools/tinker
$ cd tinker
$ checkout <VERSION>
$ cd source
$ cp ../make/Makefile .
$ pwd
> /home/user/tinker/source
```
* Modify the `Makefile`
    * change the value of variable `FFTWDIR` to the top-level `FFTW` installation;
    * locate the correct flags in the `Makefile` for your compiler and operating system,
    uncomment them, and comment out others;
    * modify the `OPTFLAGS`: use a more conservative optimization flag `-O2` instead of
    `-Ofast` or `-O3`; add `-fPIC`; keep the `OpenMP` flag;
    * add a new target to compile a shared library:
```
# make sure TABs are used for indentations in the Makefile
libtinker.so: $(OBJS)
        ${F77} -fPIC -shared -o $@ \
	<various.o files, refer to libtinker.a>
```
* Run command
```
make libtinker.so
```

<a name='libtinkergpu'></a>
# Build Tinker GPU Library

* Clone Tinker GPU from GitHub, and create a build directory
with any name anywhere.
Preferably, create `build` inside `tinker.gpu`.
```
$ pwd
> /home/user
$ git clone <GitHub/tinker.gpu>
$ cd tinker.gpu
$ mkdir build
$ cd build
```
* Copy the `Makefile` and `*.inc` files to the build directory.
```
$ pwd
> /home/user/tinker.gpu/build
$ cp ../make/Makefile ../make/*.inc .
```
* If you want compile the **CPU version** and **if you are lucky**, run:
```
make
```

Don't worry if you feel unlucky.
There are following orthogonal options avaiable and you can run
```
make <opt1=value1> <opt2=value2> <more options and values>
```
to overwrite the default values.
* `src=..` by default
    * path to the directory where you can see the Tinker GPU source directories,
    e.g. `include`, `src`, etc.
* `opt=release` by default
    * optimization options;
    * other values: `debug`, `profile`.
* `prec=real8` by default using 64-bit floating numbers
    * precisions of the floating numbers, either using 64-bit or 32-bit;
    * other values: `real4` using 32-bit floating numbers;
    `8`, `64`, `double`, `d` are equivalent to `real8`;
    `4`, `32`, `float`, `single`, `s` are equivalent to `real4`.
* `host=true` by default to compile to CPU version
    * flag to compile to CPU version or GPU version;
    * other values: `false`, to compile to GPU version;
    `on`, `1` are equivalent to `true`;
    `off`, `0` are equivalent to `false`.
* `config=<configuration file>`
    * OS and compiler specific configuration file;
    * On Linux system, `Makefile` will look for `linux.gfortran.inc` file by default.

It's very possible you need to modify the first several lines
of the configuration file as well.
Most likely, you will specify the Fortran compiler you used to compile `libtinker`,
the build directory of `libtinker`, the directory of `FFTW` library installation,
and most importantly, Fortran run-time libraries with the directory of Fortran
run-time installation.

They are documented inside `make/example.inc`.

* If you want to compile a debug GPU version with single precision, you can run
```
make opt=debug prec=s host=0
```
