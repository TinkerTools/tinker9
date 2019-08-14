# Build and Install

# Table of Contents
* [Build Tinker Library](#libtinker)
* [Build Tinker GPU Library](#libtinkergpu)
* [Install](#install)

<a name='libtinker'></a>
# Build Tinker Library
* clone Tinker from GitHub, checkout the required version (see section `Prerequisites`),
and copy the `Makefile` to the `source` directory, e.g.,
```
$ pwd
> /home/user
$ git clone https://github.com/TinkerTools/Tinker
$ cd Tinker
$ checkout <VERSION>
$ cd source
$ cp ../make/Makefile .
$ pwd
> /home/user/Tinker/source
```
* modify the `Makefile`
    * change the value of variable `FFTWDIR` to the top-level `FFTW` installation;
    * locate the correct flags in the `Makefile` for your compiler and operating system,
    uncomment them, and comment out others;
    * modify the `OPTFLAGS`: use a more conservative optimization flag `-O2` instead of
    `-Ofast` or `-O3`; add `-fPIC`; keep the `OpenMP` flag;
    * add a new target to compile a shared library:
```
# make sure TABs are used for indentations in the Makefile
# e.g. on Linux with gfortran
libtinker.so: $(OBJS)
	gfortran -shared -fPIC -o $@ \
	<various.o files, refer to libtinker.a>
```
* run command
```
make libtinker.so
```

<a name='libtinkergpu'></a>
# Build Tinker GPU Library


<a name='install'></a>
# Install
