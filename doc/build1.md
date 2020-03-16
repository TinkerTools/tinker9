# Build the Canonical Tinker


## Checkout the Required Version of Tinker
Clone Tinker from GitHub, checkout the required version (see `Prerequisites`),
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


## Make libtinker
   - Change the value of variable `FFTWDIR` to the top-level `FFTW` installation.
   - Locate the correct flags in the `Makefile` for your compiler and operating system,
     uncomment them, and comment out the others.
   - **Must disable the OpenMP flag.**
   - **Change optimization level to -O2.**
   - Run `make` command
```
make libtinker.a -j
```
