# tinker.gpu directory
src := ../tinker.gpu

# Path to the precompiled libtinker.a
libtinker := $(HOME)/git/tinker/source/libtinker.a

# FFTW
fftwflag := -L$(HOME)/local/fftw/lib -lfftw3 -lfftw3_threads

# CUDA directory
cudadir := /usr/local/cuda

# release / debug / profile
opt := release

# Compilers
CXX := nvcc
NVCC := nvcc

# Fortran runtime
fortrtflag := -L/usr/lib/gcc/x86_64-linux-gnu/5.4.0 -lgfortran -lgomp
