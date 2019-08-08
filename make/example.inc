# **********************************************************
# fortran compiler used to compiler libtinker
# set to TINKER_GFORTRAN for gnu fortran compiler
# set to TINKER_IFORT for intel fortran compiler
# or leave it commented
#
# fortran_compiler_macro__ :=


# **********************************************************
# fortran runtime needed by libtinker
# e.g.
# mac os, homebrew, gfortran 8.3.0
# -L/usr/local/Cellar/gcc/8.3.0/lib/gcc/8 -lgfortran -lgomp
# mac os, gfortran 8.1.0
# -L/usr/local/lib -lgfortran -lgomp
# ubuntu 16.04.6, apt install, gforran 5.4.0
# -L/usr/lib/gcc/x86_64-linux-gnu/5.4.0 -lgfortran -lgomp
#
# fortran_runtime_flags__ :=


# **********************************************************
# libtinker
# path to the precompiled shared tinker library, and
# file name of the the precompiled shared tinker library,
# e.g. libtinker.so or libtinker.dylib
#
# libtinker_dir__ :=
# shared_libtinker__ :=


# **********************************************************
# top-level fftw directory and fftw flags
#
# fftw_dir__ :=


# **********************************************************
# top-level cuda directory
#
# cuda_dir__ := /usr/local/cuda


# **********************************************************
# c++ compiler and related optimization flags
#
# CXX :=
# ifeq ($(opt),debug)
# 	cxx_opt_flags__ :=
# else ifeq ($(opt),release)
# 	cxx_opt_flags__ :=
# else ifeq ($(opt),profile)
# 	cxx_opt_flags__ :=
# endif
# cxx_dylib_flags__ :=
# cxx_depend_flags__ :=


# **********************************************************
# OpenACC compiler and related flags
#
# ACC :=
# acc_opt_flags__ :=
# acc_device_compile_flags__ :=
# ifeq ($(opt),debug)
# 	acc_opt_flags__ +=
# 	acc_device_compile_flags__ +=
# else ifeq ($(opt),release)
# 	acc_opt_flags__ +=
# 	acc_device_compile_flags__ +=
# else ifeq ($(opt),profile)
# 	acc_opt_flags__ +=
# 	acc_device_compile_flags__ +=
# endif
# acc_device_link_flags__ :=
# acc_host_compile_flags__ :=
# acc_host_link_flags__ :=
# acc_link_flags__ :=
# acc_depend_flags__ :=
