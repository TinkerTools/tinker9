fortran_compiler_macro__ := -DTINKER_GFORTRAN
fortran_runtime_flags__ := -L/usr/lib/gcc/x86_64-linux-gnu/5.4.0 -lgfortran -lgomp
libtinker_dir__ := $(HOME)/git/tinker/source
shared_libtinker__ := libtinker.so

fftw_dir__ := $(HOME)/local/fftw

cuda_dir__ := /usr/local/cuda

CXX := g++
ifeq ($(opt),debug)
	cxx_opt_flags__ := -O0 -g
else ifeq ($(opt),release)
	cxx_opt_flags__ := -O3
else ifeq ($(opt),profile)
	cxx_opt_flags__ := -O3 -pg
endif
cxx_dylib_flags__ := -shared
cxx_depend_flags__ := -MM

ACC := pgc++
acc_device_compile_flags__ := -acc verystrict -Minfo=accel
ifeq ($(opt),debug)
	acc_opt_flags__ := -O0 -g -traceback
	acc_device_compile_flags__ += -ta=tesla:lineinfo
else ifeq ($(opt),release)
	acc_opt_flags__ := -O3
	acc_device_compile_flags__ += -ta=tesla:fastmath
else ifeq ($(opt),profile)
	acc_opt_flags__ := -O3 -pg
	acc_device_compile_flags__ += -ta=tesla:fastmath,lineinfo
endif
acc_device_link_flags__ := -acc -Mcuda -Mcudalib=cufft
# acc_host_compile_flags__ := -acc verystrict -Minfo=accel -ta=multicore
# acc_host_link_flags__ := -ta=multicore
acc_host_compile_flags__ :=
acc_host_link_flags__ :=
acc_link_flags__ := -rdynamic -Wl,-rpath='$$ORIGIN'
acc_depend_flags__ := -MM -P
