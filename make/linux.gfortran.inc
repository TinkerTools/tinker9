fortran_compiler_macro__ := -DTINKER_GFORTRAN
fortran_runtime_flags__ := -L/usr/lib/gcc/x86_64-linux-gnu/5.4.0 -lgfortran -lgomp

ifdef libtinker_dir
libtinker_dir__ = $(libtinker_dir)
else
$(error libtinker_dir is not set)
endif

shared_libtinker__ := libtinker.so

cuda_dir__ := /home/levnon/.local/opt/pgi/linux86-64-llvm/2019/cuda/10.0

ifdef ACC
ACC := $(ACC)
else
ACC := g++
endif

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

acc_device_compile_flags__ := -acc verystrict -Minfo=accel

ifeq ($(opt),debug)
	acc_opt_flags__ := -O0 -g -traceback
	acc_device_compile_flags__ += -ta=tesla
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
