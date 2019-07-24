fortran_compiler_macro__ := TINKER_GFORTRAN
fortran_runtime_flags__ := -L/usr/lib/gcc/x86_64-linux-gnu/5.4.0 -lgfortran -lgomp
libtinker_dir__ := $(HOME)/git/tinker/source
shared_libtinker__ := libtinker.so

fftw_dir__ := $(HOME)/local/fftw

cuda_dir__ := /usr/local/cuda

CXX := g++
ifeq ($(opt),debug)
	cxx_opt_flags__ := -O0 -DDEBUG=1 -g
else ifeq ($(opt),release)
	cxx_opt_flags__ := -O3 -DNDEBUG
else ifeq ($(opt),profile)
	cxx_opt_flags__ := -O3 -DNDEBUG -pg
endif

ACC := pgc++
acc_opt_flags__ :=
acc_device_compile_flags__ := -acc verystrict -Minfo=accel
ifeq ($(opt),debug)
	acc_opt_flags__ += -O0 -DDEBUG=1 -g -traceback
	acc_device_compile_flags__ += -ta=tesla:lineinfo
else ifeq ($(opt),release)
	acc_opt_flags__ += -O3 -DNDEBUG
	acc_device_compile_flags__ += -ta=tesla:fastmath
else ifeq ($(opt),profile)
	acc_opt_flags__ += -O3 -DNDEBUG -pg
	acc_device_compile_flags__ += -ta=tesla:fastmath,lineinfo
endif
acc_device_link_flags__ := -acc -Mcuda -Mcudalib=cufft
acc_host_compile_flags__ := -acc verystrict -Minfo=accel -ta=multicore
acc_host_link_flags__ := -ta=multicore
