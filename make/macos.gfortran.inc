fortran_compiler_macro__ := TINKER_GFORTRAN
fortran_runtime_flags__ := -L/usr/local/Cellar/gcc/8.3.0/lib/gcc/8 -lgfortran -lgomp
libtinker_dir__ := $(HOME)/git/tinker/source
shared_libtinker__ := libtinker.dylib

fftw_dir__ := $(HOME)/local/fftw

cuda_dir__ :=

CXX := clang++
ifeq ($(opt),debug)
	cxx_opt_flags__ := -O0 -DDEBUG=1 -g
else ifeq ($(opt),release)
	cxx_opt_flags__ := -O3 -DNDEBUG
else ifeq ($(opt),profile)
	cxx_opt_flags__ := -O3 -DNDEBUG
endif
cxx_dylib_flags__ := -dynamiclib -undefined dynamic_lookup
cxx_depend_flags__ := -MM

ACC := clang++
acc_opt_flags__ :=
acc_device_compile_flags__ :=
ifeq ($(opt),debug)
	acc_opt_flags__ += -O0 -DDEBUG=1 -g
	acc_device_compile_flags__ +=
else ifeq ($(opt),release)
	acc_opt_flags__ += -O3 -DNDEBUG
	acc_device_compile_flags__ +=
else ifeq ($(opt),profile)
	acc_opt_flags__ += -O3 -DNDEBUG
	acc_device_compile_flags__ +=
endif
acc_device_link_flags__ :=
# acc_host_compile_flags__ := -acc verystrict -Minfo=accel -ta=multicore
# acc_host_link_flags__ := -ta=multicore
acc_host_compile_flags__ :=
acc_host_link_flags__ :=
acc_link_flags__ :=
acc_depend_flags__ := -MM
