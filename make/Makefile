############################################################
#                                                          #
#  Copy Makefile and inc files to the build directory      #
#                                                          #
############################################################

# options following the `make' command
# opt :=
# src :=
# prec :=
# hostonly :=
# config :=

############################################################

# Darwin -> darwin; Linux -> linux
os__ := $(shell uname -s 2>/dev/null | tr "[:upper:]" "[:lower:]")
ifeq ($(os__),linux)
	shared_lib_suffix__ := so
else ifeq ($(os__),darwin)
	shared_lib_suffix__ := dylib
endif

# release / debug / profile
opt := release

# tinker.gpu directory
src := ..

# GPU Precision
# real4 / 4 / 32 / float / single / s -> real4__
# real8 / 8 / 64 / double         / d -> real8__
prec := real4
ifeq ($(prec),$(filter $(prec),real4 4 32 float single s))
	prec__ := real4__
	precision_macro__ := -DTINKER_GPU_SINGLE
else ifeq ($(prec),$(filter $(prec),real8 8 64 double d))
	prec__ := real8__
	precision_macro__ := -DTINKER_GPU_DOUBLE
endif

hostonly := false

config := default_config__
ifeq ($(config),default_config__)
	ifeq ($(os__),linux)
		config := linux.gfortran.inc
	endif
endif
include $(config)
ifeq ($(prec__),real4__)
	acc_host_link_flags__ += -L$(fftw_dir__)/lib -lfftw3f
endif
fftw_flags__ := -L$(fftw_dir__)/lib -lfftw3 -lfftw3_threads

############################################################
#                                                          #
# default target: all                                      #
#                                                          #
############################################################

tinker_gpu_exec__ := tinker.gpu
all_tests_exec__ := all.tests
libtinkergpu_stem__ := tinkergpu
src_obj_dir__ := src
ifeq ($(prec__),real4__)
	precision_suffix__ :=
else ifeq ($(prec__),real8__)
	# -dp: double precision
	precision_suffix__ := -dp
	tinker_gpu_exec__ := $(tinker_gpu_exec__)$(precision_suffix__)
	all_tests_exec__ := $(all_tests_exec__)$(precision_suffix__)
	libtinkergpu_stem__ := $(libtinkergpu_stem__)$(precision_suffix__)
	src_obj_dir__ := $(src_obj_dir__)$(precision_suffix__)
endif
libtinkergpu_name__ := lib$(libtinkergpu_stem__)

.PHONY: all
all: $(tinker_gpu_exec__)

############################################################
#                                                          #
# some utilities                                           #
#                                                          #
############################################################

.PHONY: print-config
print-config:
	@echo 'config file       =' $(config)
	@echo 'opt level         =' $(opt)
	@echo 'source directory  =' $(top_dir__)
	@echo 'precision         =' $(prec)
	@echo 'host-only         =' $(hostonly)
	@echo 'OS                =' $(os__)
	@echo 'shared library    =' $(libtinkergpu_name__).$(shared_lib_suffix__)
	@echo 'executables       =' $(tinker_gpu_exec__) $(all_tests_exec__)
	@echo 'cxx compile flags =' $(cxx_compile_flags__)
	@echo 'acc compile flags =' $(acc_compile_flags__)
	@echo 'acc link flags    =' $(acc_link_flags__)

.PHONY: dirs
dirs:
	@mkdir -p $(ext_dir__)/fmt
	@mkdir -p $(ext_dir__)/tinker
	@mkdir -p $(src_obj_dir__)/test
	@mkdir -p $(src_obj_dir__)/gpu

.PHONY: clean
clean:
	rm -rf ext src* tinker.gpu* all.tests* libtinker* *dSYM

########################################################################
#                                                                      #
# libraries                                                            #
#                                                                      #
########################################################################

top_dir__ := $(realpath $(src))
ext_dir__ := ext/ext

include_flags__ := -I$(top_dir__)/include
include_flags__ += -I$(top_dir__)/ext
include_flags__ += -I$(top_dir__)/ext/ext
ifeq ($(hostonly),false)
	include_flags__ += -I$(cuda_dir__)/include
else
	include_flags__ += -I$(fftw_dir__)/include -DTINKER_HOSTONLY
endif

common_flags__ := -D$(fortran_compiler_macro__) \
$(precision_macro__) $(include_flags__)
cxx_common_flags__ := $(common_flags__) $(cxx_opt_flags__)
acc_common_flags__ := $(common_flags__) $(acc_opt_flags__)

cxx_compile_flags__ := $(cxx_common_flags__) -std=c++11 -c -fpic
acc_compile_flags__ := $(acc_common_flags__) -std=c++11 -c -fpic
ifeq ($(hostonly),false)
	acc_compile_flags__ += $(acc_device_compile_flags__)
else
	acc_compile_flags__ += $(acc_host_compile_flags__)
endif

acc_link_flags__ := -rdynamic -Wl,-rpath='$$ORIGIN'
ifeq ($(hostonly),false)
	acc_link_flags__ += $(acc_device_link_flags__)
else
	acc_link_flags__ += $(acc_host_link_flags__)
endif
acc_link_flags__ += $(acc_common_flags__) -L. -ltinker \
-l$(libtinkergpu_stem__) $(fortran_runtime_flags__) $(fftw_flags__)

########################################################################

ext_cc_objs__ := $(ext_dir__)/fmt/format.o $(ext_dir__)/fmt/posix.o
ext_cpp_objs__ := $(ext_dir__)/tinker/tinker_mod.o \
$(ext_dir__)/tinker/tinker_rt.o

main_cc_objs__ := $(src_obj_dir__)/test/main_test_all.o \
$(src_obj_dir__)/main_tinker_gpu.o

# $(top_dir__)/src/foo.cpp -> src-dp/foo.o
lib_cpp_files__ := $(shell find $(top_dir__)/src -type f -name '*.cpp' \
-not -path '$(top_dir__)/$(ext_dir__)/*' \
-not -path '$(top_dir__)/src/test/*')
lib_cpp_objs__ := $(subst $(top_dir__)/src,$(src_obj_dir__),$(lib_cpp_files__))
lib_cpp_objs__ := $(patsubst %.cpp,%.o,$(lib_cpp_objs__))

# $(top_dir__)/src/gpu/foo.cc -> src-dp/gpu/foo.o
lib_acc_files__ := $(wildcard $(top_dir__)/src/gpu/*.cc)
lib_acc_objs__ := $(subst $(top_dir__)/src,$(src_obj_dir__),$(lib_acc_files__))
lib_acc_objs__ := $(patsubst %.cc,%.o,$(lib_acc_objs__))

# $(top_dir__)/src/test/foo.cpp -> src-dp/test/foo.o
test_cpp_files__ := $(wildcard $(top_dir__)/src/test/*.cpp)
test_cpp_objs__ := $(subst $(top_dir__)/src,$(src_obj_dir__),$(test_cpp_files__))
test_cpp_objs__ := $(patsubst %.cpp,%.o,$(test_cpp_objs__))

all_cc_objs__ := $(main_cc_objs__) $(ext_cc_objs__)
all_cpp_objs__ := $(ext_cpp_objs__) $(lib_cpp_objs__) $(test_cpp_objs__)

# targets: normal-prerequisites | order-only-prerequisites
$(all_cc_objs__): %.o: $(top_dir__)/%.cc | dirs
	$(CXX) $(cxx_compile_flags__) $< -o $@
$(all_cpp_objs__): %.o: $(top_dir__)/%.cpp | dirs
	$(CXX) $(cxx_compile_flags__) $< -o $@
$(lib_acc_objs__): %.o: $(top_dir__)/%.cc | dirs
	$(ACC) $(acc_compile_flags__) $< -o $@

########################################################################

# shared libtinkergpu does not include lib_acc_objs__
$(libtinkergpu_name__).$(shared_lib_suffix__): \
$(ext_cc_objs__) $(ext_cpp_objs__) $(lib_cpp_objs__)
	$(CXX) -shared -fpic -o $@ $^

# shared libtinker
$(shared_libtinker__): $(libtinker_dir__)/$(shared_libtinker__)
	cp $^ $@


########################################################################
#                                                                      #
# tinker.gpu                                                           #
#                                                                      #
########################################################################

$(tinker_gpu_exec__): $(src_obj_dir__)/main_tinker_gpu.o \
$(lib_acc_objs__) \
| $(libtinkergpu_name__).$(shared_lib_suffix__) $(shared_libtinker__)
	$(ACC) $^ $(acc_link_flags__) -o $@


########################################################################
#                                                                      #
# tests                                                                #
#                                                                      #
########################################################################

.PHONY: test
test: $(all_tests_exec__)
	$(PWD)/$(all_tests_exec__) [ff] --durations yes

.PHONY: unittest
unittest: $(all_tests_exec__)

$(all_tests_exec__): $(src_obj_dir__)/test/main_test_all.o \
$(test_cpp_objs__) $(lib_acc_objs__) \
| $(libtinkergpu_name__).$(shared_lib_suffix__) $(shared_libtinker__)
	$(ACC) $^ $(acc_link_flags__) -o $@
