############################################################
#                                                          #
#  Copy Makefile and inc files to the build directory      #
#                                                          #
############################################################

# options following the `make' command
# opt :=
# src :=
# prec :=
# host :=
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
prec := real8
ifeq ($(prec),$(filter $(prec),real4 4 32 float single s))
	prec__ := real4__
	precision_macro__ := -DTINKER_SINGLE_PRECISION
else ifeq ($(prec),$(filter $(prec),real8 8 64 double d))
	prec__ := real8__
	precision_macro__ := -DTINKER_DOUBLE_PRECISION
endif

host := true

config := default_config__
ifeq ($(config),default_config__)
	ifeq ($(os__),linux)
		config := linux.gfortran.inc
	else ifeq ($(os__),darwin)
		config := macos.gfortran.inc
	endif
endif
include $(config)
ifeq ($(prec__),real4__)
	acc_host_link_flags__ += -L$(fftw_dir__)/lib -lfftw3f
endif
fftw_flags__ := -L$(fftw_dir__)/lib -lfftw3 -lfftw3_threads

############################################################
#                                                          #
# default target and all                                   #
#                                                          #
############################################################

tinker_gpu_exec__ := tinker.gpu
all_tests_exec__ := all.tests
libtinkergpu_stem__ := tinkergpu
src_obj_dir__ := src
libtinkergpu_name__ := lib$(libtinkergpu_stem__)

.PHONY: default_target__
default_target__: $(tinker_gpu_exec__)

.PHONY: all
all: $(tinker_gpu_exec__) unittest


############################################################
#                                                          #
# some utilities                                           #
#                                                          #
############################################################

host_only__ :=
ifeq ($(host),$(filter $(host),false off 0))
	host_only__ := false
else ifeq ($(host),$(filter $(host),true on 1))
	host_only__ := true
endif

.PHONY: print-config
print-config:
	@echo 'config file       =' $(config)
	@echo 'opt level         =' $(opt)
	@echo 'source directory  =' $(top_dir__)
	@echo 'precision         =' $(prec)
	@echo 'host only         =' $(host_only__)
	@echo 'OS                =' $(os__)
	@echo 'shared library    =' $(libtinkergpu_name__).$(shared_lib_suffix__)
	@echo 'executables       =' $(tinker_gpu_exec__) $(all_tests_exec__)
	@echo 'cxx compile flags =' $(cxx_compile_flags__)
	@echo 'acc compile flags =' $(acc_compile_flags__)
	@echo 'acc link flags    =' $(acc_link_flags__)

.PHONY: dirs
dirs:
	@mkdir -p $(ext_dir__)/fmt
	@mkdir -p $(src_obj_dir__)/test

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
ifeq ($(host_only__),false)
	include_flags__ += -I$(cuda_dir__)/include
else
	include_flags__ += -I$(fftw_dir__)/include -DTINKER_HOST
endif

common_flags__ := -D$(fortran_compiler_macro__) \
$(precision_macro__) $(include_flags__) -std=c++11
cxx_common_flags__ := $(common_flags__) $(cxx_opt_flags__)
acc_common_flags__ := $(common_flags__) $(acc_opt_flags__)

cxx_compile_flags__ := $(cxx_common_flags__) -c -fpic
acc_compile_flags__ := $(acc_common_flags__) -c -fpic
ifeq ($(host_only__),false)
	acc_compile_flags__ += $(acc_device_compile_flags__)
else
	acc_compile_flags__ += $(acc_host_compile_flags__)
endif

ifeq ($(host_only__),false)
	acc_link_flags__ += $(acc_device_link_flags__)
else
	acc_link_flags__ += $(acc_host_link_flags__)
endif
acc_link_flags__ += $(acc_common_flags__) -L. -ltinker \
-l$(libtinkergpu_stem__) $(fortran_runtime_flags__) $(fftw_flags__)

########################################################################

ext_cc_objs__ := $(ext_dir__)/fmt/format.o $(ext_dir__)/fmt/posix.o

main_cpp_files__ := $(top_dir__)/src/main_tinker_gpu.cpp \
$(top_dir__)/src/test/main_all_tests.cpp
main_cpp_objs__ := $(subst $(top_dir__)/src,$(src_obj_dir__),$(main_cpp_files__))
main_cpp_objs__ := $(patsubst %.cpp,%.o,$(main_cpp_objs__))

# $(top_dir__)/src/acc_foo.cpp -> src/acc_foo.o
lib_acc_files__ := $(wildcard $(top_dir__)/src/acc_*.cpp)
lib_acc_objs__ := $(subst $(top_dir__)/src,$(src_obj_dir__),$(lib_acc_files__))
lib_acc_objs__ := $(patsubst %.cpp,%.o,$(lib_acc_objs__))

# $(top_dir__)/src/test/foo.cpp -> src/test/foo.o
test_cpp_files__ := $(shell find $(top_dir__)/src/test -type f -name '*.cpp' \
-not -path '$(top_dir__)/src/test/main_*.cpp')
test_cpp_objs__ := $(subst $(top_dir__)/src,$(src_obj_dir__),$(test_cpp_files__))
test_cpp_objs__ := $(patsubst %.cpp,%.o,$(test_cpp_objs__))

# $(top_dir__)/src/foo.cpp -> src/foo.o
lib_cpp_files__ := $(shell find $(top_dir__)/src -type f -name '*.cpp' \
-not -path '$(top_dir__)/src/main_*.cpp' \
-not -path '$(top_dir__)/src/acc_*.cpp' \
-not -path '$(top_dir__)/src/test/*')
lib_cpp_objs__ := $(subst $(top_dir__)/src,$(src_obj_dir__),$(lib_cpp_files__))
lib_cpp_objs__ := $(patsubst %.cpp,%.o,$(lib_cpp_objs__))

all_cc_objs__ := $(ext_cc_objs__)
all_cpp_objs__ := $(main_cpp_objs__) $(lib_cpp_objs__) $(test_cpp_objs__)

# targets: normal-prerequisites | order-only-prerequisites
$(all_cc_objs__): %.o: $(top_dir__)/%.cc | dirs
	$(CXX) $(cxx_compile_flags__) $< -o $@
$(all_cpp_objs__): %.o: $(top_dir__)/%.cpp | dirs
	$(CXX) $(cxx_compile_flags__) $< -o $@
$(lib_acc_objs__): %.o: $(top_dir__)/%.cpp | dirs
	$(ACC) $(acc_compile_flags__) $< -o $@

########################################################################

# shared tinkergpu library does not include lib_acc_objs__
# lib_acc_objs__ should be statically linked
$(libtinkergpu_name__).$(shared_lib_suffix__): $(ext_cc_objs__) $(lib_cpp_objs__)
	$(CXX) $(cxx_dylib_flags__) -o $@ $^

# copy the shared tinker library to the build directory
$(shared_libtinker__): $(libtinker_dir__)/$(shared_libtinker__)
	cp $^ $@


########################################################################
#                                                                      #
# tinker.gpu                                                           #
#                                                                      #
########################################################################

# tinker.gpu: main.o acc_objs.o | libtinkergpu.so libtinker.so
# 	link tinker.gpu
$(tinker_gpu_exec__): $(src_obj_dir__)/main_tinker_gpu.o $(lib_acc_objs__) \
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

# all.tests: main.o test_objs.o acc_objs.o | libtinkergpu.so libtinker.so
# 	link all.tests
$(all_tests_exec__): $(src_obj_dir__)/test/main_all_tests.o \
$(test_cpp_objs__) $(lib_acc_objs__) \
| $(libtinkergpu_name__).$(shared_lib_suffix__) $(shared_libtinker__)
	$(ACC) $^ $(acc_link_flags__) -o $@


########################################################################
#                                                                      #
# doc                                                                  #
#                                                                      #
########################################################################

.PHONY: doc
doc:
	doxygen $(top_dir__)/doc/Doxyfile


########################################################################
#                                                                      #
# dependency                                                           #
#                                                                      #
########################################################################

cc_dependency__ := $(all_cc_objs__:%.o=%.d)
cpp_dependency__ := $(all_cpp_objs__:%.o=%.d)
acc_dependency__ := $(lib_acc_objs__:%.o=%.d)

$(cc_dependency__): %.d: $(top_dir__)/%.cc | dirs
	printf %s $(@D)/ > $@
	$(CXX) $(common_flags__) $(cxx_opt_flags__) $< $(cxx_depend_flags__) >> $@
$(cpp_dependency__): %.d: $(top_dir__)/%.cpp | dirs
	printf %s $(@D)/ > $@
	$(CXX) $(common_flags__) $(cxx_opt_flags__) $< $(cxx_depend_flags__) >> $@
$(acc_dependency__): %.d: $(top_dir__)/%.cpp | dirs
	printf %s $(@D)/ > $@
	$(ACC) $(common_flags__) $(acc_opt_flags__) $< $(acc_depend_flags__) >> $@

-include $(cc_dependency__)
-include $(cpp_dependency__)
-include $(acc_dependency__)
