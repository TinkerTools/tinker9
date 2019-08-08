############################################################
#                                                          #
#  Copy Makefile and inc files to the build directory      #
#                                                          #
############################################################

# options following the `make' command
# src :=
# opt :=
# prec :=
# host :=
# config :=

############################################################

# tinker.gpu directory
src := ..
top_dir__ := $(realpath $(src))
include $(top_dir__)/make/name.make


# default optimization
# var opt; value release / debug / profile
opt := release
ifeq ($(opt),debug)
	debug_macro__ := -DDEBUG=1
else ifeq ($(opt),release)
	debug_macro__ := -DNDEBUG
else ifeq ($(opt),profile)
	debug_macro__ := -DNDEBUG
endif


# var os__; value Darwin -> darwin / Linux -> linux
# var shared_lib_suffix__; value so / dylib
include $(top_dir__)/make/os.make


# default precision
# var prec
# value
# real4 / 4 / 32 / float / single / s -> real4__
# real8 / 8 / 64 / double         / d -> real8__
prec := real8
# var prec__; value real4__ / real8__
# var precision_macro__; value -D...
include $(top_dir__)/make/prec.make


# default host
# var host; value true on 1 / false off 0
host := true
# var host_only__; value true / false
# var host_only_macro__
include $(top_dir__)/make/host_only.make


config := default_config__
ifeq ($(config),default_config__)
	ifeq ($(os__),linux)
		config := linux.gfortran.inc
	else ifeq ($(os__),darwin)
		config := macos.clang.gfortran.inc
	endif
endif
include $(config)
fftw_flags__ := -L$(fftw_dir__)/lib -lfftw3 -lfftw3_threads


############################################################
#                                                          #
# default target and all                                   #
#                                                          #
############################################################


.PHONY: default_target__
default_target__: $(tinker_gpu_exec__)


.PHONY: all
all: $(tinker_gpu_exec__) unittest


.PHONY: unittest
unittest: $(all_tests_exec__)


.PHONY: test
test: $(all_tests_exec__)
	$(PWD)/$(all_tests_exec__) [ff] --durations yes


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


.PHONY: doc
doc:
	doxygen $(top_dir__)/doc/Doxyfile


.PHONY: clean
clean:
	rm -rf ext src* tinker.gpu* all.tests* libtinker* *dSYM


########################################################################
#                                                                      #
# flags                                                                #
#                                                                      #
########################################################################


# macros
all_macros__ := $(fortran_compiler_macro__) $(precision_macro__) \
$(host_only_macro__) $(precision_macro__)

# include directories
include_flags__ := -I$(top_dir__)/include
include_flags__ += -I$(top_dir__)/ext
include_flags__ += -I$(top_dir__)/ext/ext
ifeq ($(host_only__),false)
	include_flags__ += -I$(cuda_dir__)/include
else
	include_flags__ += -I$(fftw_dir__)/include
endif

common_flags__ := $(all_macros__) $(include_flags__) -std=c++11
acc_common_flags__ := $(common_flags__) $(acc_opt_flags__)
cxx_common_flags__ := $(common_flags__) $(cxx_opt_flags__)

# compile flags
acc_compile_flags__ := $(acc_common_flags__) -c -fpic
ifeq ($(host_only__),false)
	acc_compile_flags__ += $(acc_device_compile_flags__)
else
	acc_compile_flags__ += $(acc_host_compile_flags__)
endif
cxx_compile_flags__ := $(cxx_common_flags__) -c -fpic

# acc link flags
ifeq ($(prec__),real4__)
	fftw_flags__ += -lfftw3f
endif
ifeq ($(host_only__),false)
	acc_link_flags__ += $(acc_device_link_flags__)
else
	acc_link_flags__ += $(acc_host_link_flags__)
endif
acc_link_flags__ += $(acc_opt_flags__) -L. -ltinker \
-l$(libtinkergpu_stem__) $(fortran_runtime_flags__) $(fftw_flags__)


########################################################################
#                                                                      #
# dependencies, object files, and libraries                            #
#                                                                      #
########################################################################


# targets: normal-prerequisites | order-only-prerequisites
$(acc_dependency__): %.d: $(top_dir__)/%.cpp | dirs
	printf %s $(@D)/ > $@
	$(ACC) $(acc_common_flags__) $< $(acc_depend_flags__) >> $@
$(cc_dependency__): %.d: $(top_dir__)/%.cc | dirs
	printf %s $(@D)/ > $@
	$(CXX) $(cxx_common_flags__) $< $(cxx_depend_flags__) >> $@
$(cpp_dependency__): %.d: $(top_dir__)/%.cpp | dirs
	printf %s $(@D)/ > $@
	$(CXX) $(cxx_common_flags__) $< $(cxx_depend_flags__) >> $@
-include $(cc_dependency__)
-include $(cpp_dependency__)
-include $(acc_dependency__)


$(all_acc_objs__): %.o: $(top_dir__)/%.cpp | dirs
	$(ACC) $(acc_compile_flags__) $< -o $@
$(all_cc_objs__): %.o: $(top_dir__)/%.cc | dirs
	$(CXX) $(cxx_compile_flags__) $< -o $@
$(all_cpp_objs__): %.o: $(top_dir__)/%.cpp | dirs
	$(CXX) $(cxx_compile_flags__) $< -o $@


# shared tinkergpu library does not include lib_acc_objs__
# lib_acc_objs__ should be statically linked
$(libtinkergpu_name__).$(shared_lib_suffix__): $(ext_cc_objs__) $(lib_cpp_objs__)
	$(CXX) $(cxx_dylib_flags__) -o $@ $^
# copy the shared tinker library to the build directory
$(shared_libtinker__): $(libtinker_dir__)/$(shared_libtinker__)
	cp $^ $@


########################################################################
#                                                                      #
# executables                                                          #
#                                                                      #
########################################################################


# tinker.gpu: main.o acc_objs.o | libtinkergpu.so libtinker.so
# 	link tinker.gpu
$(tinker_gpu_exec__): $(src_obj_dir__)/main_tinker_gpu.o $(lib_acc_objs__) \
| $(libtinkergpu_name__).$(shared_lib_suffix__) $(shared_libtinker__)
	$(ACC) $^ $(acc_link_flags__) -o $@


# all.tests: main.o test_objs.o acc_objs.o | libtinkergpu.so libtinker.so
# 	link all.tests
$(all_tests_exec__): $(src_obj_dir__)/test/main_all_tests.o \
$(test_cpp_objs__) $(lib_acc_objs__) \
| $(libtinkergpu_name__).$(shared_lib_suffix__) $(shared_libtinker__)
	$(ACC) $^ $(acc_link_flags__) -o $@



