############################################################
#                                                          #
#  Copy Makefile and inc files to the build directory      #
#                                                          #
############################################################

# release / debug / profile
opt := release

# tinker.gpu directory
src := ..

# Path to the precompiled libtinker.a
libtinker := $(HOME)/git/tinker/source/libtinker.a

# FFTW
fftwflag := -L$(HOME)/local/fftw/lib -lfftw3 -lfftw3_threads

# Precision -DTINKER_GPU_DOUBLE or leave it empty
precision_flag ?=
#precision_flag ?= -DTINKER_GPU_DOUBLE

# Darwin -> darwin
os := $(shell uname -s 2>/dev/null | tr "[:upper:]" "[:lower:]")

# CUDA directory
cudadir := /usr/local/cuda

# Compilers
CXX := nvcc
NVCC := nvcc
ACC := pgc++

config := default_config
ifeq ($(config),default_config)
	ifeq ($(os),linux)
		config := linux.gfortran.inc
	else ifeq ($(os),darwin)
		config := macos.gfortran.inc
	endif
endif
include $(config)

############################################################

ifeq ($(os),linux)
	ifeq ($(opt),debug)
		cxx_opt_flags := -O0 -g -DDEBUG=1
	else ifeq ($(opt),release)
		cxx_opt_flags := -O3 -DNDEBUG
	else ifeq ($(opt),profile)
		cxx_opt_flags := -O3 -DNDEBUG -g
	endif
	ar_flags := ar crusv
	ranlib_flags := ranlib
	os_link_flags := -Xlinker --export-dynamic
else ifeq ($(os),darwin)
	ifeq ($(opt),debug)
                cxx_opt_flags := -O0 -g -DDEBUG=1
        else ifeq ($(opt),release)
                cxx_opt_flags := -O3 -DNDEBUG
        else ifeq ($(opt),profile)
                cxx_opt_flags := -O3 -DNDEBUG -g
        endif
	ar_flags := ar crusv
	ranlib_flags := ranlib -c
	os_link_flags :=
endif

############################################################
#                                                          #
# default target: all                                      #
#                                                          #
############################################################

tinker_gpu_name := tinker.gpu
all_tests_name := all.tests
libtinkergpu_name := libtinkergpu
src_objdir := src
ifeq ($(precision_flag),)
	precision_suffix :=
else
	# -dp: double precision
	precision_suffix := -dp
	tinker_gpu_name := $(tinker_gpu_name)$(precision_suffix)
	all_tests_name := $(all_tests_name)$(precision_suffix)
	libtinkergpu_name := $(libtinkergpu_name)$(precision_suffix)
	src_objdir := $(src_objdir)$(precision_suffix)
endif

.PHONY: all
all: $(tinker_gpu_name)

############################################################
#                                                          #
# some utilities                                           #
#                                                          #
############################################################

.PHONY: print-config
print-config:
	@echo 'config file      =' $(config)
	@echo 'opt level        =' $(opt)
	@echo 'OS               =' $(os)
	@echo 'source directory =' $(topdir)
	@echo 'executable       =' $(tinker_gpu_name) $(all_tests_name)

.PHONY: dirs
dirs:
	@mkdir -p $(extdir)/fmt
	@mkdir -p $(extdir)/tinker
	@mkdir -p $(src_objdir)/test
	@mkdir -p $(src_objdir)/util
	@mkdir -p $(src_objdir)/gpu

.PHONY: clean
clean:
	rm -rf ext src* tinker.gpu* all.tests* *.a *dSYM

########################################################################
#                                                                      #
# libtinkergpu.a                                                       #
#                                                                      #
########################################################################

topdir := $(realpath $(src))
extdir := ext/ext

include_flags := -I$(topdir)/include
include_flags += -I$(topdir)/ext
include_flags += -I$(topdir)/ext/ext
include_flags += -I$(cudadir)/include

cxx_common_flags := -std=c++11 -D$(fortran_name_macro) $(precision_flag)
cxx_common_flags += $(include_flags) $(cxx_opt_flags)
cxx_compile_flags := -x c++ -c $(cxx_common_flags)

acc_common_flags := $(cxx_common_flags)
acc_compile_flags += -acc verystrict -c -Minfo=accel $(acc_common_flags)

cxx_link_flags := $(cxx_common_flags)
cxx_link_flags += $(os_link_flags) $(libtinker) $(fortrtflag) $(fftwflag)

acc_link_flags := -acc -Mcuda -Mcudalib=cufft
acc_link_flags += $(acc_common_flags) $(libtinker) $(fortrtflag) $(fftwflag)

ext_cc_objs := $(extdir)/fmt/format.o $(extdir)/fmt/posix.o
ext_cpp_objs := $(extdir)/tinker/tinker_mod.o $(extdir)/tinker/tinker_rt.o

main_cc_objs := $(src_objdir)/test/test_all.o $(src_objdir)/tinker_gpu.o

# $(topdir)/src/util/foo.cpp -> src-dp/util/foo.o
lib_cpp_files := $(shell find $(topdir)/src -type f -name '*.cpp' \
-not -path '$(topdir)/$(extdir)/*' \
-not -path '$(topdir)/src/test/*')
lib_cpp_objs := $(subst $(topdir)/src,$(src_objdir),$(lib_cpp_files))
lib_cpp_objs := $(patsubst %.cpp,%.o,$(lib_cpp_objs))

# $(topdir)/src/gpu/foo.cc -> src-dp/gpu/foo.o
lib_acc_files := $(wildcard $(topdir)/src/gpu/*.cc)
lib_acc_objs := $(subst $(topdir)/src,$(src_objdir),$(lib_acc_files))
lib_acc_objs := $(patsubst %.cc,%.o,$(lib_acc_objs))

# $(topdir)/src/test/foo.cpp -> src-dp/test/foo.o
test_cpp_files := $(wildcard $(topdir)/src/test/*.cpp)
test_cpp_objs := $(subst $(topdir)/src,$(src_objdir),$(test_cpp_files))
test_cpp_objs := $(patsubst %.cpp,%.o,$(test_cpp_objs))

all_cc_objs := $(main_cc_objs) $(ext_cc_objs)
all_cpp_objs := $(ext_cpp_objs) $(lib_cpp_objs) $(test_cpp_objs)

# targets: normal-prerequisites | order-only-prerequisites
$(all_cc_objs): %.o: $(topdir)/%.cc | dirs
	$(CXX) $(cxx_compile_flags) $< -o $@
$(all_cpp_objs): %.o: $(topdir)/%.cpp | dirs
	$(CXX) $(cxx_compile_flags) $< -o $@
$(lib_acc_objs): %.o: $(topdir)/%.cc | dirs
	$(ACC) $(acc_compile_flags) $< -o $@

$(libtinkergpu_name).a: $(ext_cc_objs) $(ext_cpp_objs) $(lib_cpp_objs) $(lib_acc_objs)
	$(ar_flags) $@ $^
	$(ranlib_flags) $@

########################################################################
#                                                                      #
# tinker.gpu                                                           #
#                                                                      #
########################################################################

$(tinker_gpu_name): $(src_objdir)/tinker_gpu.o $(libtinkergpu_name).a
#	$(CXX) $(cxx_link_flags) $^ -o $@
	$(ACC) $^ $(acc_link_flags) -o $@

########################################################################
#                                                                      #
# tests                                                                #
#                                                                      #
########################################################################

.PHONY: test
test: $(all_tests_name)
	$(PWD)/$(all_tests_name) [ff] --durations yes

$(all_tests_name): $(src_objdir)/test/test_all.o $(test_cpp_objs) $(libtinkergpu_name).a
#	$(CXX) $(cxx_link_flags) $^ -o $@
	$(ACC) $^ $(acc_link_flags) -o $@
