############################################################
#                                                          #
#  Copy Makefile and inc files to the build directory      #
#                                                          #
############################################################

# release / debug / profile
opt := release

# tinker.gpu directory
src := ../tinker.gpu

# Path to the precompiled libtinker.a
libtinker := $(HOME)/git/tinker/source/libtinker.a

# FFTW
fftwflag := -L$(HOME)/local/fftw/lib -lfftw3 -lfftw3_threads

# Darwin -> darwin
os := $(shell uname -s 2>/dev/null | tr "[:upper:]" "[:lower:]")

# CUDA directory
cudadir := /usr/local/cuda

# Compilers
CXX := nvcc
NVCC := nvcc

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
endif

############################################################
#                                                          #
# default target: all                                      #
#                                                          #
############################################################

.PHONY: all
all: tinker.gpu test.all

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

.PHONY: dirs
dirs:
	@mkdir -p $(extdir)/fmt
	@mkdir -p $(extdir)/tinker
	@mkdir -p src/test
	@mkdir -p src/util

.PHONY: clean
clean:
	rm -rf ext src tinker.gpu test.all *.a *dSYM

########################################################################
#                                                                      #
# libtinkergpu.a                                                       #
#                                                                      #
########################################################################

topdir := $(realpath $(src))
extdir := ext/ext

include_flags := -I$(topdir)/include
include_flags += -isystem $(topdir)/ext
include_flags += -isystem $(topdir)/ext/ext
include_flags += -isystem $(cudadir)/include

cxx_common_flags := -std=c++11 -D$(fortran_name_macro)
cxx_common_flags += $(include_flags) $(cxx_opt_flags)

cxx_compile_flags := -x c++ -c $(cxx_common_flags)

link_flags := $(cxx_common_flags) $(libtinker) $(fortrtflag) $(fftwflag)

ext_cc_objs := $(extdir)/fmt/format.o $(extdir)/fmt/posix.o
ext_cpp_objs := $(extdir)/tinker/tinker.mod.o $(extdir)/tinker/tinker.rt.o

main_cc_objs := src/test/test.all.o
main_cpp_objs := src/tinker.gpu.o

# $(topdir)/src/util/foo.cpp -> src/util/foo.o
lib_cpp_files := $(shell find $(topdir)/src -type f -name '*.cpp' \
-not -path $(topdir)/src/tinker.gpu.cpp \
-not -path '$(topdir)/$(extdir)/*' \
-not -path '$(topdir)/src/test/*')
lib_cpp_objs := $(subst $(topdir)/,,$(lib_cpp_files))
lib_cpp_objs := $(patsubst %.cpp,%.o,$(lib_cpp_objs))

# $(topdir)/src/test/test.foo.cpp -> src/test/test.foo.o
test_cpp_files := $(wildcard $(topdir)/src/test/test.*.cpp)
test_cpp_objs := $(subst $(topdir)/,,$(test_cpp_files))
test_cpp_objs := $(patsubst %.cpp,%.o,$(test_cpp_objs))

all_cc_objs := $(main_cc_objs)
all_cpp_objs := $(main_cpp_objs)
all_cc_objs += $(ext_cc_objs)
all_cpp_objs += $(ext_cpp_objs) $(lib_cpp_objs) $(test_cpp_objs)

# targets: normal-prerequisites | order-only-prerequisites
$(all_cc_objs): %.o: $(topdir)/%.cc | dirs
	$(CXX) $(cxx_compile_flags) $< -o $@
$(all_cpp_objs): %.o: $(topdir)/%.cpp | dirs
	$(CXX) $(cxx_compile_flags) $< -o $@

libtinkergpu.a: $(ext_cc_objs) $(ext_cpp_objs) $(lib_cpp_objs)
	$(ar_flags) $@ $^
	$(ranlib_flags) $@

########################################################################
#                                                                      #
# tinker.gpu                                                           #
#                                                                      #
########################################################################

tinker.gpu: src/tinker.gpu.o libtinkergpu.a
	$(CXX) $(link_flags) $^ -o $@

########################################################################
#                                                                      #
# tests                                                                #
#                                                                      #
########################################################################

.PHONY: test
test: test.all
	@$(PWD)/test.all ~[act]

test.all: src/test/test.all.o $(test_cpp_objs) libtinkergpu.a
	$(CXX) $(link_flags) $^ -o $@
