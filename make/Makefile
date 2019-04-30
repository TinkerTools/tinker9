############################################################
#                                                          #
#  Copy this Makefile and inc file to the build directory  #
#                                                          #
############################################################

config := default_config
ifeq ($(config),default_config)
	config := gfortran.inc
endif
include $(config)

f77_default := gfortran

############################################################

ifeq ($(F77),f77)
	F77 := $(f77_default)
endif
found_f77 := false

use_gfortran := false
f77_compiler := $(shell echo $(F77) | cut -c 1-8)
ifeq ($(f77_compiler),gfortran)
	use_gfortran := true
	found_f77 := true
	fortran_name_macro := TINKER_GFORTRAN
endif

############################################################
#                                                          #
# OS                                                       #
#                                                          #
############################################################

# Darwin -> darwin
os := $(shell uname -s 2>/dev/null | tr "[:upper:]" "[:lower:]")

ifeq ($(os),linux)
	ifeq ($(opt),debug)
		cxx_opt_flags := -O0 -g -DDEBUG=1
	else ifeq ($(opt),release)
		cxx_opt_flags := -O3 -DNDEBUG
	else ifeq ($(opt),profile)
		cxx_opt_flags := -O3 -DNDEBUG -g
	endif
endif

############################################################
#                                                          #
# default target: all                                      #
#                                                          #
############################################################

.PHONY: all
all: tinker.gpu

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
	@mkdir -p src/util

.PHONY: clean
clean:
	rm -rf ext src tinker.gpu *.a *dSYM

########################################################################
#                                                                      #
# libtinkergpu.a                                                       #
#                                                                      #
########################################################################

topdir := $(src)
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

main_cpp_objs := src/tinker.gpu.o

# $(topdir)/src/util/foo.cpp -> src/util/foo.o
lib_cpp_files := $(shell find $(topdir)/src -type f -name '*.cpp' \
-not -path $(topdir)/src/tinker.gpu.cpp \
-not -path '$(topdir)/$(extdir)/*' \
-not -path '$(topdir)/src/test/*')
lib_cpp_objs := $(subst $(topdir)/,,$(lib_cpp_files))
lib_cpp_objs := $(patsubst %.cpp,%.o,$(lib_cpp_objs))

all_cc_objs := $(ext_cc_objs)
all_cpp_objs += $(main_cpp_objs)
all_cpp_objs += $(ext_cpp_objs) $(lib_cpp_objs)

# targets: normal-prerequisites | order-only-prerequisites
$(all_cc_objs): %.o: $(topdir)/%.cc | dirs
	$(CXX) $(cxx_compile_flags) $< -o $@
$(all_cpp_objs): %.o: $(topdir)/%.cpp | dirs
	$(CXX) $(cxx_compile_flags) $< -o $@

libtinkergpu.a: $(ext_cc_objs) $(ext_cpp_objs) $(lib_cpp_objs)
	ar qc $@ $^
	ranlib $@

########################################################################
#                                                                      #
# tinker.gpu                                                           #
#                                                                      #
########################################################################

tinker.gpu: src/tinker.gpu.o libtinkergpu.a
	$(CXX) $(link_flags) $^ -o $@
