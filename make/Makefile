##/ # Build Tinker GPU: Makefile, Options, and Targets
##/
##/ ## Introduction
##/ This is a self-documenting Makefile. After being processed by `doxygen`,
##/ the documentation will appear in the `Related Pages` section. Simply run
##/ ```
##/ $> make help
##/ ```
##/ can also print the documentation in terminal window. An empty line is
##/ mandatory at the end of each paragraph of the documentation.
##/
##/ Every option has a default value and is documented as
##/ `option=default_value`. In order to override the default value(s), set
##/ them in the command line as follows:
##/ `$> make option1=value1 option2=value2 targets`
##/
##/ ------
##/
##/ ## Options
##/


define errline

**************************************************

endef


# top-level directory
makefile_path__ := $(abspath $(lastword $(MAKEFILE_LIST)))
makefile_dir__ := $(dir $(makefile_path__))
# `make' will not complain if file `top_dir.mk' does not exist
-include top_dir.mk
ifeq (,$(top_dir__))
    top_dir__ := $(realpath $(makefile_dir__)/..)
endif


# linux (Linux) / darwin (macOS)
os__ := $(shell uname -s 2>/dev/null | tr "[:upper:]" "[:lower:]")
# so (Linux) / dylib (macOS)
ifeq ($(os__),linux)
    shared_lib_suffix__ := so
else ifeq ($(os__),darwin)
    shared_lib_suffix__ := dylib
endif


shared_flags__ :=
link_flags__ :=


##/ ### opt=release
##/ Optimization level. Other valid values are `debug` and `profile`.
##/
opt := release
ifeq ($(opt),debug)
    shared_flags__ += -DDEBUG=1
else ifeq ($(opt),release)
    shared_flags__ += -DNDEBUG
else ifeq ($(opt),profile)
    shared_flags__ += -DNDEBUG
endif


##/ ### prec=real8
##/ Precision of the floating point numbers.
##/   - These equivalent values all compile the code to double precision:
##/     `real8`, `8`, `64`, `double`, `d`.
##/   - These equivalent values all compile the code to single precision:
##/     `real4`, `4`, `32`, `float`, `single`, `s`.
##/
prec := real8
ifeq ($(prec),$(filter $(prec),real4 4 32 float single s))
    prec__ := real4__
    shared_flags__ += -DTINKER_SINGLE_PRECISION
else ifeq ($(prec),$(filter $(prec),real8 8 64 double d))
    prec__ := real8__
    shared_flags__ += -DTINKER_DOUBLE_PRECISION
else
    $(error Error! Flag prec=$(prec) was set incorrectly)
endif


##/ ### host=1
##/ Flag to compile to GPU version or CPU version.
##/   - These equivalent values all compile the code to GPU version:
##/     `false`, `off`, `0`.
##/   - These equivalent values all compile the code to CPU version:
##/     `true`, `on`, `1`.
##/
host := 1
ifeq ($(host),$(filter $(host),false off 0))
    host_only__ := false
else ifeq ($(host),$(filter $(host),true on 1))
    host_only__ := true
    shared_flags__ += -DTINKER_HOST
else
    $(error Error! Flag host=$(host) was set incorrectly)
endif


##/ ### build=build
##/ The Makefile is capable of creating a directory with name `$build` in the
##/ current working directory. This option enables the user to change its
##/ default name.
##/
build := build


##/ ### tinker_dir (NO DEFAULT)
##/ The directory in which user compiled the `libtinker`.
##/ If this value is not set by a command line argument, the Makefile will
##/ attempt to locate `libtinker` under `$HOME/tinker/source` and emit a fatal
##/ error when it fails to find it. No check will be performed if this option
##/ is set via command line arguments.
##/
shared_libtinker__ := libtinker.$(shared_lib_suffix__)
ifdef tinker_dir
    tinker_dir__ := $(tinker_dir)
else
    ifneq (,$(wildcard $(HOME)/tinker/source/$(shared_libtinker__)))
        tinker_dir__ := $(HOME)/tinker/source
    endif
endif
ifndef tinker_dir__
    $(warning $(errline)Warning! tinker_dir must be supplied.$(errline))
endif
link_flags__ += -L. -ltinker


##/ ### fftw_dir / fftw_include / fftw_lib (NO DEFAULT)
##/ No default values are set for these three options.
##/
##/ `fftw_dir` is the top-level `FFTW` installation, under which
##/ `include/fftw3.h` and `lib/libfftw3` are expected to be found.
##/ If this value is not set by a command line argument, the Makefile will
##/ attempt to locate them under `/usr/local`. However, no fatal error will
##/ be emitted should the Makefile fail to find them.
##/
##/ `fftw_include` and `fftw_lib` will override the values set by `fftw_dir`.
##/ Fatal errors will be emitted if `fftw_x` (`x = include/lib`) is neither set
##/ by `fftw_x` nor `fftw_dir`.
##/
ifdef fftw_dir
    fftw_include__ := $(fftw_dir)/include
    fftw_lib__ := $(fftw_dir)/lib
else
    ifneq (,$(wildcard /usr/local/include/fftw3.h))
        fftw_include__ := /usr/local/include
    endif
    ifneq (,$(wildcard /usr/local/lib/libfftw3.*))
        fftw_lib__ := /usr/local/lib
    endif
endif


ifdef fftw_include
    fftw_include__ := $(fftw_include)
endif
ifndef fftw_include__
    $(warning $(errline)Warning! fftw_include or fftw_dir must be supplied.$(errline))
endif


ifdef fftw_lib
    fftw_lib__ := $(fftw_lib)
endif
ifndef fftw_lib__
    $(warning $(errline)Warning! fftw_lib or fftw_dir must be supplied.$(errline))
endif
link_flags__ += -L$(fftw_lib__) -lfftw3 -lfftw3_threads
ifeq ($(prec__),real4__)
    link_flags__ += -lfftw3f
endif


##/ ### cuda_dir=/usr/local/cuda
##/ Top-level `CUDA` installation directory, under which `include` and `lib`
##/ directories are expected to be found.
##/
cuda_dir := /usr/local/cuda


# include directories
include_flags__ := -I$(top_dir__)/include
include_flags__ += -I$(top_dir__)/ext
include_flags__ += -I$(top_dir__)/ext/ext
ifeq ($(host_only__),false)
    include_flags__ += -I$(cuda_dir)/include
else
    include_flags__ += -I$(fftw_include__)
endif
shared_flags__ += $(include_flags__)


src_obj_dir__ := src
ext_obj_dir__ := ext/ext


# cpp
# main
main_cpp_objs__ := $(top_dir__)/src/main_tinker_gpu.cpp $(top_dir__)/src/test/main_all_tests.cpp
main_cpp_objs__ := $(subst $(top_dir__)/src,$(src_obj_dir__),$(main_cpp_objs__))
main_cpp_objs__ := $(main_cpp_objs__:%.cpp=%.o)
# ext
ext_cpp_objs__ := $(ext_obj_dir__)/fmt/format.o $(ext_obj_dir__)/fmt/posix.o
# test
# $(top_dir__)/src/test/foo.cpp -> src/test/foo.o
test_cpp_objs__ := $(shell find $(top_dir__)/src/test -type f -name '*.cpp' \
-not -path '$(top_dir__)/src/test/main_*.cpp')
test_cpp_objs__ := $(subst $(top_dir__)/src,$(src_obj_dir__),$(test_cpp_objs__))
test_cpp_objs__ := $(test_cpp_objs__:%.cpp=%.o)
# misc.
# $(top_dir__)/src/foo.cpp -> src/foo.o
cpp_objs__ := $(shell find $(top_dir__)/src -type f -name '*.cpp' \
-not -path '$(top_dir__)/src/main_*.cpp' \
-not -path '$(top_dir__)/src/acc_*.cpp' \
-not -path '$(top_dir__)/src/test/*')
cpp_objs__ := $(subst $(top_dir__)/src,$(src_obj_dir__),$(cpp_objs__))
cpp_objs__ := $(cpp_objs__:%.cpp=%.o)
# acc
# $(top_dir__)/src/acc_foo.cpp -> src/acc_foo.o
acc_objs__ := $(wildcard $(top_dir__)/src/acc_*.cpp)
acc_objs__ := $(subst $(top_dir__)/src,$(src_obj_dir__),$(acc_objs__))
acc_objs__ := $(acc_objs__:%.cpp=%.o)
# all cpp acc
lib_cpp_objs__ := $(ext_cpp_objs__) $(cpp_objs__)
all_cpp_objs__ := $(main_cpp_objs__) $(test_cpp_objs__) $(lib_cpp_objs__)
all_acc_objs__ := $(acc_objs__)


##/ ### fortran_config / cxx_config / acc_config / link_config (UNSPECIFIED)
##/ Default configuration files that contain Fortran/C++/OpenAcc/linker
##/ related flags. This is the mechanism we adopt to extend the multi-compiler
##/ support. Tested compilers are tabulated below.
##/
##/ | OS    | Fortran                 | C++                | OpenACC         | Linker          |
##/ |:-----:|:-----------------------:|:------------------:|:---------------:|:---------------:|
##/ | Linux | gfortran 5.4.0          | g++                | pgc++ -ta=telsa | pgc++ -ta=telsa |
##/ | Linux | gfortran 5.4.0          | g++                | g++             | g++             |
##/ | macOS | gfortran 8.3.0 homebrew | clang++ xcode 10.1 | clang++         | clang++         |
##/
ifeq ($(os__),linux)
    fortran_config := linux.fort.gfortran.5.4.0.mk
    cxx_config := linux.cxx.g++.mk
    ifeq ($(host_only__),false)
        acc_config := linux.acc.device.pgc++.mk
        link_config := linux.link.device.pgc++.mk
    else ifeq ($(host_only__),true)
        acc_config := linux.acc.g++.mk
        link_config := linux.link.g++.mk
    endif
else ifeq ($(os__),darwin)
    fortran_config := macos.fort.gfortran.8.3.0.brew.mk
    cxx_config := macos.cxx.clang++.xcode.10.1.mk
    acc_config := macos.acc.clang++.xcode.10.1.mk
    link_config := macos.link.clang++.xcode.10.1.mk
endif
-include $(fortran_config)
-include $(cxx_config)
-include $(acc_config)
-include $(link_config)
link_flags__ += -ltinkergpu


##/ ------
##/
##/ ## Targets
##/


##/ ### default
##/ The `tinker.gpu` executable.
##/
.PHONY: default
default: tinker.gpu


##/ ### unittest
##/ The `all.tests` executable.
##/
.PHONY: unittest
unittest: all.tests


##/ ### all
##/ Both `default` and `unittest` targets.
##/
.PHONY: all
all: default unittest


##/ ### test
##/ Run the unit tests in a random order.
.PHONY: test
test: all.tests
	$(PWD)/$^ info
	$(PWD)/$^ [ff],[util] --durations yes --order rand --rng-seed time


##/ ### dirs
##/ Create sub-directories inside the `$build` directory.
##/ ### copy_files (depends on dirs)
##/ Copy files to the `$build` directory.
##/ ### create_build (depends on copy_files)
##/ Setup the `$build` directory, including copying the necessary files.
##/ For example, setup a new build directory inside the `tinker.gpu` directory
##/ `build-v3`:
##/
##/ ```bash
##/ $> pwd
##/ #> /home/developer/tinker.gpu
##/ $> make -f make/Makefile create_build build=build-v3
##/ ```
##/
dirs__ := $(build)/$(ext_obj_dir__)/fmt
dirs__ += $(build)/$(src_obj_dir__)/test
$(dirs__):
	@echo TINKER.GPU: Creating directory \`$@\'...
	@mkdir -p $@
.PHONY: dirs
dirs: $(dirs__)
.PHONY: copy_files
copy_files: | dirs
	@cd $(build); ln -s $(top_dir__)/make/Makefile .; cd ..
	@for x in $(top_dir__)/make/*.mk; do cd $(build); ln -s $$x .; cd ..; done
	@echo top_dir__ := $(top_dir__) > $(build)/top_dir.mk
	@echo TINKER.GPU: Installation files copied to \`$(build)\' directory.
	@echo TINKER.GPU: Now run \`cd $(build)\' then \`make\' to compile \`tinker.gpu\'\; \
Run \`make help\' for help.
.PHONY: create_build
create_build: | copy_files


##/ ### info
##/ List some of the Makefile information.
##/
.PHONY: info
info:
	@echo 'Optimization      =' $(opt)
	@echo 'Source Directory  =' $(top_dir__)
	@echo 'Precision         =' $(prec)
	@echo 'Host Only         =' $(host_only__)
	@echo 'OS                =' $(os__)
	@echo '-------------------'
	@echo 'fortran_config    =' $(fortran_config)
	@echo 'cxx_config        =' $(cxx_config)
	@echo 'acc_config        =' $(acc_config)
	@echo 'link_config       =' $(link_config)
	@echo '-------------------'
	@echo 'tinker_dir        =' $(tinker_dir__)
	@echo 'fftw_include      =' $(fftw_include__)
	@echo 'fftw_lib          =' $(fftw_lib__)
	@echo '-------------------'
	@echo 'cxx compile flags =' $(cxx_compile_flags__)
	@echo '-------------------'
	@echo 'acc compile flags =' $(acc_compile_flags__)
	@echo '-------------------'
	@echo 'link flags        =' $(link_flags__)
	@echo '-------------------'


##/ ### doc
##/ Generate documentation with `doxygen`.
##/
git_head__ := $(shell git log -1 --format="%h (%cd)")
.PHONY: doc
doc:
	cp -r $(top_dir__)/README.md $(top_dir__)/doc $(top_dir__)/include .
	grep '^\##/' $(top_dir__)/make/Makefile | sed -e 's|##/ ||' | sed -e 's|##/||' > doc/build2.md
	cp doc/build2.md $(top_dir__)/doc/build2.md
	doxygen -w html doc/header.html doc/footer.html doc/customdoxygen.css-orig
	sed -e "s|Roboto,sans-serif|'Lucida Grande',Geneva,Helvetica,Arial,sans-serif|" doc/customdoxygen.css-orig > doc/customdoxygen.css
	ENV_GIT_HEAD="$(git_head__)" doxygen $(top_dir__)/Doxyfile


##/ ### help
##/ Print the Makefile documentation in the terminal window.
##/
.PHONY: help
help:
	@grep '^\##/' $(top_dir__)/make/Makefile | sed -e 's|##/ ||' | sed -e 's|##/||'

##/ ### clean
##/ Clean up the current build directory.
##/
.PHONY: clean
clean:
	rm -rf all.tests tinker.gpu libtinker* *dSYM *.out
	rm -rf *.md doc include
	rm -rf html/*
	find ext/ -type f -name '*.o' -delete
	find src/ -type f -name '*.o' -delete


$(shared_libtinker__): $(tinker_dir__)/$(shared_libtinker__)
	cp $^ $@
tinker.gpu: $(src_obj_dir__)/main_tinker_gpu.o $(all_acc_objs__) \
libtinkergpu.$(shared_lib_suffix__) | $(shared_libtinker__)
	$(LINK) $^ $(link_flags__) -o $@
all.tests: $(src_obj_dir__)/test/main_all_tests.o \
$(test_cpp_objs__) $(all_acc_objs__) \
libtinkergpu.$(shared_lib_suffix__) | $(shared_libtinker__)
	$(LINK) $^ $(link_flags__) -o $@

##/ ------
##/
##/ ## Maintaining Multiple Build Directories
##/
