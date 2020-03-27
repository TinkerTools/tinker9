# -fp-model=consistent is important to expf() calls in PME
# otherwise a lot of numbers would become NAN
# details
# https://software.intel.com/en-us/forums/intel-fortran-compiler/topic/749252

CXX := icpc
ifeq ($(opt),debug)
   cxx_opt_flags__ += -O0 -g
else ifeq ($(opt),release)
   cxx_opt_flags__ += -O3 -fp-model=consistent
else ifeq ($(opt),profile)
   cxx_opt_flags__ += -O3 -fp-model=consistent -g
endif
cxx_flags__ += -std=c++11 $(shared_flags__) $(cxx_opt_flags__)
cxx_flags__ += $(include_default_platform__)

cxx_compile_flags__ := $(cxx_flags__) -c
#161: unrecognized #pragma
cxx_compile_flags__ += -diag-disable 161

$(all_cpp_objs__): %.o: $(top_dir__)/%.cpp
	$(CXX) $< $(cxx_compile_flags__) -o $@
