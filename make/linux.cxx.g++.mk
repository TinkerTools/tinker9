CXX := g++
ifeq ($(opt),debug)
   cxx_opt_flags__ += -O0 -g
else ifeq ($(opt),release)
   cxx_opt_flags__ += -O3
else ifeq ($(opt),profile)
   cxx_opt_flags__ += -O3 -g
endif
cxx_flags__ += -std=c++11 $(shared_flags__) $(cxx_opt_flags__)

cxx_compile_flags__ := $(cxx_flags__) -c

$(all_cpp_objs__): %.o: $(top_dir__)/%.cpp
	$(CXX) $< $(cxx_compile_flags__) -o $@
