CXX := g++
ifeq ($(opt),debug)
    cxx_opt_flags__ += -O0 -g
else ifeq ($(opt),release)
    cxx_opt_flags__ += -O3
else ifeq ($(opt),profile)
    cxx_opt_flags__ += -O3 -g
endif
cxx_flags__ := -std=c++11 $(shared_flags__) $(cxx_opt_flags__)

cxx_depend_flags__ := $(cxx_flags__) -MM
cxx_compile_flags__ := $(cxx_flags__) -c -fpic

$(cpp_dependency__): %.d: $(top_dir__)/%.cpp
	printf %s $(@D)/ > $@
	$(CXX) $< $(cxx_depend_flags__) >> $@
$(all_cpp_objs__): %.o: $(top_dir__)/%.cpp
	$(CXX) $< $(cxx_compile_flags__) -o $@
libtinkergpu.$(shared_lib_suffix__): $(lib_cpp_objs__)
	$(CXX) -shared -o $@ $^
