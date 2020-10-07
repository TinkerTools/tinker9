ACC := $(CXX)
acc_opt_flags__ := $(cxx_opt_flags__)
acc_compile_flags__ := $(cxx_compile_flags__)

$(all_acc_objs__): %.o: $(top_dir__)/%.cpp
	$(ACC) $< $(acc_compile_flags__) -o $@
