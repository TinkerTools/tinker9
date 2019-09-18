ACC := $(CXX)
acc_opt_flags__ := $(cxx_opt_flags__)
acc_depend_flags__ := $(cxx_depend_flags__)
acc_compile_flags__ := $(cxx_compile_flags__)

$(acc_dependency__): %.d: $(top_dir__)/%.cpp
	printf %s $(@D)/ > $@
	$(ACC) $< $(acc_depend_flags__) >> $@
$(all_acc_objs__): %.o: $(top_dir__)/%.cpp
	$(ACC) $< $(acc_compile_flags__) -o $@
