ACC := pgc++
ifeq ($(opt),debug)
    acc_opt_flags__ += -O0 -g -traceback -ta=tesla
else ifeq ($(opt),release)
    acc_opt_flags__ += -O3 -ta=tesla:fastmath
else ifeq ($(opt),profile)
    acc_opt_flags__ += -O3 -g -ta=tesla:fastmath
endif
acc_flags__ := -std=c++11 -acc verystrict -Minfo=accel $(shared_flags__) $(acc_opt_flags__)

acc_depend_flags__ := $(acc_flags__) -MM
acc_compile_flags__ := $(acc_flags__) -c -fpic

$(acc_dependency__): %.d: $(top_dir__)/%.cpp
	printf %s $(@D)/ > $@
	$(ACC) $< $(acc_depend_flags__) >> $@
$(all_acc_objs__): %.o: $(top_dir__)/%.cpp
	$(ACC) $< $(acc_compile_flags__) -o $@
