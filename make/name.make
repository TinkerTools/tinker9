tinker_gpu_exec__ := tinker.gpu
all_tests_exec__ := all.tests
libtinkergpu_stem__ := tinkergpu
libtinkergpu_name__ := lib$(libtinkergpu_stem__)
src_obj_dir__ := src
ext_dir__ := ext/ext


# acc file
# $(top_dir__)/src/acc_foo.cpp -> src/acc_foo.o
lib_acc_files__ := $(wildcard $(top_dir__)/src/acc_*.cpp)
lib_acc_objs__ := $(subst $(top_dir__)/src,$(src_obj_dir__),$(lib_acc_files__))
lib_acc_objs__ := $(lib_acc_objs__:%.cpp=%.o)
# lib_acc_objs__ := $(patsubst %.cpp,%.o,$(lib_acc_objs__))
# all acc
all_acc_objs__ := $(lib_acc_objs__)
acc_dependency__ := $(all_acc_objs__:%.o=%.d)


# cc files
ext_cc_objs__ := $(ext_dir__)/fmt/format.o $(ext_dir__)/fmt/posix.o
# all cc
all_cc_objs__ := $(ext_cc_objs__)
cc_dependency__ := $(all_cc_objs__:%.o=%.d)


# cpp files
# main
main_cpp_files__ := $(top_dir__)/src/main_tinker_gpu.cpp $(top_dir__)/src/test/main_all_tests.cpp
main_cpp_objs__ := $(subst $(top_dir__)/src,$(src_obj_dir__),$(main_cpp_files__))
main_cpp_objs__ := $(main_cpp_objs__:%.cpp=%.o)
# test
# $(top_dir__)/src/test/foo.cpp -> src/test/foo.o
test_cpp_files__ := $(shell find $(top_dir__)/src/test -type f -name '*.cpp' \
-not -path '$(top_dir__)/src/test/main_*.cpp')
test_cpp_objs__ := $(subst $(top_dir__)/src,$(src_obj_dir__),$(test_cpp_files__))
test_cpp_objs__ := $(test_cpp_objs__:%.cpp=%.o)
# misc.
# $(top_dir__)/src/foo.cpp -> src/foo.o
lib_cpp_files__ := $(shell find $(top_dir__)/src -type f -name '*.cpp' \
-not -path '$(top_dir__)/src/main_*.cpp' \
-not -path '$(top_dir__)/src/acc_*.cpp' \
-not -path '$(top_dir__)/src/test/*')
lib_cpp_objs__ := $(subst $(top_dir__)/src,$(src_obj_dir__),$(lib_cpp_files__))
lib_cpp_objs__ := $(lib_cpp_objs__:%.cpp=%.o)
# all cpp
all_cpp_objs__ := $(main_cpp_objs__) $(test_cpp_objs__) $(lib_cpp_objs__)
cpp_dependency__ := $(all_cpp_objs__:%.o=%.d)
