# printf "set (__T9_TEST_SRC\n%s\n)\n" "$(LC_COLLATE=C ls -1a *.cpp)" > cmakesrc.txt
cmake_minimum_required (VERSION 3.15)


## test files
include (cmakesrc.txt)
add_library (__t9_all_tests_o OBJECT ${__T9_TEST_SRC})
target_include_directories (__t9_all_tests_o SYSTEM PRIVATE "${PROJECT_SOURCE_DIR}/ext/ext")
target_link_libraries (__t9_all_tests_o PRIVATE __t9_intf PRIVATE tinkerFToCpp)
