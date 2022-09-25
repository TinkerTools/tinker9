add_executable (tinker9)
target_link_libraries (tinker9
   __t9_main_o
   tinker9_acc
   tinker9_cpp
   tinker9_version
   tinker9_f
   tinkerFToCpp
   pthread
)


add_executable (all.tests)
target_link_libraries (all.tests
   __t9_all_tests_o
   tinker9_acc
   tinker9_cpp
   tinker9_version
   tinker9_f
   tinkerFToCpp
   pthread
)
