add_executable (tinker9 src/main/tinker9.cpp)
set_target_properties (tinker9 PROPERTIES
   CXX_STANDARD
      ${T9_CPPSTD}
)
target_compile_definitions (tinker9 PRIVATE "${T9_DEFS}")
target_include_directories (tinker9 SYSTEM PRIVATE "${T9_SYS_INCPATH}")
target_include_directories (tinker9 PRIVATE "${T9_INCPATH}")
set (__T9_EXTLIBS pthread t9_lfftw t9_lfftw_threads)
target_link_libraries (tinker9
   "-Wl,--start-group"
   tinker9_acc
   tinker9_cpp
   tinker9_f
   "-Wl,--end-group"
   tinkerFToCPP
   "${__T9_EXTLIBS}"
)


########################################################################


add_executable (all.tests)
target_link_libraries (all.tests
   "-Wl,--start-group"
   __t9_all_tests_o
   tinker9_acc
   tinker9_cpp
   tinker9_f
   "-Wl,--end-group"
   tinkerFToCPP
   "${__T9_EXTLIBS}"
)
