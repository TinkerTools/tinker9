add_executable (tinker9 src-main/main_tinker9.cpp)
add_dependencies (tinker9 src-acc src-libtinker)
set_target_properties (tinker9 PROPERTIES
   CXX_STANDARD
      ${T9_CPPSTD}
)
target_compile_definitions (tinker9 PRIVATE "${T9_DEFS}")
target_include_directories (tinker9 SYSTEM PRIVATE "${T9_SYS_INCPATH}")
target_include_directories (tinker9 PRIVATE "${T9_INCPATH}")
set (T9_EXTLIBS pthread t9_ltinker t9_lfftw t9_lfftw_threads)
target_link_libraries (tinker9
   "-Wl,--start-group"
   tinker9_EP_acc
   tinker9_cpp
   tinker9_f
   "-Wl,--end-group"
   "${T9_EXTLIBS}"
)


########################################################################


add_executable (all.tests)
add_dependencies (all.tests src-acc src-libtinker)
target_link_libraries (all.tests
   "-Wl,--start-group"
   __t9_all_tests_o
   tinker9_EP_acc
   tinker9_cpp
   tinker9_f
   "-Wl,--end-group"
   "${T9_EXTLIBS}"
)
