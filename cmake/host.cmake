add_executable (tinker9 src-main/main_tinker9.cpp)
add_dependencies (tinker9 cmtinker9acc)
set_target_properties (tinker9 PROPERTIES
   CXX_STANDARD
      11
)
target_compile_definitions (tinker9 PRIVATE "${T9_DEFS}")
target_include_directories (tinker9 SYSTEM PRIVATE "${T9_SYS_INCPATH}")
target_include_directories (tinker9 PRIVATE "${T9_INCPATH}")
set (T9_EXTLIBS pthread LIBTINKER LIBFFTW LIBFFTW_THREADS)
if (PREC STREQUAL "m" OR PREC STREQUAL "s")
   list (APPEND T9_EXTLIBS LIBFFTWF LIBFFTWF_THREADS)
endif ()
target_link_libraries (tinker9
   "-Wl,--start-group"
   $<TARGET_FILE:tinker9_EP_acc>
   tinker9_cpp
   tinker9_f
   "-Wl,--end-group"
   "${T9_EXTLIBS}"
)


########################################################################


add_executable (all.tests)
add_dependencies (all.tests cmtinker9acc)
target_link_libraries (all.tests
   "-Wl,--start-group"
   all_tests_o
   $<TARGET_FILE:tinker9_EP_acc>
   tinker9_cpp
   tinker9_f
   "-Wl,--end-group"
   "${T9_EXTLIBS}"
)
