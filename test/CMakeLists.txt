# printf "set (__T9_TEST_SRC\n%s\n)\n" "$(LC_COLLATE=C ls -1a *.cpp)" > cmakesrc.txt
cmake_minimum_required (VERSION 3.15)


if ("${CMAKE_CXX_COMPILE_FEATURES}" MATCHES cxx_std_14)
   if (14 LESS ${T9_CPPSTD})
      set (__T9TEST_STD ${T9_CPPSTD})
   else ()
      set (__T9TEST_STD 14)
   endif ()

   # Catch2 v3
   set (__T9_USE_CATCH2_V3 ON)
   include (ExternalProject)
   ExternalProject_Add (Catch2
      SOURCE_DIR  "${PROJECT_SOURCE_DIR}/ext/ext/catch2v3"
      CMAKE_ARGS  -DCMAKE_BUILD_TYPE=release
                  -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/catch2
                  -DCMAKE_CXX_STANDARD=14
   )
else()
   # Catch2 v2
   set (__T9TEST_STD 11)
   set (__T9_USE_CATCH2_V3 OFF)
endif()


## test files
include (cmakesrc.txt)
add_library (__t9_all_tests_o OBJECT ${__T9_TEST_SRC})
target_link_libraries (__t9_all_tests_o PRIVATE __t9_intf PRIVATE tinkerFToCpp)
set_target_properties (__t9_all_tests_o PROPERTIES CXX_STANDARD ${__T9TEST_STD})
if (__T9_USE_CATCH2_V3)
   target_include_directories (__t9_all_tests_o
      SYSTEM PRIVATE
         "${CMAKE_BINARY_DIR}/catch2/include")
   target_link_libraries (__t9_all_tests_o
      PUBLIC
         "${CMAKE_BINARY_DIR}/catch2/lib/libCatch2Main.a"
      PUBLIC
         "${CMAKE_BINARY_DIR}/catch2/lib/libCatch2.a"
   )
else ()
   target_include_directories (__t9_all_tests_o SYSTEM
      PRIVATE
         "${PROJECT_SOURCE_DIR}/ext/ext")
endif ()
