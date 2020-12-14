## test cpp files
file (GLOB TEST_CPP "${PROJECT_SOURCE_DIR}/src/test/*.cpp")


add_library (all_tests_o OBJECT ${MAINTEST_CPP} ${TEST_CPP})
target_compile_definitions (all_tests_o PRIVATE ${macro_defs})
set_target_properties (all_tests_o PROPERTIES
   CXX_STANDARD 11
   CXX_EXTENSIONS OFF
)
target_include_directories (all_tests_o SYSTEM PRIVATE ${comm_sys_inc_path})
target_include_directories (all_tests_o PRIVATE
   ${proj_internal_inc_path}
   "${PROJECT_SOURCE_DIR}/src/test"
)
if (NOT HOST)
   target_include_directories (all_tests_o PRIVATE
      "${PROJECT_SOURCE_DIR}/include/syntax/acc"
   )
endif ()


if (HOST)
add_executable (all.tests)
target_link_libraries (all.tests all_tests_o tinker9_host ${EXT_LIBS})
else ()
add_custom_target (all.tests ALL
   DEPENDS
      all_tests_o
      tinker9_acc
      tinker9_cu
      tinker9_f
      tinker9_cpp
      LIBTINKER
      LIBFFTW
      LIBFFTW_THREADS
   COMMAND
      ${CMAKE_CXX_COMPILER} CUDA_HOME=${CUDA_DIR}
      "$<$<CONFIG:DEBUG>:${LIST_CXX_FLAGS_DEBUG}>"
      "$<$<CONFIG:RELEASE>:${LIST_CXX_FLAGS_RELEASE}>"
      -o all.tests
      $<TARGET_OBJECTS:all_tests_o>
      $<TARGET_OBJECTS:tinker9_acc>
      $<TARGET_OBJECTS:tinker9_cu>
      $<TARGET_OBJECTS:tinker9_f>
      $<TARGET_OBJECTS:tinker9_cpp>
      $<TARGET_FILE:LIBTINKER>
      $<TARGET_FILE:LIBFFTW>
      $<TARGET_FILE:LIBFFTW_THREADS>
      "-L$<JOIN:${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES},;-L>"
      "-l$<JOIN:${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES},;-l>"
      -acc -Mcudalib=cufft,cublas
      $<$<CONFIG:DEBUG>:-ta=tesla:lineinfo${CCLIST4}>
      $<$<CONFIG:RELEASE>:-ta=tesla:fastmath${CCLIST4}>
   COMMAND_EXPAND_LISTS
)
endif ()


########################################################################


add_custom_target (test
   COMMAND
      ./all.tests info
   COMMAND
      ./all.tests [ff],[util] -a --durations yes --order rand --rng-seed time
   DEPENDS
      all.tests
)
