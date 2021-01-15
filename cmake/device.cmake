add_library (tinker9_main OBJECT src-main/main_tinker9.cpp)
set_target_properties (tinker9_main PROPERTIES
   CXX_STANDARD
      11
)
target_compile_definitions (tinker9_main PRIVATE "${macro_defs}")
target_include_directories (tinker9_main SYSTEM PRIVATE "${TINKER9_SYS_INC_PATH}")
target_include_directories (tinker9_main PRIVATE "${TINKER9_INTERNAL_INC_PATH}")


separate_arguments (LIST_CXX_FLAGS_DEBUG NATIVE_COMMAND ${CMAKE_CXX_FLAGS_DEBUG})
separate_arguments (LIST_CXX_FLAGS_RELEASE NATIVE_COMMAND ${CMAKE_CXX_FLAGS_RELEASE})


########################################################################


add_custom_target (tinker9 ALL
   DEPENDS
      tinker9_main
      cm9a
      tinker9_cu
      tinker9_cpp
      tinker9_f
      LIBTINKER
      LIBFFTW
      LIBFFTW_THREADS
   COMMAND
      "${TINKER9_ACC_COMPILER}"
      CUDA_HOME=${CUDA_DIR}
      "$<$<CONFIG:DEBUG>:${LIST_CXX_FLAGS_DEBUG}>"
      "$<$<CONFIG:RELEASE>:${LIST_CXX_FLAGS_RELEASE}>"
      -o tinker9
      $<TARGET_OBJECTS:tinker9_main>
      "-Wl,--start-group"
      $<TARGET_FILE:tinker9_acc>
      $<TARGET_FILE:tinker9_cu>
      $<TARGET_FILE:tinker9_cpp>
      $<TARGET_FILE:tinker9_f>
      "-Wl,--end-group"
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


########################################################################


add_custom_target (all.tests ALL
   DEPENDS
      all_tests_o
      cm9a
      tinker9_cu
      tinker9_cpp
      tinker9_f
      LIBTINKER
      LIBFFTW
      LIBFFTW_THREADS
   COMMAND
      "${TINKER9_ACC_COMPILER}"
      CUDA_HOME=${CUDA_DIR}
      "$<$<CONFIG:DEBUG>:${LIST_CXX_FLAGS_DEBUG}>"
      "$<$<CONFIG:RELEASE>:${LIST_CXX_FLAGS_RELEASE}>"
      -o all.tests
      $<TARGET_OBJECTS:all_tests_o>
      "-Wl,--start-group"
      $<TARGET_FILE:tinker9_acc>
      $<TARGET_FILE:tinker9_cu>
      $<TARGET_FILE:tinker9_cpp>
      $<TARGET_FILE:tinker9_f>
      "-Wl,--end-group"
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
