add_library (tinker9_main OBJECT src-main/main_tinker9.cpp)
set_target_properties (tinker9_main PROPERTIES
   CXX_STANDARD
      11
)
target_compile_definitions (tinker9_main PRIVATE "${T9_DEFS}")
target_include_directories (tinker9_main SYSTEM PRIVATE "${T9_SYS_INCPATH}")
target_include_directories (tinker9_main PRIVATE "${T9_INCPATH}")


separate_arguments (__T9_DEBUG_FLAGS NATIVE_COMMAND ${CMAKE_CXX_FLAGS_DEBUG})
separate_arguments (__T9_RELEASE_FLAGS NATIVE_COMMAND ${CMAKE_CXX_FLAGS_RELEASE})


## Compute Capability 60,70 -> ,cc60,cc70
set (__T9_CC4) # ""
foreach (var ${T9_CUCCLIST})
   string (APPEND __T9_CC4 ",cc${var}")
endforeach () # ,cc60,cc70


########################################################################


add_custom_target (tinker9 ALL
   DEPENDS
      tinker9_main
      cmtinker9acc
      tinker9_cu
      tinker9_cpp
      tinker9_f
      t9_ltinker
      t9_lfftw
      t9_lfftw_threads
   COMMAND
      "${T9_ACC_COMPILER}"
      CUDA_HOME=${CUDA_DIR}
      "$<$<CONFIG:DEBUG>:${__T9_DEBUG_FLAGS}>"
      "$<$<CONFIG:RELEASE>:${__T9_RELEASE_FLAGS}>"
      -o tinker9
      $<TARGET_OBJECTS:tinker9_main>
      "-Wl,--start-group"
      $<TARGET_FILE:tinker9_EP_acc>
      $<TARGET_FILE:tinker9_cu>
      $<TARGET_FILE:tinker9_cpp>
      $<TARGET_FILE:tinker9_f>
      "-Wl,--end-group"
      $<TARGET_FILE:t9_ltinker>
      $<TARGET_FILE:t9_lfftw>
      $<TARGET_FILE:t9_lfftw_threads>
      "-L$<JOIN:${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES},;-L>"
      "-l$<JOIN:${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES},;-l>"
      -acc -Mcudalib=cufft,cublas
      $<$<CONFIG:DEBUG>:-ta=tesla:lineinfo${__T9_CC4}>
      $<$<CONFIG:RELEASE>:-ta=tesla:fastmath${__T9_CC4}>
   COMMAND_EXPAND_LISTS
)


########################################################################


add_custom_target (all.tests ALL
   DEPENDS
      __t9_all_tests_o
      cmtinker9acc
      tinker9_cu
      tinker9_cpp
      tinker9_f
      t9_ltinker
      t9_lfftw
      t9_lfftw_threads
   COMMAND
      "${T9_ACC_COMPILER}"
      CUDA_HOME=${CUDA_DIR}
      "$<$<CONFIG:DEBUG>:${__T9_DEBUG_FLAGS}>"
      "$<$<CONFIG:RELEASE>:${__T9_RELEASE_FLAGS}>"
      -o all.tests
      $<TARGET_OBJECTS:__t9_all_tests_o>
      "-Wl,--start-group"
      $<TARGET_FILE:tinker9_EP_acc>
      $<TARGET_FILE:tinker9_cu>
      $<TARGET_FILE:tinker9_cpp>
      $<TARGET_FILE:tinker9_f>
      "-Wl,--end-group"
      $<TARGET_FILE:t9_ltinker>
      $<TARGET_FILE:t9_lfftw>
      $<TARGET_FILE:t9_lfftw_threads>
      "-L$<JOIN:${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES},;-L>"
      "-l$<JOIN:${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES},;-l>"
      -acc -Mcudalib=cufft,cublas
      $<$<CONFIG:DEBUG>:-ta=tesla:lineinfo${__T9_CC4}>
      $<$<CONFIG:RELEASE>:-ta=tesla:fastmath${__T9_CC4}>
   COMMAND_EXPAND_LISTS
)
