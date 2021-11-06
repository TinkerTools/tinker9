add_library (tinker9_main OBJECT src-main/main_tinker9.cpp)
set_target_properties (tinker9_main PROPERTIES
   CXX_STANDARD
      ${T9_CPPSTD}
)
target_compile_definitions (tinker9_main PRIVATE "${T9_DEFS}")
target_include_directories (tinker9_main SYSTEM PRIVATE "${T9_SYS_INCPATH}")
target_include_directories (tinker9_main PRIVATE "${T9_INCPATH}")


## Compute Capability 60,70 -> ,cc60,cc70
set (__T9_ACC_CCLST4) # ""
foreach (var ${T9_CUCCLIST})
   string (APPEND __T9_ACC_CCLST4 ",cc${var}")
endforeach () # ,cc60,cc70


separate_arguments (__T9_DEVICE_LINK_FLAGS_DEBUG           NATIVE_COMMAND "${CMAKE_CXX_FLAGS_DEBUG}          -ta=tesla:lineinfo${__T9_ACC_CCLST4}")
separate_arguments (__T9_DEVICE_LINK_FLAGS_RELWITHDEBINFO  NATIVE_COMMAND "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -ta=tesla:lineinfo,fastmath${__T9_ACC_CCLST4}")
separate_arguments (__T9_DEVICE_LINK_FLAGS_RELEASE         NATIVE_COMMAND "${CMAKE_CXX_FLAGS_RELEASE}        -ta=tesla:fastmath${__T9_ACC_CCLST4}")
separate_arguments (__T9_DEVICE_LINK_FLAGS_MINSIZEREL      NATIVE_COMMAND "${CMAKE_CXX_FLAGS_MINSIZEREL}     -ta=tesla:fastmath${__T9_ACC_CCLST4}")
## Remove "-Wno-unknown-pragmas"
list (REMOVE_ITEM __T9_DEVICE_LINK_FLAGS_DEBUG "-Wno-unknown-pragmas")
## Replace -Os flag by -O2 flag for MinSizeRel
list (TRANSFORM __T9_DEVICE_LINK_FLAGS_MINSIZEREL REPLACE "-Os" "-O2")


########################################################################


add_custom_target (tinker9 ALL
   BYPRODUCTS
      tinker9
   DEPENDS
      tinker9_main
      tinker9_acc
      tinker9_cu
      tinker9_cpp
      tinker9_f
      tinkerFToCPP
   COMMAND
      "${T9_ACC_COMPILER}"
      CUDA_HOME=${CUDA_DIR}
      -o tinker9
      $<TARGET_OBJECTS:tinker9_main>
      "-Wl,--start-group"
      $<TARGET_FILE:tinker9_acc>
      $<TARGET_FILE:tinker9_cu>
      $<TARGET_FILE:tinker9_cpp>
      $<TARGET_FILE:tinker9_f>
      "-Wl,--end-group"
      $<TARGET_FILE:tinkerFToCPP>
      "-L$<JOIN:${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES},;-L>"
      "-l$<JOIN:${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES},;-l>"
      -acc -Mcudalib=cufft,cublas
      ${__T9_DEVICE_LINK_FLAGS_${__T9_BUILD_TYPE}}
   COMMAND_EXPAND_LISTS
)


########################################################################


add_custom_target (all.tests ALL
   BYPRODUCTS
      all.tests
   DEPENDS
      __t9_all_tests_o
      tinker9_acc
      tinker9_cu
      tinker9_cpp
      tinker9_f
      tinkerFToCPP
   COMMAND
      "${T9_ACC_COMPILER}"
      CUDA_HOME=${CUDA_DIR}
      -o all.tests
      $<TARGET_OBJECTS:__t9_all_tests_o>
      "-Wl,--start-group"
      $<TARGET_FILE:tinker9_acc>
      $<TARGET_FILE:tinker9_cu>
      $<TARGET_FILE:tinker9_cpp>
      $<TARGET_FILE:tinker9_f>
      "-Wl,--end-group"
      $<TARGET_FILE:tinkerFToCPP>
      "-L$<JOIN:${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES},;-L>"
      "-l$<JOIN:${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES},;-l>"
      -acc -Mcudalib=cufft,cublas
      ${__T9_DEVICE_LINK_FLAGS_${__T9_BUILD_TYPE}}
   COMMAND_EXPAND_LISTS
)
