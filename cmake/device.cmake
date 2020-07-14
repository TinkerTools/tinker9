list (APPEND macro_defs TINKER_CUDART)
list (APPEND comm_sys_inc_path "${CUDA_DIR}/include")
file (GLOB PLATFORM_CPP "${PROJECT_SOURCE_DIR}/src/cudart/*.cpp")
list (APPEND LIB_CPP ${PLATFORM_CPP})


## CUDA
add_library (tinkergpu_cu OBJECT ${LIB_CU})
target_compile_definitions (tinkergpu_cu PRIVATE ${macro_defs})
enable_language (CUDA)
set_target_properties (tinkergpu_cu PROPERTIES
   CUDA_STANDARD 11
   CUDA_EXTENSIONS OFF
)
target_include_directories (tinkergpu_cu SYSTEM PRIVATE ${comm_sys_inc_path})
target_include_directories (tinkergpu_cu PRIVATE
   ${proj_internal_inc_path}
   "${PROJECT_SOURCE_DIR}/include/syntax/cu"
)
## Compute Capability 60,70 ->
## -gencode arch=compute_60,code=sm_60 -gencode arch=compute_70,code=sm_70
set (CCLIST2 ${COMPUTE_CAPABILITY}) # 60,70
string (REPLACE "," ";" CCLIST2 ${CCLIST2}) # 60;70
set (GENCODE_CC) # ""
foreach (var ${CCLIST2})
   # available since cmake 3.12
   target_compile_options (tinkergpu_cu PRIVATE
     "SHELL:-gencode arch=compute_${var},code=sm_${var}"
   )
endforeach () # -gencode arch=compute_60,code=sm_60 -gencode arch=compute_70,code=sm_70
## Debug add flag: -lineinfo
## Release add flag: --use_fast_math
target_compile_options (tinkergpu_cu PRIVATE
   $<$<CONFIG:DEBUG>:-lineinfo>
   $<$<CONFIG:RELEASE>:--use_fast_math -lineinfo>
)


## OpenACC
add_library (tinkergpu_acc OBJECT ${LIB_ACC})
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "PGI")
   if (CMAKE_CXX11_STANDARD_COMPILE_OPTION)
      # remove the -A flag coded in CMake/Modules/Compiler/PGI-CXX.cmake
      set (CMAKE_CXX11_STANDARD_COMPILE_OPTION -std=c++11)
   endif ()
else ()
   message (FATAL_ERROR "CXX was not set to PGI compiler.")
endif ()
target_compile_definitions (tinkergpu_acc PRIVATE ${macro_defs})
set_target_properties (tinkergpu_acc PROPERTIES
   CXX_STANDARD 11
   CXX_EXTENSIONS OFF
)
target_include_directories (tinkergpu_acc SYSTEM PRIVATE ${comm_sys_inc_path})
target_include_directories (tinkergpu_acc PRIVATE
   ${proj_internal_inc_path}
   "${PROJECT_SOURCE_DIR}/include/syntax/acc"
)
## Compute Capability 60,70 -> ,cc60,cc70
set (CCLIST4) # ""
foreach (var ${CCLIST2})
   string (APPEND CCLIST4 ,cc${var})
endforeach () # ,cc60,cc70
target_compile_options (tinkergpu_acc PRIVATE
   CUDA_HOME=${CUDA_DIR} -acc verystrict -Minfo=accel
)
## Debug add flag: -ta=tesla:lineinfo,cc60,cc70
## Release add flag: -ta=tesla:fastmath,cc60,cc70
target_compile_options (tinkergpu_acc PRIVATE
   $<$<CONFIG:DEBUG>:-ta=tesla:lineinfo${CCLIST4}>
   $<$<CONFIG:RELEASE>:-ta=tesla:fastmath${CCLIST4}>
)


add_library (tinkergpu_f OBJECT ${LIB_F})
add_library (tinkergpu_cpp OBJECT ${LIB_CPP})
target_compile_definitions (tinkergpu_cpp PRIVATE ${macro_defs})
set_target_properties (tinkergpu_cpp PROPERTIES
   CXX_STANDARD 11
   CXX_EXTENSIONS OFF
)
target_include_directories (tinkergpu_cpp SYSTEM PRIVATE ${comm_sys_inc_path})
target_include_directories (tinkergpu_cpp PRIVATE
   ${proj_internal_inc_path}
   "${PROJECT_SOURCE_DIR}/include/syntax/acc"
)


add_library (tinkergpu_main OBJECT ${MAIN_CPP})
target_compile_definitions (tinkergpu_main PRIVATE ${macro_defs})
set_target_properties (tinkergpu_main PROPERTIES
   CXX_STANDARD 11
   CXX_EXTENSIONS OFF
)
target_include_directories (tinkergpu_main SYSTEM PRIVATE ${comm_sys_inc_path})
target_include_directories (tinkergpu_main PRIVATE
   ${proj_internal_inc_path}
   "${PROJECT_SOURCE_DIR}/include/syntax/acc"
)


separate_arguments (LIST_CXX_FLAGS_DEBUG NATIVE_COMMAND ${CMAKE_CXX_FLAGS_DEBUG})
separate_arguments (LIST_CXX_FLAGS_RELEASE NATIVE_COMMAND ${CMAKE_CXX_FLAGS_RELEASE})


########################################################################


add_custom_target (tinker.gpu ALL
   DEPENDS
      tinkergpu_main
      tinkergpu_acc
      tinkergpu_cu
      tinkergpu_f
      tinkergpu_cpp
      LIBTINKER
      LIBFFTW
      LIBFFTW_THREADS
   COMMAND
      ${CMAKE_CXX_COMPILER} CUDA_HOME=${CUDA_DIR}
      "$<$<CONFIG:DEBUG>:${LIST_CXX_FLAGS_DEBUG}>"
      "$<$<CONFIG:RELEASE>:${LIST_CXX_FLAGS_RELEASE}>"
      -o tinker.gpu
      $<TARGET_OBJECTS:tinkergpu_main>
      $<TARGET_OBJECTS:tinkergpu_acc>
      $<TARGET_OBJECTS:tinkergpu_cu>
      $<TARGET_OBJECTS:tinkergpu_f>
      $<TARGET_OBJECTS:tinkergpu_cpp>
      $<TARGET_FILE:LIBTINKER>
      $<TARGET_FILE:LIBFFTW>
      $<TARGET_FILE:LIBFFTW_THREADS>
      "-L$<JOIN:${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES},;-L>"
      "-l$<JOIN:${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES},;-l>"
      -L${CUDA_DIR}/lib64/stubs -lnvidia-ml
      -acc -Mcudalib=cufft,cublas
      $<$<CONFIG:DEBUG>:-ta=tesla:lineinfo${CCLIST4}>
      $<$<CONFIG:RELEASE>:-ta=tesla:fastmath${CCLIST4}>
   COMMAND_EXPAND_LISTS
)
