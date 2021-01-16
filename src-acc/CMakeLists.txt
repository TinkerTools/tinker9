cmake_minimum_required (VERSION 3.11)


project (Tinker9_OpenACC LANGUAGES CXX)


## libtinker9_acc.a
file (GLOB TINKER9_OPENACC_SRC "*_acc.cpp")
add_library (tinker9_acc STATIC "${TINKER9_OPENACC_SRC}")
set_target_properties (tinker9_acc PROPERTIES
   CXX_STANDARD
      11
)
target_compile_definitions (tinker9_acc PRIVATE "${TINKER9_MACRO}")
target_include_directories (tinker9_acc SYSTEM PRIVATE "${TINKER9_SYS_INC_PATH}")
target_include_directories (tinker9_acc PRIVATE "${TINKER9_INTERNAL_INC_PATH}")
install (TARGETS tinker9_acc DESTINATION "${CMAKE_INSTALL_PREFIX}")



########################################################################
if (T9_ACC_ONHOST)
## host


set (T9_IGNORE_ME "${T9_ACC_CCLIST}${T9_ACC_CUDADIR}" )


########################################################################
else ()
## device


## Compute Capability 60,70 -> ,cc60,cc70
set (T9_ACC_CCLST4) # ""
foreach (var ${T9_ACC_CCLIST})
   string (APPEND T9_ACC_CCLST4 ,cc${var})
endforeach () # ,cc60,cc70


target_compile_options (tinker9_acc PRIVATE
   CUDA_HOME=${T9_ACC_CUDADIR}
   -acc verystrict
   -Minfo=accel
)
## Debug add flag: -ta=tesla:lineinfo,cc60,cc70
## Release add flag: -ta=tesla:fastmath,cc60,cc70
target_compile_options (tinker9_acc PRIVATE
   $<$<CONFIG:DEBUG>:-ta=tesla:lineinfo${T9_ACC_CCLST4}>
   $<$<CONFIG:RELEASE>:-ta=tesla:fastmath${T9_ACC_CCLST4}>
)


endif ()
