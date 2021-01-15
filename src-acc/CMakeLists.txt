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
if (HOST)
## host


########################################################################
else ()
## device


## Compute Capability 60,70 -> ,cc60,cc70
set (TINKER9_CCLST4) # ""
foreach (var ${TINKER9_CCLIST})
   string (APPEND TINKER9_CCLST4 ,cc${var})
endforeach () # ,cc60,cc70


target_compile_options (tinker9_acc PRIVATE
   CUDA_HOME=${CUDA_DIR}
   -acc verystrict
   -Minfo=accel
)
## Debug add flag: -ta=tesla:lineinfo,cc60,cc70
## Release add flag: -ta=tesla:fastmath,cc60,cc70
target_compile_options (tinker9_acc PRIVATE
   $<$<CONFIG:DEBUG>:-ta=tesla:lineinfo${TINKER9_CCLST4}>
   $<$<CONFIG:RELEASE>:-ta=tesla:fastmath${TINKER9_CCLST4}>
)


endif ()
