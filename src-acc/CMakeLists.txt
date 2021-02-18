cmake_minimum_required (VERSION 3.12)


## Configured as
##
## cmake \
## -DCMAKE_CXX_COMPILER=... \
## -DCMAKE_BUILD_TYPE=... \
## -DCMAKE_INSTALL_PREFIX=... \
## -DT9_SYS_INCPATH=... \
## -DT9_INCPATH=... \
## -DT9_ACC_DEFS=... \
## -DT9_ACC_ONHOST=... \
## -DT9_ACC_CUDADIR=... \
## -DT9_ACC_CCLIST=...


project (Tinker9_OpenACC LANGUAGES CXX)


## libtinker9_acc.a
file (GLOB __T9_ACC_SRC "*_acc.cpp")
add_library (tinker9_acc STATIC "${__T9_ACC_SRC}")
set_target_properties (tinker9_acc PROPERTIES
   CXX_STANDARD
      ${T9_CPPSTD}
)
target_compile_definitions (tinker9_acc PRIVATE "${T9_ACC_DEFS}")
target_include_directories (tinker9_acc SYSTEM PRIVATE "${T9_SYS_INCPATH}")
target_include_directories (tinker9_acc PRIVATE "${T9_INCPATH}")
install (TARGETS tinker9_acc DESTINATION "${CMAKE_INSTALL_PREFIX}")



########################################################################
if (T9_ACC_ONHOST)
## host


if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
   string (APPEND CMAKE_CXX_FLAGS_DEBUG " -Wall -Wno-unknown-pragmas")
endif ()


set (__T9_IGNORE_ME "${T9_ACC_CCLIST}${T9_ACC_CUDADIR}" )


########################################################################
else ()
## device


## Compute Capability 60,70 -> ,cc60,cc70
set (__T9_ACC_CCLST4) # ""
foreach (var ${T9_ACC_CCLIST})
   string (APPEND __T9_ACC_CCLST4 ,cc${var})
endforeach () # ,cc60,cc70


target_compile_options (tinker9_acc PRIVATE
   CUDA_HOME=${T9_ACC_CUDADIR}
   -acc verystrict
   -Minfo=accel
)
## Debug add flag: -ta=tesla:lineinfo,cc60,cc70
## Release add flag: -ta=tesla:fastmath,cc60,cc70
string (APPEND CMAKE_CXX_FLAGS_DEBUG          " -ta=tesla:lineinfo${__T9_ACC_CCLST4}")
string (APPEND CMAKE_CXX_FLAGS_RELWITHDEBINFO " -ta=tesla:lineinfo,fastmath${__T9_ACC_CCLST4}")
string (APPEND CMAKE_CXX_FLAGS_RELEASE        " -ta=tesla:fastmath${__T9_ACC_CCLST4}")
string (APPEND CMAKE_CXX_FLAGS_MINSIZEREL     " -ta=tesla:fastmath${__T9_ACC_CCLST4}")


endif ()
