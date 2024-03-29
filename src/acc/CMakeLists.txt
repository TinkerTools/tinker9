# printf "set (__T9_ACC_SRC\n%s\n)\n" "$(LC_COLLATE=C ls -1a *.cpp {amoeba,aplus,hippo}/*.cpp)" > cmakesrc.txt
cmake_minimum_required (VERSION 3.18)


## libtinker9_acc.a
include (cmakesrc.txt)
add_library (tinker9_acc STATIC "${__T9_ACC_SRC}")
target_link_libraries (tinker9_acc PRIVATE __t9_intf tinkerFToCpp)
set_target_properties (tinker9_acc PROPERTIES CXX_STANDARD ${T9_CPPSTD})


########################################################################


if (HOST)
elseif (GPU_LANG STREQUAL "OPENACC")
   string (APPEND CMAKE_CXX_FLAGS " CUDA_HOME=${CUDA_DIR}")
   target_compile_options (tinker9_acc PUBLIC
      -acc verystrict
      -Minfo=accel
   )
   ## Debug add flag: -ta=tesla:lineinfo,cc60,cc70
   ## Release add flag: -ta=tesla:fastmath,cc60,cc70
   target_compile_options (tinker9_acc PUBLIC
      "$<$<CONFIG:DEBUG>:-ta=tesla:lineinfo${__T9_ACC_CCLST4}>"
      "$<$<CONFIG:RELWITHDEBINFO>:-ta=tesla:lineinfo,fastmath${__T9_ACC_CCLST4}>"
      "$<$<CONFIG:RELEASE>:-ta=tesla:fastmath${__T9_ACC_CCLST4}>"
      "$<$<CONFIG:MINSIZEREL>:-ta=tesla:fastmath${__T9_ACC_CCLST4}>"
   )
endif ()
