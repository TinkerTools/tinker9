if (GPU_LANG STREQUAL "OPENACC")
   set (__T9_EXTRA_LINK_FLAGS
      -acc
      -Mcudalib=cufft,cublas
      "$<$<CONFIG:DEBUG>:-ta=tesla:lineinfo${__T9_ACC_CCLST4}>"
      "$<$<CONFIG:RELWITHDEBINFO>:-ta=tesla:lineinfo,fastmath${__T9_ACC_CCLST4}>"
      "$<$<CONFIG:RELEASE>:-ta=tesla:fastmath${__T9_ACC_CCLST4}>"
      "$<$<CONFIG:MINSIZEREL>:-ta=tesla:fastmath${__T9_ACC_CCLST4}>"
   )
   set (__T9_ACC_LIB_STR tinker9_acc)
elseif (GPU_LANG STREQUAL "CUDA")
   set (__T9_EXTRA_LINK_FLAGS "")
   # set (__T9_ACC_LIB_STR tinker9_acc)
endif ()


add_executable (tinker9 src/main.cc)
target_link_libraries (tinker9
   ${__T9_ACC_LIB_STR}
   tinker9_cu
   tinker9_cpp
   tinker9_version
   tinker9_f
   tinkerFToCpp
   ${__T9_EXTRA_LINK_FLAGS}
)


add_executable (all.tests)
target_link_libraries (all.tests
   __t9_all_tests_o
   ${__T9_ACC_LIB_STR}
   tinker9_cu
   tinker9_cpp
   tinker9_version
   tinker9_f
   tinkerFToCpp
   ${__T9_EXTRA_LINK_FLAGS}
)
