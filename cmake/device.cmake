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


add_executable (tinker9)
target_link_libraries (tinker9
   __t9_main_o
   -Wl,--start-group
   ${__T9_ACC_LIB_STR}
   tinker9_cu
   tinker9_cpp
   -Wl,--end-group
   tinker9_version
   tinker9_f
   tinkerFToCpp
   ${__T9_EXTRA_LINK_FLAGS}
)


add_executable (all.tests)
target_link_libraries (all.tests
   __t9_all_tests_o
   -Wl,--start-group
   ${__T9_ACC_LIB_STR}
   tinker9_cu
   tinker9_cpp
   -Wl,--end-group
   tinker9_version
   tinker9_f
   tinkerFToCpp
   ${__T9_EXTRA_LINK_FLAGS}
)


if (${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
   add_subdirectory ("src/objc")
   target_link_libraries (tinker9 tinker9_objc "-framework CoreFoundation -framework IOKit")
   target_link_libraries (all.tests tinker9_objc "-framework CoreFoundation -framework IOKit")
   foreach (var tinker9 all.tests)
      add_custom_command (TARGET "${var}" POST_BUILD
         COMMAND
            ${CMAKE_INSTALL_NAME_TOOL} -add_rpath "${CUDA_DIR}/lib" "${var}"
      )
   endforeach ()
endif ()
