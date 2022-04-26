## User Manual
add_custom_target (__t9_copy_manual_dir
   COMMAND
      ${CMAKE_COMMAND} -E copy_directory
         "${PROJECT_SOURCE_DIR}/doc/manual" "${CMAKE_BINARY_DIR}/manual"
   COMMAND
      ${CMAKE_COMMAND} -E make_directory
         "${CMAKE_BINARY_DIR}/manual/_build"
         "${CMAKE_BINARY_DIR}/manual/_static"
         "${CMAKE_BINARY_DIR}/manual/_template"
   COMMAND
      grep '^\!!/' "${PROJECT_SOURCE_DIR}/CMakeLists.txt" | sed -e 's|!!/ ||' | sed -e 's|!!/||'
         > "manual/m/install/buildwithcmake.rst"
   COMMAND
      ${CMAKE_COMMAND} -E copy
         "${CMAKE_BINARY_DIR}/manual/m/install/buildwithcmake.rst"
         "${PROJECT_SOURCE_DIR}/doc/manual/m/install/buildwithcmake.rst"
   BYPRODUCTS
      "${CMAKE_BINARY_DIR}/manual"
)


add_custom_target (man
   DEPENDS
      __t9_copy_manual_dir
   COMMAND
      make -C "${CMAKE_BINARY_DIR}/manual" html latexpdf
)


add_custom_target (html
   DEPENDS
      __t9_copy_manual_dir
   COMMAND
      make -C "${CMAKE_BINARY_DIR}/manual" html
)


add_custom_target (pdf
   DEPENDS
      __t9_copy_manual_dir
   COMMAND
      make -C "${CMAKE_BINARY_DIR}/manual" latexpdf
)


## Developer Guides
add_custom_target (doc
   COMMAND
      doxygen "${PROJECT_SOURCE_DIR}/doc/Doxyfile" ENV_GIT_HEAD="${__T9_GIT_HEAD}"
   BYPRODUCTS
      "${CMAKE_BINARY_DIR}/html"
)
