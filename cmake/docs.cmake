## User Manual
add_custom_target (man
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
   COMMAND
      make -C "${CMAKE_BINARY_DIR}/manual" html latexpdf
)


## Developer Guides
add_custom_target (doc
   COMMAND
      ${CMAKE_COMMAND} -E copy
         "${PROJECT_SOURCE_DIR}/README.md" "${CMAKE_BINARY_DIR}"
   COMMAND
      ${CMAKE_COMMAND} -E copy_directory
         "${PROJECT_SOURCE_DIR}/doc" "${CMAKE_BINARY_DIR}/doc"
   COMMAND
      ${CMAKE_COMMAND} -E copy_directory
         "${PROJECT_SOURCE_DIR}/include" "${CMAKE_BINARY_DIR}/include"
   COMMAND
      doxygen "${PROJECT_SOURCE_DIR}/Doxyfile" ENV_GIT_HEAD="${STRING_GIT_HEAD}"
)
