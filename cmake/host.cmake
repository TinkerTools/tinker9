list (APPEND macro_defs TINKER_HOST)
list (APPEND TINKER9_INTERNAL_INC_PATH "${PROJECT_SOURCE_DIR}/include/syntax/acc")
file (GLOB PLATFORM_CPP "${PROJECT_SOURCE_DIR}/src/host/*.cpp")
list (APPEND LIB_CPP ${PLATFORM_CPP})


add_library (tinker9_f OBJECT ${LIB_F})
add_library (tinker9_cpp OBJECT ${LIB_ACC} ${LIB_CPP})
target_compile_definitions (tinker9_cpp PRIVATE ${macro_defs})
set_target_properties (tinker9_cpp PROPERTIES
   CXX_STANDARD 11
   CXX_EXTENSIONS OFF)
target_include_directories (tinker9_cpp SYSTEM PRIVATE "${TINKER9_SYS_INC_PATH}")
target_include_directories (tinker9_cpp PRIVATE "${TINKER9_INTERNAL_INC_PATH}")
add_library (tinker9_host STATIC
   $<TARGET_OBJECTS:tinker9_f>
   $<TARGET_OBJECTS:tinker9_cpp>
)


add_executable (tinker9 ${MAIN_CPP})
target_compile_definitions (tinker9 PRIVATE ${macro_defs})
set (EXT_LIBS pthread LIBTINKER LIBFFTW LIBFFTW_THREADS)
if (PREC STREQUAL "m" OR PREC STREQUAL "s")
   list (APPEND EXT_LIBS LIBFFTWF LIBFFTWF_THREADS)
endif ()
set_target_properties (tinker9 PROPERTIES
   CXX_STANDARD 11
   CXX_EXTENSIONS OFF)
target_include_directories (tinker9 SYSTEM PRIVATE "${TINKER9_SYS_INC_PATH}")
target_include_directories (tinker9 PRIVATE "${TINKER9_INTERNAL_INC_PATH}")
target_link_libraries (tinker9 tinker9_host ${EXT_LIBS})
