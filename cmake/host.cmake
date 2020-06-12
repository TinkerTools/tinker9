list (APPEND macro_defs TINKER_HOST)
list (APPEND proj_internal_inc_path "${PROJECT_SOURCE_DIR}/include/syntax/acc")
file (GLOB PLATFORM_CPP "${PROJECT_SOURCE_DIR}/src/host/*.cpp")
list (APPEND LIB_CPP ${PLATFORM_CPP})


add_library (tinkergpu_f OBJECT ${LIB_F})
add_library (tinkergpu_cpp OBJECT ${LIB_ACC} ${LIB_CPP})
target_compile_definitions (tinkergpu_cpp PRIVATE ${macro_defs})
set_target_properties (tinkergpu_cpp PROPERTIES
   CXX_STANDARD 11
   CXX_EXTENSIONS OFF)
target_include_directories (tinkergpu_cpp SYSTEM PRIVATE ${comm_sys_inc_path})
target_include_directories (tinkergpu_cpp PRIVATE ${proj_internal_inc_path})
add_library (tinkergpu0 STATIC
   $<TARGET_OBJECTS:tinkergpu_f>
   $<TARGET_OBJECTS:tinkergpu_cpp>
)


add_executable (tinker.gpu ${MAIN_CPP})
target_compile_definitions (tinker.gpu PRIVATE ${macro_defs})
set (EXT_LIBS pthread LIBTINKER LIBFFTW LIBFFTW_THREADS)
if (${PREC} STREQUAL m OR ${PREC} STREQUAL s)
   list (APPEND EXT_LIBS LIBFFTWF LIBFFTWF_THREADS)
endif ()
set_target_properties (tinker.gpu PROPERTIES
   CXX_STANDARD 11
   CXX_EXTENSIONS OFF)
target_include_directories (tinker.gpu SYSTEM PRIVATE ${comm_sys_inc_path})
target_include_directories (tinker.gpu PRIVATE ${proj_internal_inc_path})
target_link_libraries (tinker.gpu tinkergpu0 ${EXT_LIBS})
