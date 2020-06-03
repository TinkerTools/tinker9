add_compile_definitions (TINKER_HOST)
list (APPEND proj_internal_inc_path "${PROJECT_SOURCE_DIR}/include/syntax/acc")
file (GLOB PLATFORM_CPP "${PROJECT_SOURCE_DIR}/src/host/*.cpp")
list (APPEND LIB_CPP ${PLATFORM_CPP})


add_library (tinkergpu0 STATIC ${LIB_F} ${LIB_ACC} ${LIB_CPP})
set_target_properties (tinkergpu0 PROPERTIES
   CXX_STANDARD 11
   CXX_EXTENSIONS OFF)
target_include_directories (tinkergpu0 SYSTEM PRIVATE ${comm_sys_inc_path})
target_include_directories (tinkergpu0 PRIVATE ${proj_internal_inc_path})


add_executable (tinker.gpu ${MAIN_CPP})
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
