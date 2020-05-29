#[============================================================[
!!/ Build Tinker GPU with CMake
!!/ ===========================
!!/
!!/ **Configure CMake**
!!/
]============================================================]


########################################################################


cmake_minimum_required (VERSION 3.10)


project (TinkerGPU LANGUAGES NONE)
enable_language (C)
enable_language (CXX)
enable_language (Fortran)
list (APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
set (CMAKE_BUILD_TYPE release CACHE STRING
   "Choose the type of build, options are: debug release")


########################################################################


set (tinker_dir "$ENV{HOME}/tinker/source" CACHE STRING "Path to libtinker.a")
add_library (LIBTINKER STATIC IMPORTED)
set_target_properties (LIBTINKER PROPERTIES
   IMPORTED_LOCATION "${tinker_dir}/libtinker.a")


set (fftw_dir /usr/local CACHE STRING "FFTW3 Directory")
add_library (LIBFFTW STATIC IMPORTED)
set_target_properties (LIBFFTW PROPERTIES
   IMPORTED_LOCATION "${fftw_dir}/lib/libfftw3.a")
add_library (LIBFFTW_THREADS STATIC IMPORTED)
set_target_properties (LIBFFTW_THREADS PROPERTIES
   IMPORTED_LOCATION "${fftw_dir}/lib/libfftw3_threads.a")
add_library (LIBFFTWF STATIC IMPORTED)
set_target_properties (LIBFFTWF PROPERTIES
   IMPORTED_LOCATION "${fftw_dir}/lib/libfftw3f.a")
add_library (LIBFFTWF_THREADS STATIC IMPORTED)
set_target_properties (LIBFFTWF_THREADS PROPERTIES
   IMPORTED_LOCATION "${fftw_dir}/lib/libfftw3f_threads.a")


set (host 0 CACHE STRING "0 (GPU) / 1 (CPU)")


set (prec m CACHE STRING "d (DOUBLE) / s (SINGLE) / m (MIXED)")
if (${prec} STREQUAL m)
   add_compile_definitions (TINKER_MIXED_PRECISION)
elseif (${prec} STREQUAL s)
   add_compile_definitions (TINKER_SINGLE_PRECISION)
elseif (${prec} STREQUAL d)
   add_compile_definitions (TINKER_DOUBLE_PRECISION)
endif ()


set (deterministic_force "" CACHE STRING
   "Deterministic Force: [NOT SET] / 0 (OFF) / 1 (ON)"
)
if (deterministic_force STREQUAL "")
   # No default behavior for empty value.
elseif (${deterministic_force} STREQUAL 1)
   add_compile_definitions (TINKER_DETERMINISTIC_FORCE=1)
elseif (${deterministic_force} STREQUAL 0)
   add_compile_definitions (TINKER_DETERMINISTIC_FORCE=0)
endif ()


set (compute_capability "60,70" CACHE STRING "CUDA Compute Capability")
set (cuda_dir /usr/local/cuda CACHE STRING "CUDA Directory")


if (${CMAKE_Fortran_COMPILER_ID} STREQUAL GNU)
   add_compile_definitions (TINKER_GFORTRAN)
elseif (${CMAKE_Fortran_COMPILER_ID} STREQUAL Intel)
   add_compile_definitions (TINKER_IFORT)
endif ()


find_package (Git)
if (GIT_FOUND)
   if (EXISTS "${PROJECT_SOURCE_DIR}/.git")
      execute_process (
         COMMAND
            "${GIT_EXECUTABLE}" log -1 --format="%h"
         WORKING_DIRECTORY
            "${PROJECT_SOURCE_DIR}"
         OUTPUT_VARIABLE
            STRING_SHA1_GIT
      )
      string (REPLACE "\"" "" STRING_SHA1_GIT ${STRING_SHA1_GIT})
      string (REPLACE "\n" "" STRING_SHA1_GIT ${STRING_SHA1_GIT})
      string (REPLACE "\r" "" STRING_SHA1_GIT ${STRING_SHA1_GIT})
      add_compile_definitions (TINKER_GPU_GIT_SHORT_HASH=${STRING_SHA1_GIT})


      execute_process (
         COMMAND
            "${GIT_EXECUTABLE}" log -1 --format="%h \(%cd\)"
         WORKING_DIRECTORY
            "${PROJECT_SOURCE_DIR}"
         OUTPUT_VARIABLE
            STRING_GIT_HEAD
      )
      string (REPLACE "\"" "" STRING_GIT_HEAD ${STRING_GIT_HEAD})
      string (REPLACE "\n" "" STRING_GIT_HEAD ${STRING_GIT_HEAD})
      string (REPLACE "\r" "" STRING_GIT_HEAD ${STRING_GIT_HEAD})
   endif ()
endif ()


########################################################################


## include path
set (comm_sys_inc_path
   "${fftw_dir}/include"
   "${PROJECT_SOURCE_DIR}/ext/ext")
set (proj_internal_inc_path "${PROJECT_SOURCE_DIR}/include")


## main file
file (GLOB MAIN_CPP "${PROJECT_SOURCE_DIR}/src/main*.cpp")
## fortran file(s)
file (GLOB LIB_F "${PROJECT_SOURCE_DIR}/src/*.f")
## acc files
file (GLOB LIB_ACC "${PROJECT_SOURCE_DIR}/src/*_acc.cpp")
## cpp files
file (GLOB LIB_CPP "${PROJECT_SOURCE_DIR}/src/*.cpp")
list (REMOVE_ITEM LIB_CPP ${MAIN_CPP} ${LIB_ACC})
## cu files
file (GLOB LIB_CU "${PROJECT_SOURCE_DIR}/src/*_cu.cu")


## tinker.gpu
if (host STREQUAL 1)
   include (host) # host.cmake
elseif (host STREQUAL 0)
   include (device) # device.cmake
endif ()


## scripts
add_custom_target (SHELL_EXE COMMAND
   "${CMAKE_COMMAND}" -E copy "${PROJECT_SOURCE_DIR}/shell/*.gpu" "${CMAKE_BINARY_DIR}")
add_dependencies (tinker.gpu SHELL_EXE)


## unit tests
include (test) # test.cmake


## documentations
include (docs) # docs.cmake
