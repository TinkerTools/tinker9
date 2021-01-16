cmake_minimum_required (VERSION 3.11)


#[[
!!/ Build Tinker9 with CMake
!!/ ========================
!!/
!!/ Configure CMake
!!/ ---------------
!!/ You can skip this section if you are familar with CMake.
!!/
!!/ Suppose the current working directory is */home/tinker9* and we
!!/ want to create a build directory called *build-cmake* in
!!/ */home/tinker9*. We can do *mkdir build-cmake* then *cd build-cmake*.
!!/ Because the top-level CMakeLists.txt file is in the parent directory,
!!/ if there was nothing else to configure, command *cmake ..* would generate
!!/ the Makefile. The alternative way is to specify the build and source
!!/ directories to CMake, e.g.,
!!/
!!/ .. code-block:: bash
!!/
!!/    cmake -B /home/tinker9/build-cmake -S /home/tinker9
!!/
!!/ Some CMake installations also provide a command line gui *ccmake* and a
!!/ simple gui program *cmake-gui* that can replace *cmake* in the commands
!!/ above.
!!/
!!/ Configure Compilers
!!/ -------------------
!!/ Set *CXX=...*, *CUDACXX=...*, and *FC=...* to specify the non-default C++,
!!/ CUDA, and Fortran compilers, respectively. These environmental variables
!!/ are supported by *cmake*.
!!/
!!/ This cmake script checks a custom environmental variable *ACC=...*
!!/ *only* for the GPU code.
!!/ If not set, the building script will take a guess at the OpenACC compilers.
!!/ It will be set to the default C++ compiler for the CPU code.
!!/
!!/ In general the only situation where you have to configure the C++
!!/ compiler in command line, regardless of using *cmake* or *ccmake*,
!!/ is to compile the GPU OpenACC code with *pgc++*. The command
!!/ *(c)cmake [...]* will become *CXX=pgc++ (c)cmake [...]*.
!!/
!!/ Configure Tinker9
!!/ -----------------
!!/ The following options are passed to CMake program with their default
!!/ values (if there is one). **-D** is prefixed to the options. CMake provides
!!/ two standard ways to let users to customize the values:
!!/
!!/ - Change their values interactively in the *ccmake* command line gui;
!!/ - Pass the new value to CMake via command line arguments
!!/   *cmake -DOPTION=NewValue*.
!!/
!!/ In addition to these two canonical methods, default value can also be set
!!/ by its corresponding environmental variable, documented as **(env)** here.
!!/ Note that there is no **-D** prefix for the environmental variables.
!!/
!!/ Here are two equivalent examples to have Tinker9 configured as follows
!!/
!!/ =======================  ===================
!!/ Item                     Value
!!/ =======================  ===================
!!/ CXX                      pgc++
!!/ opt                      release
!!/ host                     0
!!/ prec                     m
!!/ cuda_dir                 /usr/local/cuda
!!/ compute_capability       75
!!/ tinker_dir               /home/tinker/source
!!/ fftw_dir                 /usr/local
!!/ CMakeLists.txt Location  /home/tinker9
!!/ =======================  ===================
!!/
!!/ .. code-block:: bash
!!/
!!/    # use environmental variables
!!/    CXX=pgc++ \
!!/    opt=release host=0 prec=m \
!!/    cuda_dir=/usr/local/cuda compute_capability=75 \
!!/    tinker_dir=/home/tinker/source fftw_dir=/usr/local \
!!/    cmake /home/tinker9
!!/
!!/    # use cmake -DOPTIONS
!!/    CXX=pgc++ cmake -S /home/tinker9 \
!!/    -DCMAKE_BUILD_TYPE=Release -DHOST=0 -DPREC=m \
!!/    -DCUDA_DIR=/usr/local/cuda -DCOMPUTE_CAPABILITY=75 \
!!/    -DTINKER_DIR=/home/tinker/source -DFFTW_DIR=/usr/local
!!/
!!/
]]


########################################################################


#[[
!!/ **-DCMAKE_BUILD_TYPE (opt) = release**
!!/
!!/ Build type is case insensitive and can either be *release* or *debug*.
!!/
]]
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
   if (DEFINED ENV{opt})
      set (optValue $ENV{opt})
   else ()
      set (optValue Release)
   endif ()
   set (CMAKE_BUILD_TYPE ${optValue} CACHE STRING
      "Choose the type of build. Options: Debug / Release"
      FORCE
   )
endif ()


#[[
!!/ **-DTINKER_DIR (tinker_dir) = ${CMAKE_BINARY_DIR}/tinker/source**
!!/
!!/ Directory in which user compiled *libtinker.a*.
!!/
]]
if (DEFINED ENV{tinker_dir})
   get_filename_component (tinkerDirValue "$ENV{tinker_dir}" ABSOLUTE
      BASE_DIR
         "${CMAKE_BINARY_DIR}"
   )
else ()
   set (tinkerDirValue "${CMAKE_BINARY_DIR}/tinker/source")
endif ()
set (TINKER_DIR ${tinkerDirValue} CACHE PATH "Path to libtinker.a")


#[[
!!/ **-DFFTW_DIR (fftw_dir) = ${CMAKE_BINARY_DIR}/fftw**
!!/
!!/ Top-level FFTW3 installation, under which
!!/ *include/fftw3.h* and *lib/libfftw3* are expected to be found.
!!/
]]
if (DEFINED ENV{fftw_dir})
   get_filename_component (fftwDirValue "$ENV{fftw_dir}" ABSOLUTE
      BASE_DIR
         "${CMAKE_BINARY_DIR}"
   )
else ()
   set (fftwDirValue "${CMAKE_BINARY_DIR}/fftw")
endif ()
set (FFTW_DIR ${fftwDirValue} CACHE PATH "FFTW3 Directory")


#[[
!!/ **-DHOST (host) = 1**
!!/
!!/ Flag to compile to GPU (with value 0 or OFF) or CPU (with value 1 or ON)
!!/ version.
!!/
]]
if (DEFINED ENV{host})
   set (hostValue $ENV{host})
else ()
   set (hostValue 1)
endif ()
set (HOST ${hostValue} CACHE BOOL "0 or OFF GPU) / 1 or ON (CPU)")


#[[
!!/ **-DPREC (prec) = d**
!!/
!!/ Precision of the floating-point numbers. With flag *d*, all of the
!!/ floating-point numbers are treated as real*8/double values,
!!/ or real*4/single values if with flag *s*. Mixed precision flag *m* will
!!/ use real*4 or real*8 numbers in different places. Note that this flag will
!!/ not change the precision of the variables hard-coded as *float* or *double*
!!/ types.
!!/
]]
if (DEFINED ENV{prec})
   set (precValue $ENV{prec})
else ()
   set (precValue d)
endif ()
set (PREC ${precValue} CACHE STRING "d (DOUBLE) / s (SINGLE) / m (MIXED)")


#[[
!!/ **-DDETERMINISTIC_FORCE (deterministic_force) = [NO DEFAULT VALUE]**
!!/
!!/ Flag to use deterministic force. There is no default value for this flag.
!!/ This feature will be implicitly enabled by mixed and single precisions, but
!!/ can be explicitly disabled by setting the flag to 0,
!!/ and can be explicitly enabled by value 1.
!!/
!!/ In general, evaluating energy, forces etc. twice, we don't expect to get
!!/ two identical answers, but we may not care as much because the difference
!!/ is usually negligible. (See
!!/ `Why is cos(x) != cos(y)? <https://isocpp.org/wiki/faq/newbie#floating-point-arith2>`_)
!!/ Whereas in MD, two simulations with the same initial configurations can
!!/ easily diverge due to the accumulated difference. If, for whatever reason,
!!/ you are willing to elongate the process of the inevitable divergence at the
!!/ cost of slightly slower simulation speed, a more "deterministic" force
!!/ (using fixed-point arithmetic) can help.
!!/
]]
if (DEFINED ENV{deterministic_force})
   set (detForceValue $ENV{deterministic_force})
else ()
   set (detForceValue)
endif ()
set (DETERMINISTIC_FORCE ${detForceValue} CACHE STRING
   "Deterministic Force: [NO DEFAULT] / 0 (OFF) / 1 (ON)"
)


#[[
!!/ **-DCOMPUTE_CAPABILITY (compute_capability) = 60,70**
!!/
!!/ CUDA compute capability (multiplied by 10) of GPU.
!!/ Valid values (noninclusive) are 35, 50, 60, 70, 75 etc., and can be
!!/ comma-separated, e.g. 35,60.
!!/ Multiple compute capabilites will increase the size of executables.
!!/
!!/ The full list of compute capabilities can be found on the
!!/ `Nvidia website. <https://developer.nvidia.com/cuda-gpus>`_
!!/
]]
if (DEFINED ENV{compute_capability})
   set (compCapValue $ENV{compute_capability})
else ()
   set (compCapValue 60,70)
endif ()
set (COMPUTE_CAPABILITY ${compCapValue} CACHE STRING "CUDA Compute Capability")
## Compute Capability 60,70 -> 60;70
string (REPLACE "," ";" T9_CUCCLIST ${COMPUTE_CAPABILITY}) # 60;70


#[[
!!/ **-DCUDA_DIR (cuda_dir) = /usr/local/cuda**
!!/
!!/ Top-level CUDA installation directory, under which directories *include*,
!!/ *lib* or *lib64* can be found.
!!/ This option will supersede the CUDA installation identified by the official
!!/ *CUDACXX* environmental variable.
!!/
!!/ Sometimes the PGI compiler and the NVCC compiler are not "compatible." For
!!/ instance, although PGI 19.4 supports CUDA 9.2, 10.0, 10.1, but the default
!!/ CUDA version configured in PGI 19.4 may be 9.2 and the external NVCC version
!!/ is 10.1. One solution is to pass *CUDA_HOME=${cuda_dir}* to the PGI
!!/ compiler, in which case, **cuda_dir** should be set to
!!/ */usr/local/cuda-10.1*.
!!/
]]
if (DEFINED ENV{cuda_dir})
   get_filename_component (cudaDirValue "$ENV{cuda_dir}" ABSOLUTE
      BASE_DIR
         "${CMAKE_BINARY_DIR}"
   )
elseif (DEFINED ENV{CUDACXX})
   find_program (cudaCCValue $ENV{CUDACXX})
   get_filename_component (cudaCCValue "${cudaCCValue}" DIRECTORY)
else ()
   set (cudaDirValue /usr/local/cuda)
endif ()
set (CUDA_DIR ${cudaDirValue} CACHE PATH "CUDA Directory")
## cuda compiler
## set CMAKE_CUDA_COMPILER before project ()   
if (NOT HOST)
   set (CMAKE_CUDA_COMPILER "${CUDA_DIR}/bin/nvcc")
endif ()


#[[
!!/ Make Tinker9
!!/ ------------
!!/ The following Makefile targets will be generated by CMake.
!!/ Run *make -j* for the default target(s) and *make TARGET(S) -j* for others.
!!/
!!/ **tinker9**
!!/
!!/ Compile and link the *tinker9* executable.
!!/
!!/ **all.tests**
!!/
!!/ Compile and link the *all.tests* executable.
!!/
!!/ **default**
!!/
!!/ Make two targets: *tinker9* and *all.tests* executables.
!!/
!!/ **all**
!!/
!!/ Same as the default target.
!!/
!!/ **test**
!!/
!!/ Run unit tests in a random order. Exit on the first error.
!!/
!!/ **man**
!!/
!!/ Generate user manual.
!!/
!!/ **doc**
!!/
!!/ Generate developer guides.
!!/
]]


########################################################################


include (ExternalProject)


project (Tinker9 LANGUAGES NONE)
enable_language (C)
enable_language (CXX)
enable_language (Fortran)


## macro: an empty list
set (TINKER9_DEFS)


## libtinker.a
add_library (LIBTINKER STATIC IMPORTED)
set_target_properties (LIBTINKER PROPERTIES
   IMPORTED_LOCATION
      "${TINKER_DIR}/libtinker.a"
)


## fftw
add_library (LIBFFTW STATIC IMPORTED)
set_target_properties (LIBFFTW PROPERTIES
   IMPORTED_LOCATION "${FFTW_DIR}/lib/libfftw3.a")
add_library (LIBFFTW_THREADS STATIC IMPORTED)
set_target_properties (LIBFFTW_THREADS PROPERTIES
   IMPORTED_LOCATION "${FFTW_DIR}/lib/libfftw3_threads.a"
)
if (EXISTS "${FFTW_DIR}/lib/libfftw3f.a")
   add_library (LIBFFTWF STATIC IMPORTED)
   set_target_properties (LIBFFTWF PROPERTIES
      IMPORTED_LOCATION "${FFTW_DIR}/lib/libfftw3f.a"
   )
endif ()
if (EXISTS "${FFTW_DIR}/lib/libfftw3f_threads.a")
   add_library (LIBFFTWF_THREADS STATIC IMPORTED)
   set_target_properties (LIBFFTWF_THREADS PROPERTIES
      IMPORTED_LOCATION "${FFTW_DIR}/lib/libfftw3f_threads.a"
   )
endif ()


## macro: precision
if (PREC STREQUAL "m")
   list (APPEND TINKER9_DEFS TINKER_MIXED_PRECISION)
elseif (PREC STREQUAL "s")
   list (APPEND TINKER9_DEFS TINKER_SINGLE_PRECISION)
elseif (PREC STREQUAL "d")
   list (APPEND TINKER9_DEFS TINKER_DOUBLE_PRECISION)
endif ()


## macro: deterministic force
if (DETERMINISTIC_FORCE STREQUAL "")
   # No default behavior for empty value.
elseif (NOT DETERMINISTIC_FORCE)
   list (APPEND TINKER9_DEFS TINKER_DETERMINISTIC_FORCE=0)
elseif (DETERMINISTIC_FORCE)
   list (APPEND TINKER9_DEFS TINKER_DETERMINISTIC_FORCE=1)
endif ()


## macro: fortran compiler
if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
   list (APPEND TINKER9_DEFS TINKER_GFORTRAN)
elseif (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
   list (APPEND TINKER9_DEFS TINKER_IFORT)
endif ()


## macro: src directory
list (APPEND TINKER9_DEFS TINKER9_DIR=${CMAKE_SOURCE_DIR})


## openacc compiler
## set TINKER9_ACC_COMPILER after project (), otherwise CMAKE_CXX_COMPILER is empty.
if (HOST)
   set (TINKER9_ACC_COMPILER "${CMAKE_CXX_COMPILER}")
else ()
   if (DEFINED ENV{ACC})
      find_program (TINKER9_ACC_COMPILER
         NAMES
            "$ENV{ACC}"
            "nvc++"
            "pgc++"
         REQUIRED
      )
   else ()
      find_program (TINKER9_ACC_COMPILER
         NAMES
            "nvc++"
            "pgc++"
         REQUIRED
      )
   endif ()
endif ()


########################################################################


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
      list (APPEND TINKER9_DEFS TINKER9_GIT_SHORT_HASH=${STRING_SHA1_GIT})


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
set (TINKER9_SYS_INC_PATH
   "${FFTW_DIR}/include"
   "${PROJECT_SOURCE_DIR}/ext/ext")
set (TINKER9_INTERNAL_INC_PATH "${PROJECT_SOURCE_DIR}/include")
if (HOST)
   list (APPEND TINKER9_DEFS TINKER_HOST)
   list (APPEND TINKER9_INTERNAL_INC_PATH "${PROJECT_SOURCE_DIR}/include/syntax/acc")
else ()
   list (APPEND TINKER9_DEFS TINKER_CUDART)
   list (APPEND TINKER9_INTERNAL_INC_PATH "${PROJECT_SOURCE_DIR}/include/syntax/acc")
   list (APPEND TINKER9_SYS_INC_PATH "${CUDA_DIR}/include")
   # remove "include/syntax/acc" from the include paths for src-cu
   # add "include/syntax/cu" to the include paths for src-cu
endif ()


## openacc
string (REPLACE ";" "|" TINKER9_EP_OPENACC_MACRO "${TINKER9_DEFS}")
string (REPLACE ";" "|" TINKER9_EP_OPENACC_SYS_INC_PATH "${TINKER9_SYS_INC_PATH}")
string (REPLACE ";" "|" TINKER9_EP_OPENACC_INTERNAL_INC_PATH "${TINKER9_INTERNAL_INC_PATH}")
string (REPLACE ";" "|" T9_ACCCCLIST "${T9_CUCCLIST}")
ExternalProject_Add (cm9a
   SOURCE_DIR
      "${PROJECT_SOURCE_DIR}/src-acc"
   BINARY_DIR
      "${CMAKE_BINARY_DIR}/cm9a"
   PREFIX
      "${CMAKE_BINARY_DIR}/cm9a"
   LIST_SEPARATOR
      |
   CMAKE_ARGS
      "-DCMAKE_CXX_COMPILER=${TINKER9_ACC_COMPILER}"
      "-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}"
      "-DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/cm9a"
      "-DTINKER9_MACRO=${TINKER9_EP_OPENACC_MACRO}"
      "-DTINKER9_SYS_INC_PATH=${TINKER9_EP_OPENACC_SYS_INC_PATH}"
      "-DTINKER9_INTERNAL_INC_PATH=${TINKER9_EP_OPENACC_INTERNAL_INC_PATH}"
      "-DT9_ACCONHOST=${HOST}"
      "-DT9_ACC_CUDADIR=${CUDA_DIR}"
      "-DT9_ACCCCLIST=${T9_ACCCCLIST}"
)
set (TINKER9_ACC_STATIC_LIB_LOC "${CMAKE_BINARY_DIR}/cm9a/libtinker9_acc.a")
add_library (tinker9_EP_acc STATIC IMPORTED)
set_target_properties (tinker9_EP_acc PROPERTIES
   IMPORTED_LOCATION
      "${TINKER9_ACC_STATIC_LIB_LOC}"
)
set_source_files_properties ("${TINKER9_ACC_STATIC_LIB_LOC}" PROPERTIES
   EXTERNAL_OBJECT
      TRUE
   GENERATED
      TRUE
)


## test files
file (GLOB TINKER9_TEST_SRC
   "src-main/main_all_tests.cpp"
   "src/test/*.cpp"
)
add_library (all_tests_o OBJECT ${TINKER9_TEST_SRC})
target_compile_definitions (all_tests_o PRIVATE ${TINKER9_DEFS})
set_target_properties (all_tests_o PROPERTIES
   CXX_STANDARD
      11
)
target_include_directories (all_tests_o SYSTEM PRIVATE "${TINKER9_SYS_INC_PATH}")
target_include_directories (all_tests_o PRIVATE
   "${TINKER9_INTERNAL_INC_PATH}"
   "${PROJECT_SOURCE_DIR}/src/test"
)


## executables
## tinker9 and all.tests
list (APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
if (HOST)
   add_subdirectory (src)
   add_subdirectory (src-f)
   include (host) # host.cmake
else ()
   add_subdirectory (src-cu)
   add_subdirectory (src)
   add_subdirectory (src-f)
   include (device) # device.cmake
endif ()


## scripts
add_custom_target (SHELL_EXE COMMAND
   "${CMAKE_COMMAND}" -E copy "${PROJECT_SOURCE_DIR}/shell/*9.sh" "${CMAKE_BINARY_DIR}")
add_dependencies (tinker9 SHELL_EXE)


## unit tests
add_custom_target (test
   COMMAND
      ./all.tests info
   COMMAND
      ./all.tests [ff],[util] -a --durations yes --order rand --rng-seed time
   DEPENDS
      all.tests
)


## documentations
include (docs) # docs.cmake


## more files to clean
set_property (DIRECTORY
   PROPERTY
      ADDITIONAL_MAKE_CLEAN_FILES
         "${CMAKE_BINARY_DIR}/README.md"
         "${CMAKE_BINARY_DIR}/doc"
         "${CMAKE_BINARY_DIR}/html"
         "${CMAKE_BINARY_DIR}/include"
         "${CMAKE_BINARY_DIR}/latex"
         "${CMAKE_BINARY_DIR}/manual"
)
