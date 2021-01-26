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
!!/ If not set, the building script will take a guess at the OpenACC compiler.
!!/ It will be set to the default C++ compiler for the CPU code. For instance,
!!/ command *(c)cmake [...]* will become *ACC=pgc++ (c)cmake [...]*.
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
!!/ ACC                      pgc++
!!/ opt                      release
!!/ host                     0
!!/ prec                     m
!!/ cuda_dir                 /usr/local/cuda
!!/ compute_capability       75
!!/ fftw_dir                 /usr/local
!!/ CMakeLists.txt Location  /home/tinker9
!!/ =======================  ===================
!!/
!!/ .. code-block:: bash
!!/
!!/    # use environmental variables
!!/    ACC=pgc++ \
!!/    opt=release host=0 prec=m \
!!/    cuda_dir=/usr/local/cuda compute_capability=75 \
!!/    fftw_dir=/usr/local \
!!/    cmake /home/tinker9
!!/
!!/    # use cmake -DOPTIONS
!!/    ACC=pgc++ cmake -S /home/tinker9 \
!!/    -DCMAKE_BUILD_TYPE=Release -DHOST=0 -DPREC=m \
!!/    -DCUDA_DIR=/usr/local/cuda -DCOMPUTE_CAPABILITY=75 \
!!/    -DFFTW_DIR=/usr/local
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
!!/ **-DCMAKE_INSTALL_PREFIX (no env var) = [NO DEFAULT VALUE]**
!!/
!!/ Install the executables under *${CMAKE_INSTALL_PREFIX}*. If this option is
!!/ not set, *make install* is configured not to install anything, which is
!!/ different from the default cmake behavior to install the program under */usr/local*.
!!/
]]
if (DEFINED CMAKE_INSTALL_PREFIX)
   set (__T9_INSTALL_DEFINED 1)
else ()
   set (__T9_INSTALL_DEFINED 0)
endif ()


#[[
!!/ **-DFFTW_DIR (fftw_dir) = ${CMAKE_BINARY_DIR}/fftw**
!!/
!!/ Top-level FFTW3 installation, under which
!!/ *include/fftw3.h* and *lib/libfftw3* static libraries are expected to be found.
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
!!/ `NVIDIA website. <https://developer.nvidia.com/cuda-gpus>`_
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
set (T9_DEFS)


## libtinker.a
ExternalProject_Add (src-libtinker
   SOURCE_DIR
      "${PROJECT_SOURCE_DIR}/src-libtinker"
   BINARY_DIR
      "${CMAKE_BINARY_DIR}/src-libtinker"
   PREFIX
      "${CMAKE_BINARY_DIR}/src-libtinker"
   CMAKE_ARGS
      "-DCMAKE_BUILD_TYPE=RelWithDebInfo"
      "-DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/src-libtinker"
)
add_library (t9_ltinker STATIC IMPORTED)
set_target_properties (t9_ltinker PROPERTIES
   IMPORTED_LOCATION
      "${CMAKE_BINARY_DIR}/src-libtinker/libtinker.a"
)
set_source_files_properties ("${CMAKE_BINARY_DIR}/src-libtinker/libtinker.a" PROPERTIES
   EXTERNAL_OBJECT
      TRUE
   GENERATED
      TRUE
)


## fftw
add_library (t9_lfftw STATIC IMPORTED)
set_target_properties (t9_lfftw PROPERTIES
   IMPORTED_LOCATION "${FFTW_DIR}/lib/libfftw3.a")
add_library (t9_lfftw_threads STATIC IMPORTED)
set_target_properties (t9_lfftw_threads PROPERTIES
   IMPORTED_LOCATION "${FFTW_DIR}/lib/libfftw3_threads.a"
)
if (EXISTS "${FFTW_DIR}/lib/libfftw3f.a")
   add_library (t9_lfftwf STATIC IMPORTED)
   set_target_properties (t9_lfftwf PROPERTIES
      IMPORTED_LOCATION "${FFTW_DIR}/lib/libfftw3f.a"
   )
endif ()
if (EXISTS "${FFTW_DIR}/lib/libfftw3f_threads.a")
   add_library (t9_lfftwf_threads STATIC IMPORTED)
   set_target_properties (t9_lfftwf_threads PROPERTIES
      IMPORTED_LOCATION "${FFTW_DIR}/lib/libfftw3f_threads.a"
   )
endif ()


## macro: precision
if (PREC STREQUAL "m")
   list (APPEND T9_DEFS TINKER_MIXED_PRECISION)
elseif (PREC STREQUAL "s")
   list (APPEND T9_DEFS TINKER_SINGLE_PRECISION)
elseif (PREC STREQUAL "d")
   list (APPEND T9_DEFS TINKER_DOUBLE_PRECISION)
endif ()


## macro: deterministic force
if (DETERMINISTIC_FORCE STREQUAL "")
   # No default behavior for empty value.
elseif (NOT DETERMINISTIC_FORCE)
   list (APPEND T9_DEFS TINKER_DETERMINISTIC_FORCE=0)
elseif (DETERMINISTIC_FORCE)
   list (APPEND T9_DEFS TINKER_DETERMINISTIC_FORCE=1)
endif ()


## macro: fortran compiler
if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
   list (APPEND T9_DEFS TINKER_GFORTRAN)
elseif (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
   list (APPEND T9_DEFS TINKER_IFORT)
endif ()


## macro: project directory
## TINKER9_DIR is needed in the source code.
list (APPEND T9_DEFS TINKER9_DIR=${PROJECT_SOURCE_DIR})


## openacc compiler
## set T9_ACC_COMPILER after project (), otherwise CMAKE_CXX_COMPILER is empty.
if (HOST)
   set (T9_ACC_COMPILER "${CMAKE_CXX_COMPILER}")
else ()
   if (DEFINED ENV{ACC})
      find_program (T9_ACC_COMPILER
         NAMES
            "$ENV{ACC}"
            "nvc++"
            "pgc++"
         REQUIRED
      )
   else ()
      find_program (T9_ACC_COMPILER
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
            __T9_GIT_SHA1
      )
      string (REPLACE "\"" "" __T9_GIT_SHA1 ${__T9_GIT_SHA1})
      string (REPLACE "\n" "" __T9_GIT_SHA1 ${__T9_GIT_SHA1})
      string (REPLACE "\r" "" __T9_GIT_SHA1 ${__T9_GIT_SHA1})


      execute_process (
         COMMAND
            "${GIT_EXECUTABLE}" log -1 --format=%cd
         WORKING_DIRECTORY
            "${PROJECT_SOURCE_DIR}"
         OUTPUT_VARIABLE
            __T9_GIT_DATE
      )
      string (REPLACE "\"" "" __T9_GIT_DATE ${__T9_GIT_DATE})
      string (REPLACE "\n" "" __T9_GIT_DATE ${__T9_GIT_DATE})
      string (REPLACE "\r" "" __T9_GIT_DATE ${__T9_GIT_DATE})


      execute_process (
         COMMAND
            "${GIT_EXECUTABLE}" log -1 --format="%h \(%cd\)"
         WORKING_DIRECTORY
            "${PROJECT_SOURCE_DIR}"
         OUTPUT_VARIABLE
            __T9_GIT_HEAD
      )
      string (REPLACE "\"" "" __T9_GIT_HEAD ${__T9_GIT_HEAD})
      string (REPLACE "\n" "" __T9_GIT_HEAD ${__T9_GIT_HEAD})
      string (REPLACE "\r" "" __T9_GIT_HEAD ${__T9_GIT_HEAD})
   endif ()
endif ()


########################################################################


## include path
set (T9_SYS_INCPATH
   "${FFTW_DIR}/include"
   "${PROJECT_SOURCE_DIR}/ext/ext")
set (T9_INCPATH "${PROJECT_SOURCE_DIR}/include")
if (HOST)
   list (APPEND T9_DEFS TINKER_HOST)
   list (APPEND T9_INCPATH "${PROJECT_SOURCE_DIR}/include/syntax/acc")
else ()
   list (APPEND T9_DEFS TINKER_CUDART)
   list (APPEND T9_INCPATH "${PROJECT_SOURCE_DIR}/include/syntax/acc")
   list (APPEND T9_SYS_INCPATH "${CUDA_DIR}/include")
   # remove "include/syntax/acc" from the include paths for src-cu
   # add "include/syntax/cu" to the include paths for src-cu
endif ()


## openacc
string (REPLACE ";" "|" T9_EP_ACC_MACRO "${T9_DEFS}")
string (REPLACE ";" "|" T9_EP_ACC_SYS_INCPATH "${T9_SYS_INCPATH}")
string (REPLACE ";" "|" T9_EP_ACC_INCPATH "${T9_INCPATH}")
string (REPLACE ";" "|" T9_ACC_CCLIST "${T9_CUCCLIST}")
ExternalProject_Add (src-acc
   SOURCE_DIR
      "${PROJECT_SOURCE_DIR}/src-acc"
   BINARY_DIR
      "${CMAKE_BINARY_DIR}/src-acc"
   PREFIX
      "${CMAKE_BINARY_DIR}/src-acc"
   BUILD_BYPRODUCTS
      "${CMAKE_BINARY_DIR}/src-acc"
   BUILD_ALWAYS
      1
   LIST_SEPARATOR
      |
   CMAKE_ARGS
      "-DCMAKE_CXX_COMPILER=${T9_ACC_COMPILER}"
      "-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}"
      "-DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/src-acc"
      "-DT9_SYS_INCPATH=${T9_EP_ACC_SYS_INCPATH}"
      "-DT9_INCPATH=${T9_EP_ACC_INCPATH}"
      "-DT9_ACC_DEFS=${T9_EP_ACC_MACRO}"
      "-DT9_ACC_ONHOST=${HOST}"
      "-DT9_ACC_CUDADIR=${CUDA_DIR}"
      "-DT9_ACC_CCLIST=${T9_ACC_CCLIST}"
)
set (__T9_ACC_STATIC_LIB_LOC "${CMAKE_BINARY_DIR}/src-acc/libtinker9_acc.a")
add_library (tinker9_EP_acc STATIC IMPORTED)
set_target_properties (tinker9_EP_acc PROPERTIES
   IMPORTED_LOCATION
      "${__T9_ACC_STATIC_LIB_LOC}"
)
set_source_files_properties ("${__T9_ACC_STATIC_LIB_LOC}" PROPERTIES
   EXTERNAL_OBJECT
      TRUE
   GENERATED
      TRUE
)


## test files
file (GLOB T9_TEST_SRC
   "src-main/main_all_tests.cpp"
   "src/test/*.cpp"
)
add_library (__t9_all_tests_o OBJECT ${T9_TEST_SRC})
target_compile_definitions (__t9_all_tests_o PRIVATE ${T9_DEFS})
set_target_properties (__t9_all_tests_o PROPERTIES
   CXX_STANDARD
      11
)
target_include_directories (__t9_all_tests_o SYSTEM PRIVATE "${T9_SYS_INCPATH}")
target_include_directories (__t9_all_tests_o PRIVATE
   "${T9_INCPATH}"
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
if (__T9_INSTALL_DEFINED)
   set (__T9_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}")
   if (HOST)
      set (__T9_INSTALL_DIR "${__T9_INSTALL_DIR}/cpu")
   else ()
      set (__T9_INSTALL_DIR "${__T9_INSTALL_DIR}/gpu")
   endif ()
   if (PREC STREQUAL "m")
      set (__T9_INSTALL_DIR "${__T9_INSTALL_DIR}-m")
   elseif (PREC STREQUAL "s")
      set (__T9_INSTALL_DIR "${__T9_INSTALL_DIR}-s")
   elseif (PREC STREQUAL "d")
      set (__T9_INSTALL_DIR "${__T9_INSTALL_DIR}-d")
   endif ()
   install (PROGRAMS tinker9
      DESTINATION
         "${__T9_INSTALL_DIR}"
   )
endif ()


## scripts
add_custom_target (__t9_copy_shell_exe)
foreach (var
   "${PROJECT_SOURCE_DIR}/shell/analyze9.sh"   
   "${PROJECT_SOURCE_DIR}/shell/bar9.sh"
   "${PROJECT_SOURCE_DIR}/shell/dynamic9.sh"
   "${PROJECT_SOURCE_DIR}/shell/info9.sh"
   "${PROJECT_SOURCE_DIR}/shell/minimize9.sh"
   "${PROJECT_SOURCE_DIR}/shell/testgrad9.sh"
   )
   get_filename_component (__t9_fname "${var}" NAME)
   add_custom_command (TARGET __t9_copy_shell_exe POST_BUILD
      COMMAND
         "${CMAKE_COMMAND}" -E copy "${var}" "${CMAKE_BINARY_DIR}"
      BYPRODUCTS
         "${CMAKE_BINARY_DIR}/${__t9_fname}"
   )
   if (__T9_INSTALL_DEFINED)
      install (PROGRAMS "${var}"
         DESTINATION
            "${__T9_INSTALL_DIR}"
      )
   endif ()
endforeach ()
add_dependencies (tinker9 __t9_copy_shell_exe)


## unit tests
add_custom_target (test
   DEPENDS
      all.tests
   COMMAND
      ./all.tests info
   COMMAND
      ./all.tests [ff],[util] -a --durations yes --order rand --rng-seed time
)


## documentations
include (docs) # docs.cmake
