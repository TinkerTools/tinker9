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
!!/ Set *CXX=...* and *FC=...* to specify the non-default C++ and Fortran
!!/ compilers, respectively.
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


project (Tinker9 LANGUAGES NONE)
enable_language (C)
enable_language (CXX)
set (macro_defs)
enable_language (Fortran)
list (APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
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


########################################################################


#[[
!!/ **-DTINKER_DIR (tinker_dir) = ${CMAKE_BINARY_DIR}/tinker/source**
!!/
!!/ Directory in which user compiled *libtinker.a*.
!!/
]]
if (DEFINED ENV{tinker_dir})
   get_filename_component (tinkerDirValue "$ENV{tinker_dir}"
      ABSOLUTE BASE_DIR "${CMAKE_BINARY_DIR}"
   )
else ()
   set (tinkerDirValue "${CMAKE_BINARY_DIR}/tinker/source")
endif ()
set (TINKER_DIR ${tinkerDirValue} CACHE PATH "Path to libtinker.a")
##
add_library (LIBTINKER STATIC IMPORTED)
set_target_properties (LIBTINKER PROPERTIES
   IMPORTED_LOCATION "${TINKER_DIR}/libtinker.a"
)



#[[
!!/ **-DFFTW_DIR (fftw_dir) = ${CMAKE_BINARY_DIR}/fftw**
!!/
!!/ Top-level FFTW3 installation, under which
!!/ *include/fftw3.h* and *lib/libfftw3* are expected to be found.
!!/
]]
if (DEFINED ENV{fftw_dir})
   get_filename_component (fftwDirValue "$ENV{fftw_dir}"
      ABSOLUTE BASE_DIR "${CMAKE_BINARY_DIR}"
   )
else ()
   set (fftwDirValue "${CMAKE_BINARY_DIR}/fftw")
endif ()
set (FFTW_DIR ${fftwDirValue} CACHE PATH "FFTW3 Directory")
##
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
##
if (PREC STREQUAL "m")
   list (APPEND macro_defs TINKER_MIXED_PRECISION)
elseif (PREC STREQUAL "s")
   list (APPEND macro_defs TINKER_SINGLE_PRECISION)
elseif (PREC STREQUAL "d")
   list (APPEND macro_defs TINKER_DOUBLE_PRECISION)
endif ()


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
##
if (DETERMINISTIC_FORCE STREQUAL "")
   # No default behavior for empty value.
elseif (NOT DETERMINISTIC_FORCE)
   list (APPEND macro_defs TINKER_DETERMINISTIC_FORCE=0)
elseif (DETERMINISTIC_FORCE)
   list (APPEND macro_defs TINKER_DETERMINISTIC_FORCE=1)
endif ()


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
#[[
!!/ **-DCUDA_DIR (cuda_dir) = /usr/local/cuda**
!!/
!!/ Top-level CUDA installation directory, under which directories *include*,
!!/ *lib* or *lib64* can be found.
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
   get_filename_component (cudaDirValue "$ENV{cuda_dir}"
      ABSOLUTE BASE_DIR "${CMAKE_BINARY_DIR}"
   )
else ()
   set (cudaDirValue /usr/local/cuda)
endif ()
set (CUDA_DIR ${cudaDirValue} CACHE PATH "CUDA Directory")


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
!!/ Generate user's manual.
!!/
!!/ **doc**
!!/
!!/ Generate developer's manual.
!!/
]]


if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
   list (APPEND macro_defs TINKER_GFORTRAN)
elseif (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
   list (APPEND macro_defs TINKER_IFORT)
endif ()


list (APPEND macro_defs TINKER9_DIR=${CMAKE_SOURCE_DIR})


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
      list (APPEND macro_defs TINKER9_GIT_SHORT_HASH=${STRING_SHA1_GIT})


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
   "${FFTW_DIR}/include"
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


## tinker9
if (HOST)
   include (host) # host.cmake
else ()
   include (device) # device.cmake
endif ()


## scripts
add_custom_target (SHELL_EXE COMMAND
   "${CMAKE_COMMAND}" -E copy "${PROJECT_SOURCE_DIR}/shell/*9.sh" "${CMAKE_BINARY_DIR}")
add_dependencies (tinker9 SHELL_EXE)


## unit tests
include (test) # test.cmake


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
