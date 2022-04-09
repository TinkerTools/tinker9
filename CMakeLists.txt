cmake_minimum_required (VERSION 3.15)


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
!!/ *ACC* will also be used as the C++ compiler. The value of *CXX* (if set)
!!/ will be neglected.
!!/
!!/ Configure Tinker9
!!/ -----------------
!!/ The following options are passed to CMake program with their default
!!/ values (if there is one). **-D** is prefixed to the options. CMake provides
!!/ two standard ways to let users customize the values:
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
!!/ opt                      release
!!/ host                     0
!!/ prec                     m
!!/ cuda_dir                 /usr/local/cuda
!!/ compute_capability       75
!!/ CMakeLists.txt Location  /home/tinker9
!!/ =======================  ===================
!!/
!!/ .. code-block:: bash
!!/
!!/    # use environmental variables
!!/    opt=release host=0 prec=m \
!!/    cuda_dir=/usr/local/cuda compute_capability=75 \
!!/    cmake /home/tinker9
!!/
!!/    # use cmake -DOPTIONS
!!/    cmake /home/tinker9 \
!!/    -DCMAKE_BUILD_TYPE=Release -DHOST=0 -DPREC=m \
!!/    -DCUDA_DIR=/usr/local/cuda -DCOMPUTE_CAPABILITY=75
!!/
!!/
]]


########################################################################


#[[
!!/ **-DCMAKE_BUILD_TYPE (opt) = Release**
!!/
!!/ Standard *CMAKE_BUILD_TYPE* option. Build type is case insensitive and
!!/ can be *Debug*, *Release*, *RelWithDebInfo* (release with debug info),
!!/ and *MinSizeRel* (minimum size release).
!!/
]]
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
   if (DEFINED ENV{opt})
      set (optValue $ENV{opt})
   else ()
      set (optValue Release)
   endif ()
   string (TOLOWER ${optValue} optValue)
   if (optValue STREQUAL debug)
      set (optValue Debug)
   elseif (optValue STREQUAL release)
      set (optValue Release)
   elseif (optValue STREQUAL relwithdebinfo)
      set (optValue RelWithDebInfo)
   elseif (optValue STREQUAL minsizerel)
      set (optValue MinSizeRel)
   endif ()
   set (CMAKE_BUILD_TYPE ${optValue} CACHE STRING
      "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel ..."
      FORCE
   )
   set_property (CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
      Debug Release RelWithDebInfo MinSizeRel
   )
endif ()
string (TOUPPER ${CMAKE_BUILD_TYPE} __T9_BUILD_TYPE)


#[[
!!/ **-DCMAKE_INSTALL_PREFIX (install) = [NO DEFAULT VALUE]**
!!/
!!/ Install the executables under *${CMAKE_INSTALL_PREFIX}*. If this option is
!!/ not set, *make install* is configured not to install anything, which is
!!/ different from the default cmake behavior to install the program under */usr/local*.
!!/
]]
if (DEFINED CMAKE_INSTALL_PREFIX)
   set (__T9_INSTALL_DEFINED 1)
elseif (DEFINED ENV{install})
   set (installValue "$ENV{install}")
   set (CMAKE_INSTALL_PREFIX ${installValue} CACHE STRING
      "Install path prefix, prepended onto install directories."
      FORCE
   )
   if (installValue STREQUAL "")
      set (__T9_INSTALL_DEFINED 0)
   elseif ()
      set (__T9_INSTALL_DEFINED 1)
   endif ()
else ()
   set (__T9_INSTALL_DEFINED 0)
endif ()


#[[
!!/ **-DSTD (std) = 11**
!!/
!!/ C++ syntax standard. The source code is c++11-compliant, and should have no
!!/ problems compiled with c++14. If set to *14* here, users should make sure
!!/ the compilers are c++14-compliant.
!!/
]]
if (DEFINED ENV{std})
   set (stdValue $ENV{std})
else ()
   set (stdValue 11)
endif ()
set (STD ${stdValue} CACHE STRING "C++ standard: 11 and 14.")
set_property (CACHE STD PROPERTY STRINGS 11 14)
set (T9_CPPSTD ${STD})


#[[
!!/ **-DPREC (prec) = mixed**
!!/
!!/ Precision of the floating-point numbers. With flag *double*/*d*, all of the
!!/ floating-point numbers are treated as real\*8/double values,
!!/ or real\*4/single values if with flag *single*/*s*. Mixed precision flag *mixed*/*m* will
!!/ use real\*4 or real\*8 numbers in different places. Note that this flag will
!!/ not change the precision of the variables hard-coded as *float* or *double*
!!/ types.
!!/
]]
if (DEFINED ENV{prec})
   set (precValue $ENV{prec})
else ()
   set (precValue mixed)
endif ()
string (TOLOWER ${precValue} precValue)
set (PREC ${precValue} CACHE STRING "Floating-Point Precision: double/d, single/s, mixed/m")
set_property (CACHE PREC PROPERTY STRINGS double single mixed d s m)
if (PREC STREQUAL "double")
   set (PREC "d")
elseif (PREC STREQUAL "single")
   set (PREC "s")
elseif (PREC STREQUAL "mixed")
   set (PREC "m")
endif ()


#[[
!!/ **-DDETERMINISTIC_FORCE (deterministic_force) = AUTO**
!!/
!!/ Flag to use deterministic force.
!!/ This feature will be implicitly enabled by mixed and single precisions, but
!!/ can be explicitly disabled by setting the flag to *OFF* (or 0),
!!/ and can be explicitly enabled by value *ON* (or 1).
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
   set (detForceValue AUTO)
endif ()
string (TOUPPER ${detForceValue} detForceValue)
set (DETERMINISTIC_FORCE ${detForceValue} CACHE STRING
   "Deterministic Force: AUTO, OFF, ON"
)
set_property (CACHE DETERMINISTIC_FORCE PROPERTY STRINGS AUTO "OFF" "ON")


#[[
!!/ **-DHOST (host) = OFF**
!!/
!!/ Flag to compile to GPU (with value 0 or OFF) or CPU (with value 1 or ON)
!!/ version.
!!/
]]
if (DEFINED ENV{host})
   set (hostValue $ENV{host})
else ()
   set (hostValue OFF)
endif ()
set (HOST ${hostValue} CACHE BOOL "Build Executables for GPU Platform (OFF) or CPU Platform (ON)")


#[[
!!/ **-DGPU_LANG (gpu_lang) = OPENACC** (experimental, in progress)
!!/
!!/ If set to *CUDA*, the GPU code will only use the cuda source files.
!!/ And the program will crash at runtime if it falls into an OpenACC code path.
!!/
]]
if (DEFINED ENV{gpu_lang})
   set (gpuLangValue $ENV{gpu_lang})
else ()
   set (gpuLangValue openacc)
endif ()
string (TOUPPER ${gpuLangValue} gpuLangValue)
set (GPU_LANG ${gpuLangValue} CACHE STRING
   "GPU Programming Language (experimental): OPENACC, CUDA"
)


#[[
!!/ **-DCOMPUTE_CAPABILITY (compute_capability) = 60,70**
!!/
!!/ GPU code only.
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
set (COMPUTE_CAPABILITY ${compCapValue} CACHE STRING
   "[GPU ONLY] CUDA Compute Capability Multiplied by 10 (Comma-Separated)"
)
## Compute Capability 60,70 -> 60;70
string (REPLACE "," ";" T9_CUCCLIST ${COMPUTE_CAPABILITY}) # 60;70


#[[
!!/ **-DCUDA_DIR (cuda_dir) = /usr/local/cuda**
!!/
!!/ Nvidia GPU code only.
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
   find_program (T9_CUDA_COMPILER $ENV{CUDACXX})
   get_filename_component (cudaDirValue "${T9_CUDA_COMPILER}" DIRECTORY) # /usr/local/cuda/bin
   get_filename_component (cudaDirValue "${cudaDirValue}" DIRECTORY)     # /usr/local/cuda
else ()
   set (cudaDirValue /usr/local/cuda)
endif ()
set (CUDA_DIR ${cudaDirValue} CACHE PATH "[GPU ONLY] CUDA Directory")
## cuda compiler
## set CMAKE_CUDA_COMPILER before project ()
if (HOST)
else ()
   set (CMAKE_CUDA_COMPILER "${CUDA_DIR}/bin/nvcc")
endif ()


#[[
!!/ **-DFFTW_DIR (fftw_dir) = ${CMAKE_BINARY_DIR}/fftw**
!!/
!!/ CPU code only.
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
set (FFTW_DIR ${fftwDirValue} CACHE PATH "[CPU ONLY] FFTW3 Directory")


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


## openacc compiler
if (HOST)
elseif (GPU_LANG STREQUAL "OPENACC")
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
   # set this line before project()
   set (CMAKE_CXX_COMPILER "${T9_ACC_COMPILER}")
   # 11   unrecognized preprocessing directive
   # 177  variable was declared but never referenced
   string (APPEND CMAKE_CXX_FLAGS " --diag_suppress=11,177")
   # Compute Capability 60,70 -> ,cc60,cc70
   set (__T9_ACC_CCLST4) # ""
   foreach (var ${T9_CUCCLIST})
      string (APPEND __T9_ACC_CCLST4 ",cc${var}")
   endforeach () # ,cc60,cc70
endif ()


########################################################################


project (Tinker9 LANGUAGES NONE)
if (GPU_LANG STREQUAL "CUDA")
   enable_language (CUDA)
endif ()
enable_language (CXX)
enable_language (Fortran)
set (CMAKE_CXX_STANDARD ${T9_CPPSTD})
set (CMAKE_CXX_STANDARD_REQUIRED ON)
set (CMAKE_CXX_EXTENSIONS OFF)


if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
   string (APPEND CMAKE_CXX_FLAGS_DEBUG " -Wall -Wno-unknown-pragmas")
endif ()


## fftw
if (PREC STREQUAL "d")
   add_library (t9_lfftw STATIC IMPORTED)
   set_target_properties (t9_lfftw PROPERTIES
      IMPORTED_LOCATION "${FFTW_DIR}/lib/libfftw3.a"
   )
   add_library (t9_lfftw_threads STATIC IMPORTED)
   set_target_properties (t9_lfftw_threads PROPERTIES
      IMPORTED_LOCATION "${FFTW_DIR}/lib/libfftw3_threads.a"
   )
elseif (PREC STREQUAL "m" OR PREC STREQUAL "s")
   add_library (t9_lfftw STATIC IMPORTED)
   set_target_properties (t9_lfftw PROPERTIES
      IMPORTED_LOCATION "${FFTW_DIR}/lib/libfftw3f.a"
   )
   add_library (t9_lfftw_threads STATIC IMPORTED)
   set_target_properties (t9_lfftw_threads PROPERTIES
      IMPORTED_LOCATION "${FFTW_DIR}/lib/libfftw3f_threads.a"
   )
endif ()


find_package (Git)
if (GIT_FOUND)
   execute_process (
      COMMAND
         "${GIT_EXECUTABLE}" submodule update --init
      WORKING_DIRECTORY
         "${PROJECT_SOURCE_DIR}"
   )


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


## macro: an empty list
set (T9_DEFS)
## macro: precision
if (PREC STREQUAL "m")
   list (APPEND T9_DEFS TINKER_MIXED_PRECISION)
elseif (PREC STREQUAL "s")
   list (APPEND T9_DEFS TINKER_SINGLE_PRECISION)
elseif (PREC STREQUAL "d")
   list (APPEND T9_DEFS TINKER_DOUBLE_PRECISION)
endif ()
## macro: HOST
if (HOST)
   list (APPEND T9_DEFS TINKER_HOST)
else ()
   list (APPEND T9_DEFS TINKER_CUDART)
   if (GPU_LANG STREQUAL "CUDA")
      list (APPEND T9_DEFS TINKER_GPULANG_CUDA)
   elseif (GPU_LANG STREQUAL "OPENACC")
      list (APPEND T9_DEFS TINKER_GPULANG_OPENACC)
   endif ()
endif ()
## macro: deterministic force
if (DETERMINISTIC_FORCE STREQUAL "AUTO")
   # No default behavior for empty value.
elseif (NOT DETERMINISTIC_FORCE)
   list (APPEND T9_DEFS TINKER_DETERMINISTIC_FORCE=0)
elseif (DETERMINISTIC_FORCE)
   list (APPEND T9_DEFS TINKER_DETERMINISTIC_FORCE=1)
endif ()
## macro: project directory
## TINKER9_DIR is needed in the source code.
list (APPEND T9_DEFS TINKER9_DIR=${PROJECT_SOURCE_DIR})


## libs
add_subdirectory (ext/interface)
add_subdirectory (src/f)
add_subdirectory (src)
add_subdirectory (src/acc)
if (HOST)
else ()
   add_subdirectory (src/cu)
endif ()
add_subdirectory (test)


## executables
## tinker9 and all.tests
list (APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
if (HOST)
   include (host) # host.cmake
else ()
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
   install (PROGRAMS "${CMAKE_BINARY_DIR}/tinker9"
      DESTINATION
         "${__T9_INSTALL_DIR}"
   )
endif ()


## symlinks
add_custom_target (__t9_symlink)
foreach (var
   analyze9
   bar9
   dynamic9
   info9
   minimize9
   testgrad9
   )
   add_custom_command (TARGET __t9_symlink POST_BUILD
      COMMAND
         "${CMAKE_COMMAND}" -E create_symlink tinker9 "${var}"
   )
   if (__T9_INSTALL_DEFINED)
      install (PROGRAMS "${CMAKE_BINARY_DIR}/${var}"
         DESTINATION
            "${__T9_INSTALL_DIR}"
      )
   endif ()
endforeach ()
add_dependencies (tinker9 __t9_symlink)


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
