# printf "set (__T9_CPP_SRC\n%s\n)\n" "$(LC_COLLATE=C ls -1a *.cpp {ff,ffamoeba,ffhippo,ffpchg,math,md,tool}/*.cpp)" > cmakesrc.txt
cmake_minimum_required (VERSION 3.12)


include (cmakesrc.txt)
if (HOST)
   list (APPEND __T9_CPP_SRC
host/darray.cpp
host/execq.cpp
host/fft.cpp
host/gpucard.cpp
host/pmestream.cpp
   )
else ()
   list (APPEND __T9_CPP_SRC
cudart/darray.cpp
cudart/error.cpp
cudart/execq.cpp
cudart/fft.cpp
cudart/gpucard.cpp
cudart/pmestream.cpp
cudart/thrustcache.cpp
   )
endif ()
add_library (tinker9_cpp OBJECT "${__T9_CPP_SRC}")


add_library (tinker9_version OBJECT
   xinfo.cc
   promo.cc
)
target_link_libraries (tinker9_version PUBLIC tinker9_cpp)


# macro
target_compile_definitions (tinker9_cpp PUBLIC ${T9_DEFS})
if (__T9_GIT_SHA1)
   target_compile_definitions (tinker9_version PRIVATE
      TINKER9_GIT_SHORT_HASH=${__T9_GIT_SHA1})
endif()
if (__T9_GIT_DATE)
   target_compile_definitions (tinker9_version PRIVATE
      TINKER9_GIT_DATE="${__T9_GIT_DATE}")
endif ()
# include
if (HOST)
else ()
   target_include_directories (tinker9_cpp PUBLIC "${CUDA_DIR}/include")
endif ()
target_include_directories (tinker9_cpp     PUBLIC  "${PROJECT_SOURCE_DIR}/include")
target_include_directories (tinker9_cpp     PRIVATE "${PROJECT_SOURCE_DIR}/include/syntax/acc")
target_include_directories (tinker9_version PRIVATE "${PROJECT_SOURCE_DIR}/include/syntax/acc")
# lib
target_link_libraries (tinker9_cpp PUBLIC tinkerFToCpp)
# FFTW
if (HOST)
   target_include_directories (tinker9_cpp SYSTEM PUBLIC "${FFTW_DIR}/include")
   target_link_libraries (tinker9_cpp PUBLIC t9_lfftw t9_lfftw_threads)
else ()
endif ()
