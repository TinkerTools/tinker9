# printf "set (__T9_CU_SRC\n%s\n)\n" "$(LC_COLLATE=C ls -1a *.cu)" > cmakesrc.txt
cmake_minimum_required (VERSION 3.12)


enable_language (CUDA)
set (CMAKE_CUDA_STANDARD ${T9_CPPSTD})
set (CMAKE_CUDA_STANDARD_REQUIRED ON)
set (CMAKE_CUDA_EXTENSIONS OFF)


include (cmakesrc.txt)
add_library (tinker9_cu OBJECT "${__T9_CU_SRC}")
target_include_directories (tinker9_cu PRIVATE "${PROJECT_SOURCE_DIR}/include/syntax/cu")
target_link_libraries (tinker9_cu PRIVATE tinker9_cpp)


## -gencode arch=compute_60,code=sm_60 -gencode arch=compute_70,code=sm_70
foreach (var ${T9_CUCCLIST})
   target_compile_options (tinker9_cu PRIVATE
      "SHELL:-gencode arch=compute_${var},code=sm_${var}"
   )
endforeach ()
## Debug add flag: -lineinfo
## Release add flag: --use_fast_math
target_compile_options (tinker9_cu PRIVATE
   "$<$<CONFIG:DEBUG>:-lineinfo>"
   "$<$<CONFIG:RELWITHDEBINFO>:-lineinfo;--use_fast_math>"
   "$<$<CONFIG:RELEASE>:--use_fast_math>"
   "$<$<CONFIG:MINSIZEREL>:--use_fast_math>"
)
target_compile_options (tinker9_cu PRIVATE
   "$<$<CONFIG:DEBUG>:--compiler-options;-Wall>"
   "$<$<CONFIG:DEBUG>:--compiler-options;-Wno-unknown-pragmas>"
)
