cmake_minimum_required (VERSION 3.12)


enable_language (CUDA)


set (__T9_CUDA_INCPATH "${T9_INCPATH}")
list (REMOVE_ITEM __T9_CUDA_INCPATH "${PROJECT_SOURCE_DIR}/include/syntax/acc")


set (__T9_CUDA_SRC
compilers.cu
echarge.cu
echglj.cu
echgtrn.cu
edisp.cu
ehal.cu
elj.cu
emplar.cu
empole_chgpen.cu
empole.cu
epolar_chgpen.cu
epolar.cu
erepel.cu
evalence.cu
field_chgpen.cu
field.cu
induce.cu
induce_donly.cu
lpiston.cu
mathfunc_parallel.cu
mdintg.cu
mdpt.cu
pcg.cu
pcg_donly.cu
pmestuf.cu
rattle.cu
spatial.cu
)
add_library (tinker9_cu STATIC "${__T9_CUDA_SRC}")
set_target_properties (tinker9_cu PROPERTIES
   CUDA_STANDARD
      ${T9_CPPSTD}
)
target_compile_definitions (tinker9_cu PRIVATE ${T9_DEFS})
target_include_directories (tinker9_cu SYSTEM PRIVATE "${T9_SYS_INCPATH}")
target_include_directories (tinker9_cu PRIVATE
   "${__T9_CUDA_INCPATH}"
   "${PROJECT_SOURCE_DIR}/include/syntax/cu"
)
## -gencode arch=compute_60,code=sm_60 -gencode arch=compute_70,code=sm_70
foreach (var ${T9_CUCCLIST})
   target_compile_options (tinker9_cu PRIVATE
      "SHELL:-gencode arch=compute_${var},code=sm_${var}"
   )
endforeach ()
## Debug add flag: -lineinfo
## Release add flag: --use_fast_math
string (APPEND CMAKE_CUDA_FLAGS_DEBUG          " -lineinfo")
string (APPEND CMAKE_CUDA_FLAGS_RELWITHDEBINFO " -lineinfo --use_fast_math")
string (APPEND CMAKE_CUDA_FLAGS_RELEASE        " --use_fast_math")
string (APPEND CMAKE_CUDA_FLAGS_MINSIZEREL     " --use_fast_math")


string (APPEND CMAKE_CUDA_FLAGS_DEBUG " --compiler-options -Wall")
string (APPEND CMAKE_CUDA_FLAGS_DEBUG " --compiler-options -Wno-unknown-pragmas")
