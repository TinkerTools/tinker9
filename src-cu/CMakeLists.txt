cmake_minimum_required (VERSION 3.12)


enable_language (CUDA)


set (__T9_CUDA_INCPATH "${T9_INCPATH}")
list (REMOVE_ITEM __T9_CUDA_INCPATH "${PROJECT_SOURCE_DIR}/include/syntax/acc")


set (__T9_CUDA_SRC
compilers_cu.cu
echarge_cu.cu
echglj_cu.cu
echgtrn_cu.cu
edisp_cu.cu
ehal_cu.cu
elj_cu.cu
emplar_cu.cu
empole_chgpen_cu.cu
empole_cu.cu
epolar_chgpen_cu.cu
epolar_cu.cu
erepel_cu.cu
evalence_cu.cu
field_chgpen_cu.cu
field_cu.cu
induce_cu.cu
induce_donly_cu.cu
mathfunc_parallel_cu.cu
mdintg_cu.cu
mdpt_cu.cu
pcg_cu.cu
pcg_donly_cu.cu
pmestuf_cu.cu
rattle_cu.cu
spatial_cu.cu
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
