cmake_minimum_required (VERSION 3.11)


set (CMAKE_CUDA_COMPILER "${CUDA_DIR}/bin/nvcc")


enable_language (CUDA)


set (TINKER9_CUDA_INTERNAL_INC_PATH "${TINKER9_INTERNAL_INC_PATH}")
list (REMOVE_ITEM TINKER9_CUDA_INTERNAL_INC_PATH "${PROJECT_SOURCE_DIR}/include/syntax/acc")


file (GLOB TINKER9_CUDA_SRC "*_cu.cu")
add_library (tinker9_cu STATIC "${TINKER9_CUDA_SRC}")
set_target_properties (tinker9_cu PROPERTIES
   CUDA_STANDARD
      11
)
target_compile_definitions (tinker9_cu PRIVATE ${TINKER9_DEFS})
target_include_directories (tinker9_cu SYSTEM PRIVATE "${TINKER9_SYS_INC_PATH}")
target_include_directories (tinker9_cu PRIVATE
   "${TINKER9_CUDA_INTERNAL_INC_PATH}"
   "${PROJECT_SOURCE_DIR}/include/syntax/cu"
)
foreach (var ${TINKER9_CCLIST})
   target_compile_options (tinker9_cu PRIVATE
      "SHELL:-gencode arch=compute_${var},code=sm_${var}"
   )
endforeach () # -gencode arch=compute_60,code=sm_60 -gencode arch=compute_70,code=sm_70
## Debug add flag: -lineinfo
## Release add flag: --use_fast_math
target_compile_options (tinker9_cu PRIVATE
   $<$<CONFIG:DEBUG>:-lineinfo>
   $<$<CONFIG:RELEASE>:--use_fast_math>
)
