#include "platform.h"
#include "tinker9.h"
#include "tool/ioprint.h"
#include "tool/iotext.h"

namespace tinker {
void platformData(RcOp op)
{
   if (op & RcOp::DEALLOC) {
      pltfm_config = Platform::UNSET;
   }

   if (op & RcOp::INIT) {
#if TINKER_HOST
      pltfm_config = Platform::ACC;
#endif
#if TINKER_CUDART
      // Feature: If the platform has been hard-coded, do not change it.
      if (pltfm_config == Platform::UNSET) {
         std::string gpu_package = "";
         if (const char* str = std::getenv("gpu_package")) {
            gpu_package = str;
            Text::upcase(gpu_package);
         }
         if (const char* str = std::getenv("GPU_PACKAGE")) {
            gpu_package = str;
            Text::upcase(gpu_package);
         }
         if (gpu_package == "") {
            get_kv("GPU-PACKAGE", gpu_package, "CUDA");
         }
         if (gpu_package == "CUDA") {
            pltfm_config = Platform::CUDA;
            print(stdout, " Primary GPU package :  CUDA\n");
         } else if (gpu_package == "OPENACC") {
            pltfm_config = Platform::ACC;
            print(stdout, " Primary GPU package :  OpenACC\n");
         }
      }
#endif
   }
}
}
