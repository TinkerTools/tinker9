#include "platform.h"
#include "tinker_rt.h"

namespace tinker {
int pltfm_config;

void platform_data(RcOp op)
{
   if (op & rc_dealloc) {
      pltfm_config = UNSET_PLTFM;
   }

   if (op & rc_init) {
#if TINKER_HOST
      pltfm_config = ACC_PLTFM;
#endif
#if TINKER_CUDART
      // Feature: If the platform has been hard-coded, do not change it.
      if (pltfm_config == UNSET_PLTFM) {
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
            pltfm_config = CU_PLTFM;
            print(stdout, " Primary GPU package :  CUDA\n");
         } else if (gpu_package == "OPENACC") {
            pltfm_config = ACC_PLTFM;
            print(stdout, " Primary GPU package :  OpenACC\n");
         }
      }
#endif
   }
}
}
