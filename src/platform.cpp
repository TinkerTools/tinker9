#include "tool/platform.h"
#include "tool/argkey.h"
#include "tool/ioprint.h"
#include "tool/iotext.h"

namespace tinker {
void platformData(RcOp op)
{
   if (op & RcOp::DEALLOC) {
      pltfm_config = Platform::UNKNOWN;
   }

   if (op & RcOp::INIT) {
#if TINKER_GPULANG_OPENACC
      // Feature: If the platform has been hard-coded, do not change it.
      if (pltfm_config == Platform::UNKNOWN) {
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
            getKV("GPU-PACKAGE", gpu_package, "CUDA");
         }
         if (gpu_package == "CUDA") {
            pltfm_config = Platform::CUDA;
            print(stdout, " Primary GPU package :  CUDA\n");
         } else if (gpu_package == "OPENACC") {
            pltfm_config = Platform::ACC;
            print(stdout, " Primary GPU package :  OpenACC\n");
         }
      }
#elif TINKER_GPULANG_CUDA
      pltfm_config = Platform::CUDA;
#else
      pltfm_config = Platform::ACC;
#endif
   }
}
}
