#include "platform.h"
#include "tinker_rt.h"


TINKER_NAMESPACE_BEGIN
void platform_data(rc_op op)
{
   if (op & rc_dealloc) {
      platform::config = platform::NOT_SET;
   }


   if (op & rc_init) {
#if TINKER_HOST
      platform::config = platform::ACC_PLTFM;
#endif
#if TINKER_CUDART
      // Feature: If the platform has been hard-coded, do not change it.
      if (platform::config == platform::NOT_SET) {
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
            get_kv_pair("GPU-PACKAGE", gpu_package, "CUDA");
         }
         if (gpu_package == "CUDA") {
            platform::config = platform::CU_PLTFM;
            print(stdout, " Use the CUDA Platform\n");
         } else if (gpu_package == "OPENACC") {
            platform::config = platform::ACC_PLTFM;
            print(stdout, " Use the OpenACC Platform\n");
         }
      }
#endif
   }
}
TINKER_NAMESPACE_END
