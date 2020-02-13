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
         platform::config = platform::CU_PLTFM;
         std::string gpu_package;
         get_kv_pair("GPU-PACKAGE", gpu_package, "CUDA");
         if (gpu_package == "OPENACC")
            platform::config = platform::ACC_PLTFM;
      }
#endif
   }
}
TINKER_NAMESPACE_END
