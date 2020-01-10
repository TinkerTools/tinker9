#include "platform.h"
#include "tinker_rt.h"


TINKER_NAMESPACE_BEGIN
void platform_data(rc_op op)
{
   if (op & rc_dealloc) {
      platform::config = 0;
   }


   if (op & rc_init) {
#if TINKER_HOST
      platform::config = platform::acc_pltfm;
#endif
#if TINKER_CUDART
      platform::config = platform::cu_pltfm;
      std::string gpu_package;
      get_kv_pair("GPU-PACKAGE", gpu_package, "CUDA");
      if (gpu_package == "OPENACC")
         platform::config = platform::acc_pltfm;
#endif
   }
}
TINKER_NAMESPACE_END
