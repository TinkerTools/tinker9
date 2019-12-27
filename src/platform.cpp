#include "platform.h"
#include "tinker_rt.h"


TINKER_NAMESPACE_BEGIN
void platform_data(rc_op op)
{
   if (op & rc_dealloc) {
      platform::config = 0;
   }


   if (op & rc_init) {
      platform::config = platform::acc_pltfm;
      std::string gpu_package;
      get_kv_pair("GPU-PACKAGE", gpu_package, "OPENACC");
      if (gpu_package == "CUDA")
         platform::config |= platform::cu_pltfm;
   }
}
TINKER_NAMESPACE_END
