#include "evalence.h"
#include "platform.h"
#include "tool/error.h"

namespace tinker {
void evalence_cu(int vers);
void evalence(int vers)
{
#if TINKER_CUDART
   if (pltfm_config & CU_PLTFM)
      evalence_cu(vers);
   else
#endif
      TINKER_THROW("Combined valence energy term should not have been called.\n");
}
}
