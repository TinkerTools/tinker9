#include "tool/externfunc.h"
#include "tool/platform.h"

namespace tinker {
TINKER_F2EXTN(cu, 1, acc, 0, void, evalence, int);
void evalence(int vers)
{
#if TINKER_CUDART
   if (pltfm_config & Platform::CUDA)
      TINKER_F1CALL(cu, evalence, vers);
#endif
}
}
