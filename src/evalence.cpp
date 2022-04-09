#include "tool/externfunc.h"
#include "tool/platform.h"

namespace tinker {
TINKER_F2EXTN(void, evalence, cu, 1, acc, 0, int);
void evalence(int vers)
{
#if TINKER_CUDART
   TINKER_F1CALL(evalence, cu, vers);
#else
   TINKER_F1CALL(evalence, acc, vers);
#endif
}
}
