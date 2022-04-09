#include "tool/externfunc.h"
#include "tool/platform.h"

namespace tinker {
TINKER_F2EXTN(cu, 1, acc, 0, void, evalence, int);
void evalence(int vers)
{
#if TINKER_CUDART
   TINKER_F1CALL(cu, evalence, vers);
#else
   TINKER_F1CALL(acc, evalence, vers);
#endif
}
}
