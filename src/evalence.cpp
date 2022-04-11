#include "tool/externfunc.h"
#include "tool/platform.h"

namespace tinker {
TINKER_F2VOID(cu, 1, acc, 0, evalence, int);
void evalence(int vers)
{
   TINKER_F2CALL(cu, 1, acc, 0, evalence, vers);
}
}
