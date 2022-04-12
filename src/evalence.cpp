#include "tool/externfunc.h"
#include "tool/platform.h"

namespace tinker {
TINKER_FVOID2(cu, 1, acc, 0, evalence, int);
void evalence(int vers)
{
   TINKER_FCALL2(cu, 1, acc, 0, evalence, vers);
}
}
