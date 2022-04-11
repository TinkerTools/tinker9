#include "tool/externfunc.h"
#include "tool/rcman.h"

namespace tinker {
TINKER_F2VOID(cu, 0, acc, 1, cudalibData, RcOp);
void cudalibData(RcOp op)
{
   TINKER_F2CALL(cu, 0, acc, 1, cudalibData, op);
}
}
