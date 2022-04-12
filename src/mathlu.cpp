#include "tool/externfunc.h"

namespace tinker {
template <int n, class R>
void symlusolve_acc(const R* aUpRowMajor, R* b);
}

namespace tinker {
// TINKER_F2VOID(cu, 0, acc, 1, symlusolve);
template <int n, class R>
void symlusolve(const R* aUpRowMajor, R* b)
{
// TINKER_F2CALL(cu, 0, acc, 1, symlusolve);
#if TINKER_GPULANG_OPENACC
   symlusolve_acc<n, R>(aUpRowMajor, b);
#elif TINKER_GPULANG_CUDA
   throwExceptionMissingFunction("symlusolve_cu<...>");
#else
   symlusolve_acc<n, R>(aUpRowMajor, b);
#endif
}

template void symlusolve<3, float>(const float*, float*);
template void symlusolve<6, float>(const float*, float*);
template void symlusolve<3, double>(const double*, double*);
template void symlusolve<6, double>(const double*, double*);
}
