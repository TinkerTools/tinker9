#include "gpu/decl_mdstate.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
void sum_energy_acc_impl__(real* ebuf, int top) {
  #pragma acc serial deviceptr(ebuf,esum)
  {
    for (int i = 1; i < top; ++i)
      *esum += ebuf[i];
  }
}

void sum_virial_acc_impl__(real* vbuf, int top, int virlen) {
  #pragma acc serial deviceptr(vbuf,vir)
  {
    for (int i = 1; i < top; ++i)
      for (int j = 0; j < 9; ++j)
        vir[j] += vbuf[i * virlen + j];
  }
}
}
TINKER_NAMESPACE_END
