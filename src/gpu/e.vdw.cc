#include "gpu/e.vdw.h"
#include "gpu/acc.h"
#include "gpu/image.h"
#include "gpu/mdstate.h"
#include "gpu/nblist.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
template <int USE, int VDWTYP>
void evdw_tmpl() {
  constexpr int do_e = USE & use_energy;
  constexpr int do_g = USE & use_grad;
  constexpr int do_v = USE & use_virial;
  constexpr int do_a = USE & use_analyz;

  #pragma acc data deviceptr(xred,yred,zred,gx,gy,gz,vir,box,\
                             ev)
  {
    #pragma acc serial async(queue_nb)
    { *ev = 0; }
  }
}
}
TINKER_NAMESPACE_END

extern "C" {
m_tinker_using_namespace;
TINKER_NONBONDED_GEN(tinker_gpu_evdw_lj, gpu::evdw_tmpl, gpu::evdw_lj);
TINKER_NONBONDED_GEN(tinker_gpu_evdw_buck, gpu::evdw_tmpl, gpu::evdw_buck);
TINKER_NONBONDED_GEN(tinker_gpu_evdw_mm3hb, gpu::evdw_tmpl, gpu::evdw_mm3hb);
TINKER_NONBONDED_GEN(tinker_gpu_evdw_hal, gpu::evdw_tmpl, gpu::evdw_hal);
TINKER_NONBONDED_GEN(tinker_gpu_evdw_gauss, gpu::evdw_tmpl, gpu::evdw_gauss);

#define TINKER_GPU_EVDW_GEN_(ver)                                              \
  void tinker_gpu_evdw##ver() {                                                \
    if (gpu::vdwtyp == gpu::evdw_buck)                                         \
      tinker_gpu_evdw_buck##ver();                                             \
    else if (gpu::vdwtyp == gpu::evdw_mm3hb)                                   \
      tinker_gpu_evdw_mm3hb##ver();                                            \
    else if (gpu::vdwtyp == gpu::evdw_hal)                                     \
      tinker_gpu_evdw_hal##ver();                                              \
    else if (gpu::vdwtyp == gpu::evdw_gauss)                                   \
      tinker_gpu_evdw_gauss##ver();                                            \
    else                                                                       \
      tinker_gpu_evdw_lj##ver();                                               \
  }
TINKER_GPU_EVDW_GEN_(0);
TINKER_GPU_EVDW_GEN_(1);
TINKER_GPU_EVDW_GEN_(3);
TINKER_GPU_EVDW_GEN_(4);
TINKER_GPU_EVDW_GEN_(5);
TINKER_GPU_EVDW_GEN_(6);
#undef TINKER_GPU_EVDW_GEN_
}
