#include "gpu/decl_mdstate.h"
#include "gpu/e_polar.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
template <int USE>
void epolar_ewald_tmpl(const real (*gpu_uind)[3], const real (*gpu_uinp)[3]) {
  constexpr int do_e = USE & use_energy;
  constexpr int do_a = USE & use_analyz;
  constexpr int do_g = USE & use_grad;
  constexpr int do_v = USE & use_virial;
  static_assert(do_v ? do_g : true, "");
  static_assert(do_a ? do_e : true, "");

  if_constexpr(do_e || do_a || do_v) {
    #pragma acc serial deviceptr(ep,nep,vir_ep)
    {
      if_constexpr(do_e) { *ep = 0; }
      if_constexpr(do_a) { *nep = 0; }
      if_constexpr(do_v) {
        for (int i = 0; i < 9; ++i) {
          vir_ep[i] = 0;
        }
      }
    }
  }

  if_constexpr(USE == use_energy) { epolar0_dotprod(gpu_uind, udirp); }
}
}
TINKER_NAMESPACE_END

extern "C" {
m_tinker_using_namespace;
void tinker_gpu_epolar_ewald0() {
  gpu::induce(&gpu::uind[0][0], &gpu::uinp[0][0]);
  gpu::epolar_ewald_tmpl<gpu::v0>(gpu::uind, gpu::uinp);
}
void tinker_gpu_epolar_ewald1() {
  gpu::induce(&gpu::uind[0][0], &gpu::uinp[0][0]);
  gpu::epolar_ewald_tmpl<gpu::v1>(gpu::uind, gpu::uinp);
}
void tinker_gpu_epolar_ewald3() {
  gpu::induce(&gpu::uind[0][0], &gpu::uinp[0][0]);
  gpu::epolar_ewald_tmpl<gpu::v3>(gpu::uind, gpu::uinp);
}
void tinker_gpu_epolar_ewald4() {
  gpu::induce(&gpu::uind[0][0], &gpu::uinp[0][0]);
  gpu::epolar_ewald_tmpl<gpu::v4>(gpu::uind, gpu::uinp);
}
void tinker_gpu_epolar_ewald5() {
  gpu::induce(&gpu::uind[0][0], &gpu::uinp[0][0]);
  gpu::epolar_ewald_tmpl<gpu::v5>(gpu::uind, gpu::uinp);
}
void tinker_gpu_epolar_ewald6() {
  gpu::induce(&gpu::uind[0][0], &gpu::uinp[0][0]);
  gpu::epolar_ewald_tmpl<gpu::v6>(gpu::uind, gpu::uinp);
}
}
