#include "gpu/e.bond.h"
#include "gpu/mdstate.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
template <int USE, int BNDTYP>
void ebond_tmpl() {
  constexpr int do_e = USE & use_energy;
  constexpr int do_g = USE & use_grad;
  constexpr int do_v = USE & use_virial;

  real eb = 0;

// clang-format off
#pragma acc data deviceptr(x,y,z,gx,gy,gz,\
                           ibnd,bl,bk)
#pragma acc data copy(eb,vir_xx,vir_yx,vir_zx,\
                         vir_xy,vir_yy,vir_zy,\
                         vir_xz,vir_yz,vir_zz)
#pragma acc parallel loop reduction(+:eb,vir_xx,vir_yx,vir_zx,\
                                         vir_xy,vir_yy,vir_zy,\
                                         vir_xz,vir_yz,vir_zz)
  // clang-format on
  for (int i = 0; i < nbond; ++i) {
    int ia = ibnd[i][0];
    int ib = ibnd[i][1];
    real ideal = bl[i];
    real force = bk[i];

    real xab = x[ia] - x[ib];
    real yab = y[ia] - y[ib];
    real zab = z[ia] - z[ib];

    // FIXME image
    real rab = REAL_SQRT(xab * xab + yab * yab + zab * zab);
    real dt = rab - ideal;

    real e;
    real deddt;
    if_constexpr(BNDTYP & ebond_harmonic) {
      real dt2 = dt * dt;
      if_constexpr(do_e) e =
          bndunit * force * dt2 * (1 + cbnd * dt + qbnd * dt2);
      if_constexpr(do_g) deddt =
          2 * bndunit * force * dt * (1 + 1.5f * cbnd * dt + 2 * qbnd * dt2);
    }
    else if_constexpr(BNDTYP & ebond_morse) {
      real expterm = REAL_EXP(-2 * dt);
      real bde = 0.25f * bndunit * force;
      if_constexpr(do_e) e = bde * (1 - expterm) * (1 - expterm);
      if_constexpr(do_g) deddt = 4 * bde * (1 - expterm) * expterm;
    }

    if_constexpr(do_e) eb += e;

    real de;
    real dedx, dedy, dedz;
    if_constexpr(do_g) {
      de = deddt / rab;
      dedx = de * xab;
      dedy = de * yab;
      dedz = de * zab;
#pragma acc atomic update
      gx[ia] += dedx;
#pragma acc atomic update
      gy[ia] += dedy;
#pragma acc atomic update
      gz[ia] += dedz;
#pragma acc atomic update
      gx[ib] -= dedx;
#pragma acc atomic update
      gy[ib] -= dedy;
#pragma acc atomic update
      gz[ib] -= dedz;
    }

    if_constexpr(do_v) {
      real vxx = xab * dedx;
      real vyx = yab * dedx;
      real vzx = zab * dedx;
      real vyy = yab * dedy;
      real vzy = zab * dedy;
      real vzz = zab * dedz;

      vir_xx += vxx;
      vir_yx += vyx;
      vir_zx += vzx;
      vir_xy += vyx;
      vir_yy += vyy;
      vir_zy += vzy;
      vir_xz += vzx;
      vir_yz += vzy;
      vir_zz += vzz;
    }
  }
}
}
TINKER_NAMESPACE_END

m_tinker_using_namespace;
extern "C" {
TINKER_BONDED_GEN(tinker_gpu_ebond_harmonic, gpu::ebond_tmpl,
                  gpu::ebond_harmonic);
TINKER_BONDED_GEN(tinker_gpu_ebond_morse, gpu::ebond_tmpl, gpu::ebond_morse);
}
