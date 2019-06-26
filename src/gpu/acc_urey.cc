#include "gpu/acc.h"
#include "gpu/decl_mdstate.h"
#include "gpu/e_urey.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
template <int USE>
void eurey_tmpl() {
  constexpr int do_e = USE & use_energy;
  constexpr int do_g = USE & use_grad;
  constexpr int do_v = USE & use_virial;
  sanity_check<USE>();

  #pragma acc serial deviceptr(eub,vir_eub)
  {
    if_constexpr(do_e) { *eub = 0; }
    if_constexpr(do_v) {
      for (int i = 0; i < 9; ++i)
        vir_eub[i] = 0;
    }
  }

  #pragma acc parallel loop independent\
              deviceptr(x,y,z,gx,gy,gz,\
              iury,uk,ul,\
              eub,vir_eub)
  for (int i = 0; i < nurey; ++i) {
    const int ia = iury[i][0];
    const int ic = iury[i][2];
    const real ideal = ul[i];
    const real force = uk[i];

    real xac = x[ia] - x[ic];
    real yac = y[ia] - y[ic];
    real zac = z[ia] - z[ic];

    real rac = REAL_SQRT(xac * xac + yac * yac + zac * zac);
    real dt = rac - ideal;
    real dt2 = dt * dt;

    if_constexpr(do_e) {
      real e = ureyunit * force * dt2 * (1 + cury * dt + qury * dt2);
      #pragma acc atomic update
      *eub += e;
    }

    if_constexpr(do_g) {
      real deddt =
          2 * ureyunit * force * dt * (1 + 1.5f * cury * dt + 2 * qury * dt2);
      real de = deddt * REAL_RECIP(rac);
      real dedx = de * xac;
      real dedy = de * yac;
      real dedz = de * zac;

      #pragma acc atomic update
      gx[ia] += dedx;
      #pragma acc atomic update
      gy[ia] += dedy;
      #pragma acc atomic update
      gz[ia] += dedz;
      #pragma acc atomic update
      gx[ic] -= dedx;
      #pragma acc atomic update
      gy[ic] -= dedy;
      #pragma acc atomic update
      gz[ic] -= dedz;

      if_constexpr(do_v) {
        real vxx = xac * dedx;
        real vyx = yac * dedx;
        real vzx = zac * dedx;
        real vyy = yac * dedy;
        real vzy = zac * dedy;
        real vzz = zac * dedz;

        #pragma acc atomic update
        vir_eub[_xx] += vxx;
        #pragma acc atomic update
        vir_eub[_yx] += vyx;
        #pragma acc atomic update
        vir_eub[_zx] += vzx;
        #pragma acc atomic update
        vir_eub[_xy] += vyx;
        #pragma acc atomic update
        vir_eub[_yy] += vyy;
        #pragma acc atomic update
        vir_eub[_zy] += vzy;
        #pragma acc atomic update
        vir_eub[_xz] += vzx;
        #pragma acc atomic update
        vir_eub[_yz] += vzy;
        #pragma acc atomic update
        vir_eub[_zz] += vzz;
      }
    }
  } // end for (int i)
}

void eurey_acc_impl__(int vers) {
  if (vers == v0 || vers == v3)
    eurey_tmpl<v0>();
  else if (vers == v1)
    eurey_tmpl<v1>();
  else if (vers == v4)
    eurey_tmpl<v4>();
  else if (vers == v5)
    eurey_tmpl<v5>();
  else if (vers == v6)
    eurey_tmpl<v6>();
}
}
TINKER_NAMESPACE_END
