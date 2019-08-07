#include "acc_seq.h"
#include "e_urey.h"
#include "md.h"

TINKER_NAMESPACE_BEGIN
template <int USE>
void eurey_tmpl() {
  constexpr int do_e = USE & calc::energy;
  constexpr int do_g = USE & calc::grad;
  constexpr int do_v = USE & calc::virial;
  sanity_check<USE>();

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
        vir_eub[0] += vxx;
        #pragma acc atomic update
        vir_eub[1] += vyx;
        #pragma acc atomic update
        vir_eub[2] += vzx;
        #pragma acc atomic update
        vir_eub[3] += vyx;
        #pragma acc atomic update
        vir_eub[4] += vyy;
        #pragma acc atomic update
        vir_eub[5] += vzy;
        #pragma acc atomic update
        vir_eub[6] += vzx;
        #pragma acc atomic update
        vir_eub[7] += vzy;
        #pragma acc atomic update
        vir_eub[8] += vzz;
      }
    }
  } // end for (int i)
}

void eurey_acc_impl_(int vers) {
  if (vers == calc::v0 || vers == calc::v3)
    eurey_tmpl<calc::v0>();
  else if (vers == calc::v1)
    eurey_tmpl<calc::v1>();
  else if (vers == calc::v4)
    eurey_tmpl<calc::v4>();
  else if (vers == calc::v5)
    eurey_tmpl<calc::v5>();
  else if (vers == calc::v6)
    eurey_tmpl<calc::v6>();
}
TINKER_NAMESPACE_END
