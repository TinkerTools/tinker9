#include "acc_add.h"
#include "e_urey.h"
#include "md.h"

TINKER_NAMESPACE_BEGIN
template <int USE>
void eurey_tmpl() {
  constexpr int do_e = USE & calc::energy;
  constexpr int do_g = USE & calc::grad;
  constexpr int do_v = USE & calc::virial;
  sanity_check<USE>();

  auto* eub = eub_handle.e()->buffer();
  auto* vir_eub = eub_handle.vir()->buffer();
  auto bufsize = eub_handle.buffer_size();

  #pragma acc parallel num_gangs(bufsize)\
              deviceptr(x,y,z,gx,gy,gz,\
              iury,uk,ul,\
              eub,vir_eub)
  #pragma acc loop gang independent
  for (int i = 0; i < nurey; ++i) {
    int offset = i & (bufsize - 1);
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
      atomic_add_value(e, eub, offset);
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

        atomic_add_value(vxx, vyx, vzx, vyy, vzy, vzz, vir_eub, offset);
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
