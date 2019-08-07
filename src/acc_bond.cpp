#include "acc_seq.h"
#include "gpu/e_bond.h"
#include "md.h"
#include <cassert>

// TODO: test morse potential

TINKER_NAMESPACE_BEGIN
template <int USE, int BNDTYP>
void ebond_tmpl() {
  constexpr int do_e = USE & calc::energy;
  constexpr int do_g = USE & calc::grad;
  constexpr int do_v = USE & calc::virial;
  sanity_check<USE>();

  #pragma acc parallel loop independent\
              deviceptr(x,y,z,gx,gy,gz,\
              ibnd,bl,bk,\
              eb,vir_eb)
  for (int i = 0; i < nbond; ++i) {
    int ia = ibnd[i][0];
    int ib = ibnd[i][1];
    real ideal = bl[i];
    real force = bk[i];

    real xab = x[ia] - x[ib];
    real yab = y[ia] - y[ib];
    real zab = z[ia] - z[ib];

    real rab = REAL_SQRT(xab * xab + yab * yab + zab * zab);
    real dt = rab - ideal;

    MAYBE_UNUSED real e;
    MAYBE_UNUSED real deddt;
    if_constexpr(BNDTYP == bond_harmonic) {
      real dt2 = dt * dt;
      if_constexpr(do_e) e =
          bndunit * force * dt2 * (1 + cbnd * dt + qbnd * dt2);
      if_constexpr(do_g) deddt =
          2 * bndunit * force * dt * (1 + 1.5f * cbnd * dt + 2 * qbnd * dt2);
    }
    else if_constexpr(BNDTYP == bond_morse) {
      real expterm = REAL_EXP(-2 * dt);
      real bde = 0.25f * bndunit * force;
      if_constexpr(do_e) e = bde * (1 - expterm) * (1 - expterm);
      if_constexpr(do_g) deddt = 4 * bde * (1 - expterm) * expterm;
    }

    if_constexpr(do_e) {
      #pragma acc atomic update
      *eb += e;
    }

    if_constexpr(do_g) {
      real de = deddt * REAL_RECIP(rab);
      real dedx = de * xab;
      real dedy = de * yab;
      real dedz = de * zab;
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

      if_constexpr(do_v) {
        real vxx = xab * dedx;
        real vyx = yab * dedx;
        real vzx = zab * dedx;
        real vyy = yab * dedy;
        real vzy = zab * dedy;
        real vzz = zab * dedz;

        #pragma acc atomic update
        vir_eb[0] += vxx;
        #pragma acc atomic update
        vir_eb[1] += vyx;
        #pragma acc atomic update
        vir_eb[2] += vzx;
        #pragma acc atomic update
        vir_eb[3] += vyx;
        #pragma acc atomic update
        vir_eb[4] += vyy;
        #pragma acc atomic update
        vir_eb[5] += vzy;
        #pragma acc atomic update
        vir_eb[6] += vzx;
        #pragma acc atomic update
        vir_eb[7] += vzy;
        #pragma acc atomic update
        vir_eb[8] += vzz;
      }
    }
  } // end for (int i)
}

void ebond_acc_impl_(int vers) {
  if (bndtyp == bond_harmonic)
    if (vers == calc::v0 || vers == calc::v3)
      ebond_tmpl<calc::v0, bond_harmonic>();
    else if (vers == calc::v1)
      ebond_tmpl<calc::v1, bond_harmonic>();
    else if (vers == calc::v4)
      ebond_tmpl<calc::v4, bond_harmonic>();
    else if (vers == calc::v5)
      ebond_tmpl<calc::v5, bond_harmonic>();
    else if (vers == calc::v6)
      ebond_tmpl<calc::v6, bond_harmonic>();
    else
      assert(false);
  else if (bndtyp == bond_morse)
    if (vers == calc::v0 || vers == calc::v3)
      ebond_tmpl<calc::v0, bond_morse>();
    else if (vers == calc::v1)
      ebond_tmpl<calc::v1, bond_morse>();
    else if (vers == calc::v4)
      ebond_tmpl<calc::v4, bond_morse>();
    else if (vers == calc::v5)
      ebond_tmpl<calc::v5, bond_morse>();
    else if (vers == calc::v6)
      ebond_tmpl<calc::v6, bond_morse>();
    else
      assert(false);
  else
    assert(false);
}
TINKER_NAMESPACE_END
