#include "gpu/decl_mathfunc.h"
#include "gpu/decl_mdstate.h"
#include "gpu/e_angle.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
template <int USE>
void eangle_tmpl() {
  constexpr int do_e = USE & use_energy;
  constexpr int do_g = USE & use_grad;
  constexpr int do_v = USE & use_virial;
  sanity_check<USE>();

  #pragma acc serial deviceptr(ea,vir_ea)
  {
    if_constexpr(do_e) { *ea = 0; }
    if_constexpr(do_v) {
      for (int i = 0; i < 9; ++i)
        vir_ea[i] = 0;
    }
  }

  #pragma acc parallel loop independent\
              deviceptr(x,y,z,gx,gy,gz,\
              iang,anat,ak,angtyp,\
              ea,vir_ea)
  for (int i = 0; i < nangle; ++i) {
    int ia = iang[i][0];
    int ib = iang[i][1];
    int ic = iang[i][2];
    // int id = iang[i][3];
    real ideal = anat[i];
    real force = ak[i];
    int angtypi = angtyp[i];

    real xia = x[ia];
    real yia = y[ia];
    real zia = z[ia];
    real xib = x[ib];
    real yib = y[ib];
    real zib = z[ib];
    real xic = x[ic];
    real yic = y[ic];
    real zic = z[ic];

    if (angtypi != angle_in_plane) {
      real xab = xia - xib;
      real yab = yia - yib;
      real zab = zia - zib;
      real xcb = xic - xib;
      real ycb = yic - yib;
      real zcb = zic - zib;

      real rab2 = xab * xab + yab * yab + zab * zab;
      real rcb2 = xcb * xcb + ycb * ycb + zcb * zcb;

      if (rab2 != 0 && rcb2 != 0) {
        real xp = ycb * zab - zcb * yab;
        real yp = zcb * xab - xcb * zab;
        real zp = xcb * yab - ycb * xab;
        real rp = REAL_SQRT(xp * xp + yp * yp + zp * zp);
        rp = REAL_MAX(rp, (real)0.0001);
        real dot = xab * xcb + yab * ycb + zab * zcb;
        real cosine = dot * REAL_RSQRT(rab2 * rcb2);
        cosine = REAL_MIN((real)0.1, REAL_MAX(-1, cosine));
        real angle = radian * REAL_ACOS(cosine);

        MAYBE_UNUSED real e;
        MAYBE_UNUSED real deddt;
        if (angtypi == angle_harmonic) {
          real dt = angle - ideal;
          real dt2 = dt * dt;
          real dt3 = dt2 * dt;
          real dt4 = dt2 * dt2;
          if_constexpr(do_e) e = angunit * force * dt2 *
              (1 + cang * dt + qang * dt2 + pang * dt3 + sang * dt4);
          if_constexpr(do_g) deddt = angunit * force * dt * radian *
              (2 + 3 * cang * dt + 4 * qang * dt2 + 5 * pang * dt3 +
               6 * sang * dt4);
        }

        if_constexpr(do_e) {
          #pragma acc atomic update
          *ea += e;
        }

        if_constexpr(do_g) {
          real terma = -deddt * REAL_RECIP(rab2 * rp);
          real termc = deddt * REAL_RECIP(rcb2 * rp);
          real dedxia = terma * (yab * zp - zab * yp);
          real dedyia = terma * (zab * xp - xab * zp);
          real dedzia = terma * (xab * yp - yab * xp);
          real dedxic = termc * (ycb * zp - zcb * yp);
          real dedyic = termc * (zcb * xp - xcb * zp);
          real dedzic = termc * (xcb * yp - ycb * xp);
          real dedxib = -dedxia - dedxic;
          real dedyib = -dedyia - dedyic;
          real dedzib = -dedzia - dedzic;

          #pragma acc atomic update
          gx[ia] += dedxia;
          #pragma acc atomic update
          gy[ia] += dedyia;
          #pragma acc atomic update
          gz[ia] += dedzia;
          #pragma acc atomic update
          gx[ib] += dedxib;
          #pragma acc atomic update
          gy[ib] += dedyib;
          #pragma acc atomic update
          gz[ib] += dedzib;
          #pragma acc atomic update
          gx[ic] += dedxic;
          #pragma acc atomic update
          gy[ic] += dedyic;
          #pragma acc atomic update
          gz[ic] += dedzic;

          if_constexpr(do_v) {
            real vxx = xab * dedxia + xcb * dedxic;
            real vyx = yab * dedxia + ycb * dedxic;
            real vzx = zab * dedxia + zcb * dedxic;
            real vyy = yab * dedyia + ycb * dedyic;
            real vzy = zab * dedyia + zcb * dedyic;
            real vzz = zab * dedzia + zcb * dedzic;

            #pragma acc atomic update
            vir_ea[_xx] += vxx;
            #pragma acc atomic update
            vir_ea[_yx] += vyx;
            #pragma acc atomic update
            vir_ea[_zx] += vzx;
            #pragma acc atomic update
            vir_ea[_xy] += vyx;
            #pragma acc atomic update
            vir_ea[_yy] += vyy;
            #pragma acc atomic update
            vir_ea[_zy] += vzy;
            #pragma acc atomic update
            vir_ea[_xz] += vzx;
            #pragma acc atomic update
            vir_ea[_yz] += vzy;
            #pragma acc atomic update
            vir_ea[_zz] += vzz;
          }
        }
      }
    }
  } // end for (int i)
}

void eangle_acc_impl__(int vers) {
  if (vers == v0 || vers == v3)
    eangle_tmpl<v0>();
  else if (vers == v1)
    eangle_tmpl<v1>();
  else if (vers == v4)
    eangle_tmpl<v4>();
  else if (vers == v5)
    eangle_tmpl<v5>();
  else if (vers == v6)
    eangle_tmpl<v6>();
}
}
TINKER_NAMESPACE_END
