#include "acc_add.h"
#include "e_angle.h"
#include "md.h"

// TODO: test linear and fourier angles
// TODO: merge the branching

/**
 * Comment
 * Zhi Wang, July 1, 2019
 *
 * The original Tinker implementation has following expressions
 * xip = xib + xt * delta
 * xap = xia - xip
 *
 * And they were reorganized to
 * xap = xia - xib - xt * delta
 * for higher accuracy in the single precision mode.
 *
 * Consider an example where
 * xia = 33.553368, xib = 34.768604
 * xt = 0.33142909, delta = 0.0044494048,
 * the later expression gives a better numerical result.
 */

TINKER_NAMESPACE_BEGIN
template <int USE>
void eangle_tmpl() {
  constexpr int do_e = USE & calc::energy;
  constexpr int do_g = USE & calc::grad;
  constexpr int do_v = USE & calc::virial;
  sanity_check<USE>();

  auto* ea = ea_handle.e()->buffer();
  auto* vir_ea = ea_handle.vir()->buffer();
  auto bufsize = ea_handle.buffer_size();

  #pragma acc parallel loop gang(bufsize) independent\
              deviceptr(x,y,z,gx,gy,gz,\
              iang,anat,ak,angtyp,\
              ea,vir_ea)
  for (int i = 0; i < nangle; ++i) {
    int offset = i & (bufsize - 1);
    int ia = iang[i][0];
    int ib = iang[i][1];
    int ic = iang[i][2];
    int id = iang[i][3];
    real ideal = anat[i];
    real force = ak[i];
    eangle_t angtypi = angtyp[i];

    real xia = x[ia];
    real yia = y[ia];
    real zia = z[ia];
    real xib = x[ib];
    real yib = y[ib];
    real zib = z[ib];
    real xic = x[ic];
    real yic = y[ic];
    real zic = z[ic];

    if (angtypi != eangle_t::in_plane) {
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
        if (angtypi == eangle_t::harmonic) {
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

        if_constexpr(do_e) { atomic_add_value(e, ea, offset); }

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

            int offv = offset * 16;
            atomic_add_value(vxx, vir_ea, offv + 0);
            atomic_add_value(vyx, vir_ea, offv + 1);
            atomic_add_value(vzx, vir_ea, offv + 2);
            atomic_add_value(vyx, vir_ea, offv + 3);
            atomic_add_value(vyy, vir_ea, offv + 4);
            atomic_add_value(vzy, vir_ea, offv + 5);
            atomic_add_value(vzx, vir_ea, offv + 6);
            atomic_add_value(vzy, vir_ea, offv + 7);
            atomic_add_value(vzz, vir_ea, offv + 8);
          }
        }
      }
    } else {
      real xid = x[id];
      real yid = y[id];
      real zid = z[id];
      real xad = xia - xid;
      real yad = yia - yid;
      real zad = zia - zid;
      real xbd = xib - xid;
      real ybd = yib - yid;
      real zbd = zib - zid;
      real xcd = xic - xid;
      real ycd = yic - yid;
      real zcd = zic - zid;
      real xt = yad * zcd - zad * ycd;
      real yt = zad * xcd - xad * zcd;
      real zt = xad * ycd - yad * xcd;
      real rt2 = xt * xt + yt * yt + zt * zt;
      real delta = -(xt * xbd + yt * ybd + zt * zbd) * REAL_RECIP(rt2);
      // real xip = xib + xt * delta;
      // real yip = yib + yt * delta;
      // real zip = zib + zt * delta;
      real xap = xia - xib - xt * delta;
      real yap = yia - yib - yt * delta;
      real zap = zia - zib - zt * delta;
      real xcp = xic - xib - xt * delta;
      real ycp = yic - yib - yt * delta;
      real zcp = zic - zib - zt * delta;
      real rap2 = xap * xap + yap * yap + zap * zap;
      real rcp2 = xcp * xcp + ycp * ycp + zcp * zcp;
      if (rap2 != 0 && rcp2 != 0) {
        real dot = xap * xcp + yap * ycp + zap * zcp;
        real cosine = dot * REAL_RSQRT(rap2 * rcp2);
        cosine = REAL_MIN(1, REAL_MAX(-1, cosine));
        real angle = radian * REAL_ACOS(cosine);
        real dt = angle - ideal;
        real dt2 = dt * dt;
        real dt3 = dt2 * dt;
        real dt4 = dt2 * dt2;

        if_constexpr(do_e) {
          real e = angunit * force * dt2 *
              (1 + cang * dt + qang * dt2 + pang * dt3 + sang * dt4);
          atomic_add_value(e, ea, offset);
        }

        if_constexpr(do_g) {
          real deddt = angunit * force * dt * radian *
              (2 + 3 * cang * dt + 4 * qang * dt2 + 5 * pang * dt3 +
               6 * sang * dt4);
          real xm = ycp * zap - zcp * yap;
          real ym = zcp * xap - xcp * zap;
          real zm = xcp * yap - ycp * xap;
          real rm = REAL_SQRT(xm * xm + ym * ym + zm * zm);
          rm = REAL_MAX(rm, (real)0.0001);

          // chain rule terms for first derivative components

          real terma = -deddt * REAL_RECIP(rap2 * rm);
          real termc = deddt * REAL_RECIP(rcp2 * rm);
          real dedxia = terma * (yap * zm - zap * ym);
          real dedyia = terma * (zap * xm - xap * zm);
          real dedzia = terma * (xap * ym - yap * xm);
          real dedxic = termc * (ycp * zm - zcp * ym);
          real dedyic = termc * (zcp * xm - xcp * zm);
          real dedzic = termc * (xcp * ym - ycp * xm);
          real dedxip = -dedxia - dedxic;
          real dedyip = -dedyia - dedyic;
          real dedzip = -dedzia - dedzic;

          // chain rule components for the projection of the central atom

          real delta2, term;
          delta2 = 2 * delta;
          real ptrt2 =
              (dedxip * xt + dedyip * yt + dedzip * zt) * REAL_RECIP(rt2);
          term = (zcd * ybd - ycd * zbd) + delta2 * (yt * zcd - zt * ycd);
          real dpdxia = delta * (ycd * dedzip - zcd * dedyip) + term * ptrt2;
          term = (xcd * zbd - zcd * xbd) + delta2 * (zt * xcd - xt * zcd);
          real dpdyia = delta * (zcd * dedxip - xcd * dedzip) + term * ptrt2;
          term = (ycd * xbd - xcd * ybd) + delta2 * (xt * ycd - yt * xcd);
          real dpdzia = delta * (xcd * dedyip - ycd * dedxip) + term * ptrt2;
          term = (yad * zbd - zad * ybd) + delta2 * (zt * yad - yt * zad);
          real dpdxic = delta * (zad * dedyip - yad * dedzip) + term * ptrt2;
          term = (zad * xbd - xad * zbd) + delta2 * (xt * zad - zt * xad);
          real dpdyic = delta * (xad * dedzip - zad * dedxip) + term * ptrt2;
          term = (xad * ybd - yad * xbd) + delta2 * (yt * xad - xt * yad);
          real dpdzic = delta * (yad * dedxip - xad * dedyip) + term * ptrt2;

          // compute derivative components for this interaction

          dedxia = dedxia + dpdxia;
          dedyia = dedyia + dpdyia;
          dedzia = dedzia + dpdzia;
          real dedxib = dedxip;
          real dedyib = dedyip;
          real dedzib = dedzip;
          dedxic = dedxic + dpdxic;
          dedyic = dedyic + dpdyic;
          dedzic = dedzic + dpdzic;
          real dedxid = -dedxia - dedxib - dedxic;
          real dedyid = -dedyia - dedyib - dedyic;
          real dedzid = -dedzia - dedzib - dedzic;

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
          #pragma acc atomic update
          gx[id] += dedxid;
          #pragma acc atomic update
          gy[id] += dedyid;
          #pragma acc atomic update
          gz[id] += dedzid;

          if_constexpr(do_v) {
            real vxx = xad * dedxia + xbd * dedxib + xcd * dedxic;
            real vyx = yad * dedxia + ybd * dedxib + ycd * dedxic;
            real vzx = zad * dedxia + zbd * dedxib + zcd * dedxic;
            real vyy = yad * dedyia + ybd * dedyib + ycd * dedyic;
            real vzy = zad * dedyia + zbd * dedyib + zcd * dedyic;
            real vzz = zad * dedzia + zbd * dedzib + zcd * dedzic;

            int offv = offset * 16;
            atomic_add_value(vxx, vir_ea, offv + 0);
            atomic_add_value(vyx, vir_ea, offv + 1);
            atomic_add_value(vzx, vir_ea, offv + 2);
            atomic_add_value(vyx, vir_ea, offv + 3);
            atomic_add_value(vyy, vir_ea, offv + 4);
            atomic_add_value(vzy, vir_ea, offv + 5);
            atomic_add_value(vzx, vir_ea, offv + 6);
            atomic_add_value(vzy, vir_ea, offv + 7);
            atomic_add_value(vzz, vir_ea, offv + 8);
          }
        }
      }
    }
  } // end for (int i)
}

void eangle_acc_impl_(int vers) {
  if (vers == calc::v0 || vers == calc::v3)
    eangle_tmpl<calc::v0>();
  else if (vers == calc::v1)
    eangle_tmpl<calc::v1>();
  else if (vers == calc::v4)
    eangle_tmpl<calc::v4>();
  else if (vers == calc::v5)
    eangle_tmpl<calc::v5>();
  else if (vers == calc::v6)
    eangle_tmpl<calc::v6>();
}
TINKER_NAMESPACE_END
