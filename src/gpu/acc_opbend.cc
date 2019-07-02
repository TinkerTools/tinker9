#include "gpu/acc.h"
#include "gpu/decl_mdstate.h"
#include "gpu/e_angle.h"
#include "gpu/e_opbend.h"

// TODO: test W-D-C

/**
 * Comments
 * Zhi Wang, Jun 25, 2019
 *
 * The original implementation in Tinker uses ACOS(cosine) to calculate the
 * out-of-plane angle, which is the major source of error in the single
 * precision mode.
 *
 * These angles (theta) are usually very small (e.g. 0.001 rad), so the value of
 * variable cosine is very close to 1 (cosine = SQRT(1 - eps**2)). As a result,
 * it is much more accurate to use theta = ASIN(eps) instead of theta =
 * ACOS(cosine) to calculate the angles.
 */

TINKER_NAMESPACE_BEGIN
namespace gpu {
template <int USE, eopbend_t TYP>
void eopbend_tmpl() {
  constexpr int do_e = USE & use_energy;
  constexpr int do_g = USE & use_grad;
  constexpr int do_v = USE & use_virial;
  sanity_check<USE>();

  #pragma acc serial deviceptr(eopb,vir_eopb)
  {
    if_constexpr(do_e) { *eopb = 0; }
    if_constexpr(do_v) {
      for (int i = 0; i < 9; ++i)
        vir_eopb[i] = 0;
    }
  }

  #pragma acc parallel loop independent\
              deviceptr(x,y,z,gx,gy,gz,\
              iopb,opbk,iang,\
              eopb,vir_eopb)
  for (int iopbend = 0; iopbend < nopbend; ++iopbend) {
    const real force = opbk[iopbend];
    const int i = iopb[iopbend];
    const int ia = iang[i][0];
    const int ib = iang[i][1];
    const int ic = iang[i][2];
    const int id = iang[i][3];

    real xia = x[ia];
    real yia = y[ia];
    real zia = z[ia];
    real xib = x[ib];
    real yib = y[ib];
    real zib = z[ib];
    real xic = x[ic];
    real yic = y[ic];
    real zic = z[ic];
    real xid = x[id];
    real yid = y[id];
    real zid = z[id];

    real xab = xia - xib;
    real yab = yia - yib;
    real zab = zia - zib;
    real xcb = xic - xib;
    real ycb = yic - yib;
    real zcb = zic - zib;
    real xdb = xid - xib;
    real ydb = yid - yib;
    real zdb = zid - zib;
    real xad = xia - xid;
    real yad = yia - yid;
    real zad = zia - zid;
    real xcd = xic - xid;
    real ycd = yic - yid;
    real zcd = zic - zid;

    MAYBE_UNUSED real rab2, rad2;
    MAYBE_UNUSED real rcb2, rcd2;
    MAYBE_UNUSED real dot;
    real cc;
    if_constexpr(TYP == opbend_w_d_c) {

      // W-D-C angle between A-B-C plane and B-D vector for D-B<AC

      rab2 = xab * xab + yab * yab + zab * zab;
      rcb2 = xcb * xcb + ycb * ycb + zcb * zcb;
      cc = rab2 * rcb2 - REAL_SQ(xab * xcb + yab * ycb + zab * zcb);
      if_constexpr(do_g) dot = xab * xcb + yab * ycb + zab * zcb;
    }
    else if_constexpr(TYP == opbend_allinger) {

      // Allinger angle between A-C-D plane and D-B vector for D-B<AC

      rad2 = xad * xad + yad * yad + zad * zad;
      rcd2 = xcd * xcd + ycd * ycd + zcd * zcd;
      cc = rad2 * rcd2 - REAL_SQ(xad * xcd + yad * ycd + zad * zcd);
      if_constexpr(do_g) dot = xad * xcd + yad * ycd + zad * zcd;
    }

    // find the out-of-plane angle bending energy

    real ee = xdb * (yab * zcb - zab * ycb) + ydb * (zab * xcb - xab * zcb) +
        zdb * (xab * ycb - yab * xcb);
    real rdb2 = xdb * xdb + ydb * ydb + zdb * zdb;

    if (rdb2 != 0 && cc != 0) {
      real sine = REAL_ABS(ee) * REAL_RSQRT(cc * rdb2);
      sine = REAL_MIN(1, sine);
      real angle = radian * REAL_ASIN(sine);
      real dt = angle;
      real dt2 = dt * dt;
      real dt3 = dt2 * dt;
      real dt4 = dt2 * dt2;

      if_constexpr(do_e) {
        real e = opbunit * force * dt2 *
            (1 + copb * dt + qopb * dt2 + popb * dt3 + sopb * dt4);
        #pragma acc atomic update
        *eopb += e;
      }

      if_constexpr(do_g) {
        real deddt = opbunit * force * dt * radian *
            (2 + 3 * copb * dt + 4 * qopb * dt2 + 5 * popb * dt3 +
             6 * sopb * dt4);
        real dedcos =
            -deddt * REAL_SIGN(1, ee) * REAL_RSQRT(cc * rdb2 - ee * ee);
        real term = ee * REAL_RECIP(cc);
        real dccdxia, dccdyia, dccdzia;
        real dccdxic, dccdyic, dccdzic;
        real dccdxid, dccdyid, dccdzid;
        if_constexpr(TYP == opbend_w_d_c) {
          dccdxia = (xab * rcb2 - xcb * dot) * term;
          dccdyia = (yab * rcb2 - ycb * dot) * term;
          dccdzia = (zab * rcb2 - zcb * dot) * term;
          dccdxic = (xcb * rab2 - xab * dot) * term;
          dccdyic = (ycb * rab2 - yab * dot) * term;
          dccdzic = (zcb * rab2 - zab * dot) * term;
          dccdxid = 0;
          dccdyid = 0;
          dccdzid = 0;
        }
        else if_constexpr(TYP == opbend_allinger) {
          dccdxia = (xad * rcd2 - xcd * dot) * term;
          dccdyia = (yad * rcd2 - ycd * dot) * term;
          dccdzia = (zad * rcd2 - zcd * dot) * term;
          dccdxic = (xcd * rad2 - xad * dot) * term;
          dccdyic = (ycd * rad2 - yad * dot) * term;
          dccdzic = (zcd * rad2 - zad * dot) * term;
          dccdxid = -dccdxia - dccdxic;
          dccdyid = -dccdyia - dccdyic;
          dccdzid = -dccdzia - dccdzic;
        }

        term = ee * REAL_RECIP(rdb2);
        real deedxia = ydb * zcb - zdb * ycb;
        real deedyia = zdb * xcb - xdb * zcb;
        real deedzia = xdb * ycb - ydb * xcb;
        real deedxic = yab * zdb - zab * ydb;
        real deedyic = zab * xdb - xab * zdb;
        real deedzic = xab * ydb - yab * xdb;
        real deedxid = ycb * zab - zcb * yab + xdb * term;
        real deedyid = zcb * xab - xcb * zab + ydb * term;
        real deedzid = xcb * yab - ycb * xab + zdb * term;

        // compute first derivative components for this angle

        real dedxia = dedcos * (dccdxia + deedxia);
        real dedyia = dedcos * (dccdyia + deedyia);
        real dedzia = dedcos * (dccdzia + deedzia);
        real dedxic = dedcos * (dccdxic + deedxic);
        real dedyic = dedcos * (dccdyic + deedyic);
        real dedzic = dedcos * (dccdzic + deedzic);
        real dedxid = dedcos * (dccdxid + deedxid);
        real dedyid = dedcos * (dccdyid + deedyid);
        real dedzid = dedcos * (dccdzid + deedzid);
        real dedxib = -dedxia - dedxic - dedxid;
        real dedyib = -dedyia - dedyic - dedyid;
        real dedzib = -dedzia - dedzic - dedzid;

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
          real vxx = xab * dedxia + xcb * dedxic + xdb * dedxid;
          real vyx = yab * dedxia + ycb * dedxic + ydb * dedxid;
          real vzx = zab * dedxia + zcb * dedxic + zdb * dedxid;
          real vyy = yab * dedyia + ycb * dedyic + ydb * dedyid;
          real vzy = zab * dedyia + zcb * dedyic + zdb * dedyid;
          real vzz = zab * dedzia + zcb * dedzic + zdb * dedzid;

          #pragma acc atomic update
          vir_eopb[_xx] += vxx;
          #pragma acc atomic update
          vir_eopb[_yx] += vyx;
          #pragma acc atomic update
          vir_eopb[_zx] += vzx;
          #pragma acc atomic update
          vir_eopb[_xy] += vyx;
          #pragma acc atomic update
          vir_eopb[_yy] += vyy;
          #pragma acc atomic update
          vir_eopb[_zy] += vzy;
          #pragma acc atomic update
          vir_eopb[_xz] += vzx;
          #pragma acc atomic update
          vir_eopb[_yz] += vzy;
          #pragma acc atomic update
          vir_eopb[_zz] += vzz;
        }
      }
    }
  } // end for (int iopbend)
}

void eopbend_acc_impl__(int vers) {
  if (opbtyp == opbend_w_d_c) {
    if (vers == v0 || vers == v3)
      eopbend_tmpl<v0, opbend_w_d_c>();
    else if (vers == v1)
      eopbend_tmpl<v1, opbend_w_d_c>();
    else if (vers == v4)
      eopbend_tmpl<v4, opbend_w_d_c>();
    else if (vers == v5)
      eopbend_tmpl<v5, opbend_w_d_c>();
    else if (vers == v6)
      eopbend_tmpl<v6, opbend_w_d_c>();
  } else if (opbtyp == opbend_allinger) {
    if (vers == v0 || vers == v3)
      eopbend_tmpl<v0, opbend_allinger>();
    else if (vers == v1)
      eopbend_tmpl<v1, opbend_allinger>();
    else if (vers == v4)
      eopbend_tmpl<v4, opbend_allinger>();
    else if (vers == v5)
      eopbend_tmpl<v5, opbend_allinger>();
    else if (vers == v6)
      eopbend_tmpl<v6, opbend_allinger>();
  } else {
    assert(false);
  }
}
}
TINKER_NAMESPACE_END
