#include "add.h"
#include "e_angle.h"
#include "e_opbend.h"
#include "mathfunc_const.h"
#include "md.h"
#include <cassert>

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
template <int USE, eopbend_t TYP>
void eopbend_tmpl()
{
   constexpr int do_e = USE & calc::energy;
   constexpr int do_g = USE & calc::grad;
   constexpr int do_v = USE & calc::virial;
   sanity_check<USE>();

   auto bufsize = buffer_size();

   #pragma acc parallel loop independent async\
               deviceptr(x,y,z,gx,gy,gz,\
               iopb,opbk,iang,\
               eopb,vir_eopb)
   for (int iopbend = 0; iopbend < nopbend; ++iopbend) {
      int offset = iopbend & (bufsize - 1);
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
      if CONSTEXPR (TYP == eopbend_t::w_d_c) {

         // W-D-C angle between A-B-C plane and B-D vector for D-B<AC

         rab2 = xab * xab + yab * yab + zab * zab;
         rcb2 = xcb * xcb + ycb * ycb + zcb * zcb;
         cc = rab2 * rcb2 -
            (xab * xcb + yab * ycb + zab * zcb) *
               (xab * xcb + yab * ycb + zab * zcb);
         if CONSTEXPR (do_g)
            dot = xab * xcb + yab * ycb + zab * zcb;
      } else if CONSTEXPR (TYP == eopbend_t::allinger) {

         // Allinger angle between A-C-D plane and D-B vector for D-B<AC

         rad2 = xad * xad + yad * yad + zad * zad;
         rcd2 = xcd * xcd + ycd * ycd + zcd * zcd;
         cc = rad2 * rcd2 -
            (xad * xcd + yad * ycd + zad * zcd) *
               (xad * xcd + yad * ycd + zad * zcd);
         if CONSTEXPR (do_g)
            dot = xad * xcd + yad * ycd + zad * zcd;
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

         if CONSTEXPR (do_e) {
            real e = opbunit * force * dt2 *
               (1 + copb * dt + qopb * dt2 + popb * dt3 + sopb * dt4);
            atomic_add(e, eopb, offset);
         }

         if CONSTEXPR (do_g) {
            real deddt = opbunit * force * dt * radian *
               (2 + 3 * copb * dt + 4 * qopb * dt2 + 5 * popb * dt3 +
                6 * sopb * dt4);
            real dedcos =
               -deddt * REAL_SIGN(1, ee) * REAL_RSQRT(cc * rdb2 - ee * ee);
            real term = ee * REAL_RECIP(cc);
            real dccdxia, dccdyia, dccdzia;
            real dccdxic, dccdyic, dccdzic;
            real dccdxid, dccdyid, dccdzid;
            if CONSTEXPR (TYP == eopbend_t::w_d_c) {
               dccdxia = (xab * rcb2 - xcb * dot) * term;
               dccdyia = (yab * rcb2 - ycb * dot) * term;
               dccdzia = (zab * rcb2 - zcb * dot) * term;
               dccdxic = (xcb * rab2 - xab * dot) * term;
               dccdyic = (ycb * rab2 - yab * dot) * term;
               dccdzic = (zcb * rab2 - zab * dot) * term;
               dccdxid = 0;
               dccdyid = 0;
               dccdzid = 0;
            } else if CONSTEXPR (TYP == eopbend_t::allinger) {
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

            atomic_add(dedxia, gx, ia);
            atomic_add(dedyia, gy, ia);
            atomic_add(dedzia, gz, ia);
            atomic_add(dedxib, gx, ib);
            atomic_add(dedyib, gy, ib);
            atomic_add(dedzib, gz, ib);
            atomic_add(dedxic, gx, ic);
            atomic_add(dedyic, gy, ic);
            atomic_add(dedzic, gz, ic);
            atomic_add(dedxid, gx, id);
            atomic_add(dedyid, gy, id);
            atomic_add(dedzid, gz, id);

            if CONSTEXPR (do_v) {
               real vxx = xab * dedxia + xcb * dedxic + xdb * dedxid;
               real vyx = yab * dedxia + ycb * dedxic + ydb * dedxid;
               real vzx = zab * dedxia + zcb * dedxic + zdb * dedxid;
               real vyy = yab * dedyia + ycb * dedyic + ydb * dedyid;
               real vzy = zab * dedyia + zcb * dedyic + zdb * dedyid;
               real vzz = zab * dedzia + zcb * dedzic + zdb * dedzid;

               atomic_add(vxx, vyx, vzx, vyy, vzy, vzz, vir_eopb, offset);
            }
         }
      }
   } // end for (int iopbend)
}

void eopbend_acc(int vers)
{
   if (opbtyp == eopbend_t::w_d_c) {
      if (vers == calc::v0 || vers == calc::v3)
         eopbend_tmpl<calc::v0, eopbend_t::w_d_c>();
      else if (vers == calc::v1)
         eopbend_tmpl<calc::v1, eopbend_t::w_d_c>();
      else if (vers == calc::v4)
         eopbend_tmpl<calc::v4, eopbend_t::w_d_c>();
      else if (vers == calc::v5)
         eopbend_tmpl<calc::v5, eopbend_t::w_d_c>();
      else if (vers == calc::v6)
         eopbend_tmpl<calc::v6, eopbend_t::w_d_c>();
   } else if (opbtyp == eopbend_t::allinger) {
      if (vers == calc::v0 || vers == calc::v3)
         eopbend_tmpl<calc::v0, eopbend_t::allinger>();
      else if (vers == calc::v1)
         eopbend_tmpl<calc::v1, eopbend_t::allinger>();
      else if (vers == calc::v4)
         eopbend_tmpl<calc::v4, eopbend_t::allinger>();
      else if (vers == calc::v5)
         eopbend_tmpl<calc::v5, eopbend_t::allinger>();
      else if (vers == calc::v6)
         eopbend_tmpl<calc::v6, eopbend_t::allinger>();
   } else {
      assert(false);
   }
}
TINKER_NAMESPACE_END
