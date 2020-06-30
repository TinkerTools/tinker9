#include "add.h"
#include "eimptor.h"
#include "md.h"


namespace tinker {
template <class Ver>
void eimptor_acc1()
{
   constexpr bool do_e = Ver::e;
   constexpr bool do_g = Ver::g;
   constexpr bool do_v = Ver::v;


   size_t bufsize = buffer_size();


   #pragma acc parallel loop independent async\
               deviceptr(x,y,z,deitx,deity,deitz,\
               iitors,itors1,itors2,itors3,\
               eit,vir_eit)
   for (int i = 0; i < nitors; ++i) {
      int offset = i & (bufsize - 1);
      const int ia = iitors[i][0];
      const int ib = iitors[i][1];
      const int ic = iitors[i][2];
      const int id = iitors[i][3];


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
      real xba = xib - xia;
      real yba = yib - yia;
      real zba = zib - zia;
      real xcb = xic - xib;
      real ycb = yic - yib;
      real zcb = zic - zib;
      real xdc = xid - xic;
      real ydc = yid - yic;
      real zdc = zid - zic;


      real xt = yba * zcb - ycb * zba;
      real yt = zba * xcb - zcb * xba;
      real zt = xba * ycb - xcb * yba;
      real xu = ycb * zdc - ydc * zcb;
      real yu = zcb * xdc - zdc * xcb;
      real zu = xcb * ydc - xdc * ycb;
      real xtu = yt * zu - yu * zt;
      real ytu = zt * xu - zu * xt;
      real ztu = xt * yu - xu * yt;
      real rt2 = xt * xt + yt * yt + zt * zt;
      real ru2 = xu * xu + yu * yu + zu * zu;
      real rtru = REAL_SQRT(rt2 * ru2);


      if (rtru != 0) {
         real rcb = REAL_SQRT(xcb * xcb + ycb * ycb + zcb * zcb);
         real cosine = (xt * xu + yt * yu + zt * zu) * REAL_RECIP(rtru);
         real sine =
            (xcb * xtu + ycb * ytu + zcb * ztu) * REAL_RECIP(rcb * rtru);


         real v1 = itors1[i][0];
         real c1 = itors1[i][2];
         real s1 = itors1[i][3];
         real v2 = itors2[i][0];
         real c2 = itors2[i][2];
         real s2 = itors2[i][3];
         real v3 = itors3[i][0];
         real c3 = itors3[i][2];
         real s3 = itors3[i][3];


         real cosine2 = cosine * cosine - sine * sine;
         real sine2 = 2 * cosine * sine;
         real cosine3 = cosine * cosine2 - sine * sine2;
         real sine3 = cosine * sine2 + sine * cosine2;
         real phi1 = 1 + (cosine * c1 + sine * s1);
         real phi2 = 1 + (cosine2 * c2 + sine2 * s2);
         real phi3 = 1 + (cosine3 * c3 + sine3 * s3);


         if CONSTEXPR (do_e) {
            real e = itorunit * (v1 * phi1 + v2 * phi2 + v3 * phi3);
            atomic_add(e, eit, offset);
         }


         if CONSTEXPR (do_g) {
            real dphi1 = (cosine * s1 - sine * c1);
            real dphi2 = 2 * (cosine2 * s2 - sine2 * c2);
            real dphi3 = 3 * (cosine3 * s3 - sine3 * c3);
            real dedphi = itorunit * (v1 * dphi1 + v2 * dphi2 + v3 * dphi3);


            real xca = xic - xia;
            real yca = yic - yia;
            real zca = zic - zia;
            real xdb = xid - xib;
            real ydb = yid - yib;
            real zdb = zid - zib;

            real inv1 = REAL_RECIP(rt2 * rcb);
            real inv2 = REAL_RECIP(ru2 * rcb);
            real dedxt = dedphi * (yt * zcb - ycb * zt) * inv1;
            real dedyt = dedphi * (zt * xcb - zcb * xt) * inv1;
            real dedzt = dedphi * (xt * ycb - xcb * yt) * inv1;
            real dedxu = -dedphi * (yu * zcb - ycb * zu) * inv2;
            real dedyu = -dedphi * (zu * xcb - zcb * xu) * inv2;
            real dedzu = -dedphi * (xu * ycb - xcb * yu) * inv2;


            real dedxia = zcb * dedyt - ycb * dedzt;
            real dedyia = xcb * dedzt - zcb * dedxt;
            real dedzia = ycb * dedxt - xcb * dedyt;
            real dedxib = yca * dedzt - zca * dedyt + zdc * dedyu - ydc * dedzu;
            real dedyib = zca * dedxt - xca * dedzt + xdc * dedzu - zdc * dedxu;
            real dedzib = xca * dedyt - yca * dedxt + ydc * dedxu - xdc * dedyu;
            real dedxic = zba * dedyt - yba * dedzt + ydb * dedzu - zdb * dedyu;
            real dedyic = xba * dedzt - zba * dedxt + zdb * dedxu - xdb * dedzu;
            real dedzic = yba * dedxt - xba * dedyt + xdb * dedyu - ydb * dedxu;
            real dedxid = zcb * dedyu - ycb * dedzu;
            real dedyid = xcb * dedzu - zcb * dedxu;
            real dedzid = ycb * dedxu - xcb * dedyu;


            atomic_add(dedxia, deitx, ia);
            atomic_add(dedyia, deity, ia);
            atomic_add(dedzia, deitz, ia);
            atomic_add(dedxib, deitx, ib);
            atomic_add(dedyib, deity, ib);
            atomic_add(dedzib, deitz, ib);
            atomic_add(dedxic, deitx, ic);
            atomic_add(dedyic, deity, ic);
            atomic_add(dedzic, deitz, ic);
            atomic_add(dedxid, deitx, id);
            atomic_add(dedyid, deity, id);
            atomic_add(dedzid, deitz, id);


            if CONSTEXPR (do_v) {
               real vxx = xcb * (dedxic + dedxid) - xba * dedxia + xdc * dedxid;
               real vyx = ycb * (dedxic + dedxid) - yba * dedxia + ydc * dedxid;
               real vzx = zcb * (dedxic + dedxid) - zba * dedxia + zdc * dedxid;
               real vyy = ycb * (dedyic + dedyid) - yba * dedyia + ydc * dedyid;
               real vzy = zcb * (dedyic + dedyid) - zba * dedyia + zdc * dedyid;
               real vzz = zcb * (dedzic + dedzid) - zba * dedzia + zdc * dedzid;
               atomic_add(vxx, vyx, vzx, vyy, vzy, vzz, vir_eit, offset);
            }
         }
      }
   }
}


void eimptor_acc(int vers)
{
   if (vers == calc::v0 || vers == calc::v3)
      eimptor_acc1<calc::V0>();
   else if (vers == calc::v1)
      eimptor_acc1<calc::V1>();
   else if (vers == calc::v4)
      eimptor_acc1<calc::V4>();
   else if (vers == calc::v5)
      eimptor_acc1<calc::V5>();
   else if (vers == calc::v6)
      eimptor_acc1<calc::V6>();
}
}
