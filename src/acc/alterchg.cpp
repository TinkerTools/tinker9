#include "add.h"
#include "ff/amoeba/elecamoeba.h"
#include "ff/elec.h"
#include "ff/hippo/cflux.h"
#include "ff/hippo/elechippo.h"
#include "ff/pchg/evalence.h"
#include "math/const.h"
#include "math/libfunc.h"
#include "md/inc.h"
#include <cassert>

namespace tinker {
void bndchg_acc1()
{
   #pragma acc parallel loop independent async\
               deviceptr(x,y,z,ibnd,bflx,bl,\
               pdelta)
   for (int i = 0; i < nbond; ++i) {
      int ia = ibnd[i][0];
      int ib = ibnd[i][1];

      real pb = bflx[i];
      real ideal = bl[i];
      real xab = x[ia] - x[ib];
      real yab = y[ia] - y[ib];
      real zab = z[ia] - z[ib];
      real rab = REAL_SQRT(xab * xab + yab * yab + zab * zab);
      real dq = pb * (rab - ideal);

      atomic_add(-dq, pdelta, ia);
      atomic_add(dq, pdelta, ib);
   } // end for (int i)
}

void angchg_acc1()
{
   #pragma acc parallel loop independent async\
               deviceptr(x,y,z,iang,anat,angtyp,\
               aflx,abflx,bl,balist,pdelta)
   for (int i = 0; i < nangle; ++i) {
      int ia = iang[i][0];
      int ib = iang[i][1];
      int ic = iang[i][2];

      real pa1 = aflx[i][0];
      real pa2 = aflx[i][1];
      real pb1 = abflx[i][0];
      real pb2 = abflx[i][1];
      real ideal = anat[i];
      real xia = x[ia];
      real yia = y[ia];
      real zia = z[ia];
      real xib = x[ib];
      real yib = y[ib];
      real zib = z[ib];
      real xic = x[ic];
      real yic = y[ic];
      real zic = z[ic];
      real xab = xia - xib;
      real yab = yia - yib;
      real zab = zia - zib;
      real xcb = xic - xib;
      real ycb = yic - yib;
      real zcb = zic - zib;
      real rab = REAL_SQRT(xab * xab + yab * yab + zab * zab);
      real rcb = REAL_SQRT(xcb * xcb + ycb * ycb + zcb * zcb);

      real dot, angle, cosine;
      if (rab != 0 and rcb != 0) {
         dot = xab * xcb + yab * ycb + zab * zcb;
         cosine = dot / (rab * rcb);
         cosine = REAL_MIN(1.0, REAL_MAX(-1, cosine));
         angle = radian * REAL_ACOS(cosine);
      }

      int ab = balist[i][0];
      int cb = balist[i][1];
      real rab0 = bl[ab];
      real rcb0 = bl[cb];
      real dq1 = pb1 * (rcb - rcb0) + pa1 * (angle - ideal) / radian;
      real dq2 = pb2 * (rab - rab0) + pa2 * (angle - ideal) / radian;
      atomic_add(dq1, pdelta, ia);
      atomic_add(dq2, pdelta, ic);
      atomic_add(-(dq1 + dq2), pdelta, ib);
   } // end for (int i)
}

void alterchg_acc()
{
   darray::zero(g::q0, n, pdelta);

   bndchg_acc1();
   angchg_acc1();

   // alter monopoles and charge penetration
   #pragma acc parallel loop independent async\
           deviceptr(pval,pval0,pdelta,pole,mono0)
   for (int i = 0; i < n; ++i) {
      pval[i] = pval0[i] + pdelta[i];
      pole[i][0] = mono0[i] + pdelta[i];
   }
}
}
