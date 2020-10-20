#include "add.h"
#include "cflux.h"
#include "eangle.h"
#include "ebond.h"
#include "elec.h"
#include "mathfunc_const.h"
#include "md.h"
#include <cassert>
#include <tinker/detail/atmlst.hh>
#include <tinker/detail/atomid.hh>
#include <tinker/detail/mplpot.hh>


namespace tinker {
void bndchg_acc1()
{
   #pragma acc parallel loop independent async\
               deviceptr(x,y,z,ibnd,bflx,bl,\
               atomic,couple_n12,couple_i12,pdelta)
   for (int i = 0; i < nbond; ++i) {
      int ia = ibnd[i][0];
      int ib = ibnd[i][1];
      real pb = bflx[i];
      int atoma = atomic[ia];
      int atomb = atomic[ib];



      real ideal = bl[i];
      real xab = x[ia] - x[ib];
      real yab = y[ia] - y[ib];
      real zab = z[ia] - z[ib];

      real rab = REAL_SQRT(xab * xab + yab * yab + zab * zab);
      real dq = pb * (rab - ideal);

      // determine higher priority of the bonded atoms
      int priority;
      if (atoma > atomb)
         priority = 1;
      else if (atoma < atomb)
         priority = -1;
      else {
         int n12a = couple_n12[ia];
         int n12b = couple_n12[ib];

         if (n12a > n12b)
            priority = 1;
         else if (n12a < n12b)
            priority = -1;
         else {
            int nha = 0;
            int nhb = 0;

            for (int j = 0; j < n12a; ++j) {
               if (atomic[couple_i12[ia][j]] == 1)
                  nha += 1;
            }

            for (int j = 0; j < n12b; ++j) {
               if (atomic[couple_i12[ib][j]] == 1)
                  nhb += 1;
            }

            if (nha > nhb)
               priority = 1;
            else if (nha < nhb)
               priority = -1;
            else
               priority = 0;
         }
      }

      atomic_add(-dq * priority, pdelta, ia);
      atomic_add(dq * priority, pdelta, ib);
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
      // pdelta[ia] += dq1;
      // pdelta[ib] -= (dq1 + dq2);
      // pdelta[ic] += dq2;

      atomic_add(dq1, pdelta, ia);
      atomic_add(dq2, pdelta, ic);
      atomic_add(-(dq1+dq2), pdelta, ib);
   } // end for (int i)
}

void alterchg_acc()
{
   darray::zero(PROCEED_NEW_Q, n, pdelta);

   bndchg_acc1();
   angchg_acc1();
   
   // alter atomic partial charge values
   // nion, iion

   // alter monopoles and charge penetration

   // for (int i = 0; i < npole; ++i) {

   // charge penenetration
   //if (mplpot::use_chgpen)
   #pragma acc parallel loop independent async\
            deviceptr(pval,pval0,pdelta,pole,mono0)
   for (int i = 0; i < n; ++i) {
      pval[i] = pval0[i] + pdelta[i];
      pole[i][0] = mono0[i] + pdelta[i];
   }

   // pole[i][0] = mono0[i] + pdelta[k];
   // for (int i = 0; i < nion; ++i) {
   //    k = iion[i];
   //    pchg[i] = pchg0[i] + pdelta[k]
   // }
}
}
