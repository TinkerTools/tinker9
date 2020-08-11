#include "add.h"
#include "ebond.h"
#include "md.h"
#include <cassert>

namespace tinker {
template <class Ver, class BNDTYP>
void bndchg_acc1()
{
   auto bufsize = buffer_size();

   #pragma acc parallel loop independent async\
               deviceptr(x,y,z, atomic, i12, \
               ibnd,bflux,bl,\
               pdelta)
   for (int i = 0; i < nbond; ++i) {
      int offset = i & (bufsize - 1);
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
         int n12a = n12[ia];
         int n12b = n12[ib];

         if (n12a > n12b)
            priority = 1;
         else if (n12a < n12b)
            priority = -1;
         else {
            int nha = 0;
            int nhb = 0;

            for (int j = 0; j < n12a; ++j) {
                if (atomic[i12[j][ia]] == 1)
                  nha += 1;
            }

            for (int j = 0; j < n12b; ++j) {
                if (atomic[i12[j][ib]] == 1)
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

      pdelta[ia] -= dq * priority;
      pdelta[ib] += dq * priority;
   } // end for (int i)
}


template <class Ver>
void angchg_acc1()
{
   auto bufsize = buffer_size();

   #pragma acc parallel loop independent async\
               deviceptr(x,y,z,\
               iang,anat,balist,angtyp, \
               aflx, abflx,)
   for (int i = 0; i < nangle; ++i) {
      int offset = i & (bufsize - 1);
      int ia = iang[i][0];
      int ib = iang[i][1];
      int ic = iang[i][2];
      
      real pa1 = aflx[i][0];
      real pa1 = aflx[i][1];
      real pb1 = abflx[i][0];
      real pb2 = abflx[i][1];

      real ideal = anat[i];
      real force = ak[i];

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

      real rab2 = xab * xab + yab * yab + zab * zab;
      real rcb2 = xcb * xcb + ycb * ycb + zcb * zcb;

      if (rab2 != 0 && rcb2 != 0) {
         real dot = xab * xcb + yab * ycb + zab * zcb;
         real cosine = dot * REAL_RSQRT(rab2 * rcb2);
         cosine = REAL_MIN((real)0.1, REAL_MAX(-1, cosine));
         real angle = radian * REAL_ACOS(cosine);
      }

      real rab0 = bl[balist[i][0]];
      real rcb0 = bl[balist[i][1]];

      dq1 = pb1 * (rcb - rcb0) + pa1 * (angle - ideal) / radian;
      dq2 = pb2 * (rab - rab0) + pa2 * (angle - ideal) / radian;
      pdelta[ia] += dq1;
      pdelta[ib] -= (dq1 + dq2);
      pdelta[ic] += dq2;
   } // end for (int i)
}

void alterchg() {
   real pdelta[n];
   for (int i = 0; i < n; ++i)
      pdelta[i] = 0;

   bndchg_acc1();
   angchg_acc1();
   // alter atomic partial charge values
   // nion, iion 

   //alter monopoles and charge penetration

   for (int i = 0; i < npole; ++i) {
      int k = ipole[i];
      pole[i][0] = mono0[i] + pdelta[k];

      // charge penenetration
      if CONSTEXPR (use_chgpen)
         pval[i] = pval0[i] + pdelta[k];
   }

   
   for (int i = 0; i < nion; ++i) {
      k = iion[i];
      pchg[i] = pchg0[i] + pdelta[k]
   }
   
}
}
