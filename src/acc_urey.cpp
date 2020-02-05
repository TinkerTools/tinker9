#include "add.h"
#include "e_urey.h"
#include "md.h"
#include "named_struct.h"

TINKER_NAMESPACE_BEGIN
template <class Ver>
void eurey_acc1()
{
   constexpr int do_e = Ver::e;
   constexpr int do_g = Ver::g;
   constexpr int do_v = Ver::v;

   auto bufsize = buffer_size();

   #pragma acc parallel loop independent async\
               deviceptr(x,y,z,gx,gy,gz,\
               iury,uk,ul,\
               eub,vir_eub)
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

      if CONSTEXPR (do_e) {
         real e = ureyunit * force * dt2 * (1 + cury * dt + qury * dt2);
         atomic_add(e, eub, offset);
      }

      if CONSTEXPR (do_g) {
         real deddt =
            2 * ureyunit * force * dt * (1 + 1.5f * cury * dt + 2 * qury * dt2);
         real de = deddt * REAL_RECIP(rac);
         real dedx = de * xac;
         real dedy = de * yac;
         real dedz = de * zac;

         atomic_add(dedx, gx, ia);
         atomic_add(dedy, gy, ia);
         atomic_add(dedz, gz, ia);
         atomic_add(-dedx, gx, ic);
         atomic_add(-dedy, gy, ic);
         atomic_add(-dedz, gz, ic);

         if CONSTEXPR (do_v) {
            real vxx = xac * dedx;
            real vyx = yac * dedx;
            real vzx = zac * dedx;
            real vyy = yac * dedy;
            real vzy = zac * dedy;
            real vzz = zac * dedz;

            atomic_add(vxx, vyx, vzx, vyy, vzy, vzz, vir_eub, offset);
         }
      }
   } // end for (int i)
}

void eurey_acc(int vers)
{
   if (vers == calc::v0 || vers == calc::v3)
      eurey_acc1<EnergyVersion0>();
   else if (vers == calc::v1)
      eurey_acc1<EnergyVersion1>();
   else if (vers == calc::v4)
      eurey_acc1<EnergyVersion4>();
   else if (vers == calc::v5)
      eurey_acc1<EnergyVersion5>();
   else if (vers == calc::v6)
      eurey_acc1<EnergyVersion6>();
}
TINKER_NAMESPACE_END
